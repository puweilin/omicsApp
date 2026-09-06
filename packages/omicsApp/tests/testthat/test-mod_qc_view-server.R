# testServer harness for the QC view (slice 3C).
#
# Drives the module with two scenarios:
#   1. No project (current_project = NULL) → falls back to
#      example_qc_bundle() and renders the demo header.
#   2. Live project built from the tiny synthetic xlsx →
#      re-runs run_qc() against the user's experiment and
#      re-derives the bundle when the missingness slider moves.
#
# The QC view has no Run button (per slice-3 convention QC is
# cheap); a slider change must propagate to `last_bundle`.

test_that("qc view falls back to example_qc_bundle when project is NULL", {
  current_project <- shiny::reactiveVal(NULL)
  shiny::testServer(
    qc_view_server,
    args = list(current_project = current_project),
    {
      # Drive the controls to their defaults.
      session$setInputs(missing_threshold = 0.5, outlier_method = "iqr")
      bundle <- last_bundle()
      expect_s3_class(bundle, "analysis_bundle")
      expect_identical(bundle$analysis_name, "run_qc")
      # Fixture is proteomics, 80 features × 24 samples.
      expect_equal(bundle$input_info$omics_type, "proteomics")
      expect_null(last_error())
    }
  )
})

test_that("qc view re-runs run_qc against the live project", {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("readxl")

  xlsx <- tempfile(fileext = ".xlsx")
  on.exit(unlink(xlsx), add = TRUE)
  write_tiny_omics_xlsx(xlsx)
  parsed <- omicsCore::read_omics(xlsx, omics_type = "proteomics", assay_type = "normalized_intensity")
  inp <- parsed$input
  proj <- omicsCore::omics_project(
    name        = "test",
    experiments = list(proteomics = inp)
  )
  current_project <- shiny::reactiveVal(proj)

  shiny::testServer(
    qc_view_server,
    args = list(current_project = current_project),
    {
      session$setInputs(missing_threshold = 0.5, outlier_method = "iqr")
      bundle <- last_bundle()
      expect_s3_class(bundle, "analysis_bundle")
      # 5 features × 6 samples (no NAs in the fixture, no outliers
      # at the loose default threshold).
      expect_equal(bundle$input_info$n_features_in, 5L)
      expect_equal(bundle$input_info$n_samples_in,  6L)
      expect_equal(bundle$input_info$n_features_out, 5L)

      # Tighten the missing threshold to 0 — features with any NA
      # would be flagged; the fixture has none, so the count stays
      # at 5 but the param round-trips to the bundle.
      session$setInputs(missing_threshold = 0.0)
      expect_equal(last_bundle()$params$missing_threshold, 0.0)
    }
  )
})

test_that("qc view surfaces run_qc errors instead of crashing", {
  # Build a degenerate input that run_qc will reject (a single
  # sample after subsetting). We achieve this by feeding a real
  # input but cranking the sample_missing_threshold path via a
  # custom run_qc shim is overkill; the simpler route is to mock
  # active() by handing in a project whose experiment has just one
  # column — run_qc will refuse outlier detection on it.
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("readxl")

  xlsx <- tempfile(fileext = ".xlsx")
  on.exit(unlink(xlsx), add = TRUE)
  write_tiny_omics_xlsx(xlsx, n_samples = 6L)
  parsed <- omicsCore::read_omics(xlsx, omics_type = "proteomics", assay_type = "normalized_intensity")
  inp <- parsed$input
  proj <- omicsCore::omics_project(
    name        = "test",
    experiments = list(proteomics = inp)
  )
  current_project <- shiny::reactiveVal(proj)

  shiny::testServer(
    qc_view_server,
    args = list(current_project = current_project),
    {
      # missing_threshold = 0 with a synthetic NA-free input keeps
      # things sane; flip it to a value that drives all features
      # out instead — run_qc stops with "QC would remove all
      # samples or features".  We inject NAs by mutating the input
      # in place is awkward in testServer; the simpler way to hit
      # the error path is to pass an unsupported outlier method
      # via setInputs (radio values are validated by the choices
      # list at the UI level, but omicsCore refuses anything else by
      # name server-side).
      session$setInputs(missing_threshold = 0.5,
                        outlier_method    = "unknown_method")
      expect_false(is.null(last_error()))
      expect_match(last_error(), "`outlier_method` must be one or more of", fixed = TRUE)
    }
  )
})
