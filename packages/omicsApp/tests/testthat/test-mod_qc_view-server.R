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

# ---- a saved result -------------------------------------------------

# The tutorial project with a QC result saved at settings no control
# opens on, as Open or Restore would bring it.
qc_saved_project <- function() {
  p <- tutorial_project()
  saved <- omicsCore::run_qc(p$experiments$proteomics, missing_threshold = 0.3,
                             outlier_method = "iqr", impute_method = "none",
                             missing_filter = "any_group", group_col = "group")
  p$bundles <- list(qc = saved)
  p
}

count_run_qc <- function(env = parent.frame()) {
  calls <- new.env()
  calls$n <- 0L
  real <- omicsCore::run_qc
  testthat::local_mocked_bindings(
    run_qc = function(...) { calls$n <- calls$n + 1L; real(...) },
    .package = "omicsCore", .env = env)
  calls
}

test_that("a saved QC result is shown with its settings, not recomputed", {
  p <- qc_saved_project()
  saved <- p$bundles$qc
  calls <- count_run_qc()
  restored <- NULL
  testthat::local_mocked_bindings(
    qc_restore_controls = function(session, vals) restored <<- vals, .package = "omicsApp")
  shiny::testServer(qc_view_server, args = list(current_project = shiny::reactiveVal(p)), {
    session$flushReact()
    expect_identical(calls$n, 0L)
    expect_identical(last_bundle(), saved)
    # Its settings go back into the controls.
    expect_identical(restored$layer, "proteomics")
    expect_equal(restored$missing_threshold, 0.3)
    expect_identical(restored$outlier_method, "iqr")
    expect_identical(restored$impute_method, "none")
    expect_identical(restored$missing_filter, "any_group")
    expect_identical(restored$missing_group_col, "group")
    html <- paste(unlist(output$stats), collapse = " ")
    expect_match(html, "30% missing in every group", fixed = TRUE)
    expect_match(html, "IQR", fixed = TRUE)

    # The browser applying those values is not a change.
    session$setInputs(missing_threshold = 0.3, outlier_method = "iqr",
                      impute_method = "none", missing_filter = "any_group",
                      missing_group_col = "group")
    expect_identical(calls$n, 0L)
    expect_identical(last_bundle(), saved)

    # Changing one is: QC runs, with the other saved settings kept.
    session$setInputs(missing_threshold = 0.6)
    expect_identical(calls$n, 1L)
    b <- last_bundle()
    expect_equal(b$params$missing_threshold, 0.6)
    expect_identical(b$params$outlier_method, "iqr")
    expect_identical(b$params$missing_filter, "any_group")
    expect_identical(b$params$group_col, "group")
  })
})

test_that("a saved result waits for controls that have not caught up yet", {
  # The controls still hold the previous project's values for a moment
  # after a project arrives; that moment is not a request to recompute.
  p <- qc_saved_project()
  calls <- count_run_qc()
  testthat::local_mocked_bindings(qc_restore_controls = function(session, vals) NULL,
                                  .package = "omicsApp")
  proj <- shiny::reactiveVal(NULL)
  shiny::testServer(qc_view_server, args = list(current_project = proj), {
    session$setInputs(missing_threshold = 0.5, outlier_method = "all")
    n_demo <- calls$n
    proj(p)
    session$flushReact()
    expect_identical(calls$n, n_demo)
    expect_identical(last_bundle(), p$bundles$qc)
  })
})

test_that("another layer is computed; the saved layer shows the saved result again", {
  p <- qc_saved_project()
  calls <- count_run_qc()
  testthat::local_mocked_bindings(qc_restore_controls = function(session, vals) NULL,
                                  .package = "omicsApp")
  shiny::testServer(qc_view_server, args = list(current_project = shiny::reactiveVal(p)), {
    session$flushReact()
    expect_identical(calls$n, 0L)
    session$setInputs(layer = "rnaseq")
    expect_identical(calls$n, 1L)
    expect_identical(last_bundle()$input_info$omics_type, "rnaseq")
    session$setInputs(layer = "proteomics")
    expect_identical(calls$n, 1L)
    expect_identical(last_bundle(), p$bundles$qc)
  })
})

test_that("a change to the project that leaves the layer alone does not rerun QC", {
  p <- tutorial_project()
  calls <- count_run_qc()
  proj <- shiny::reactiveVal(p)
  shiny::testServer(qc_view_server, args = list(current_project = proj), {
    session$setInputs(missing_threshold = 0.5, outlier_method = "pca")
    n <- calls$n
    expect_gte(n, 1L)
    q <- proj()
    q$bundles$diff <- "another view's result"
    proj(q)
    session$flushReact()
    expect_identical(calls$n, n)
    # The view publishing its own result is not a saved result to adopt.
    q$bundles$qc <- last_bundle()
    proj(q)
    session$flushReact()
    expect_identical(calls$n, n)
  })
})

test_that("a restored session keeps its saved QC result, in the app and on disk", {
  # The whole wiring: the generation bump a restore causes, the view
  # adopting the saved result, and the bundle-attach observer and the
  # autosave seeing that result rather than a recomputed one.
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("readxl")
  skip_if_not_installed("withr")
  store <- file.path(withr::local_tempdir(), "store")
  dir.create(store)
  withr::local_envvar(OMICSAPP_DATA_DIR = store)
  xlsx <- withr::local_tempfile(fileext = ".xlsx")
  write_tiny_omics_xlsx(xlsx, n_features = 40L, n_samples = 8L, seed = 1)

  shiny::testServer(app_server, {
    suppressWarnings(session$setInputs(
      `import-omics_type` = "proteomics",
      `import-file` = list(datapath = xlsx, name = basename(xlsx),
                           size = file.info(xlsx)$size)))
    session$setInputs(`import-confirm` = 1)
    session$setInputs(`qc-missing_threshold` = 0.3, `qc-outlier_method` = "iqr")
    expect_equal(current_project()$bundles$qc$params$missing_threshold, 0.3)
  })
  expect_equal(omicsCore::load_project(autosave_file(store))$bundles$qc$params$missing_threshold, 0.3)

  calls <- count_run_qc()
  shiny::testServer(app_server, {
    session$flushReact()
    expect_identical(calls$n, 0L)
    expect_equal(qc_bundle()$params$missing_threshold, 0.3)
    expect_identical(qc_bundle()$params$outlier_method, "iqr")
    expect_equal(current_project()$bundles$qc$params$missing_threshold, 0.3)
  })
  expect_equal(omicsCore::load_project(autosave_file(store))$bundles$qc$params$missing_threshold, 0.3)
})

# ---- group-wise missing filter -----------------------------------------

test_that("the missing filter offers the recorded group column and reaches run_qc", {
  p <- tutorial_project()
  shiny::testServer(qc_view_server, args = list(current_project = shiny::reactiveVal(p)), {
    session$setInputs(layer = "proteomics", missing_threshold = 0.5, outlier_method = "pca")
    expect_identical(last_bundle()$params$missing_filter, "global")
    html <- paste(unlist(output$ui_missing_filter), collapse = " ")
    expect_match(html, "In at least one group", fixed = TRUE)
    expect_match(html, '<option value="group" selected>', fixed = TRUE)

    session$setInputs(missing_filter = "any_group")
    b <- last_bundle()
    expect_identical(b$params$missing_filter, "any_group")
    # Unset, the group column is the layer's recorded one.
    expect_identical(b$params$group_col, "group")

    session$setInputs(missing_filter = "all_groups", missing_group_col = "group")
    expect_identical(last_bundle()$params$missing_filter, "all_groups")
  })
})

test_that("the missing filter is not offered when nothing splits the samples", {
  set.seed(1)
  x <- matrix(stats::rnorm(60, 20), 10, dimnames = list(paste0("f", 1:10), paste0("s", 1:6)))
  inp <- omicsCore::omics_input(x, data.frame(sample = colnames(x), row.names = colnames(x)),
                                data.frame(feature_id = rownames(x), row.names = rownames(x)),
                                omics_type = "proteomics", assay_type = "normalized_intensity")
  p <- omicsCore::omics_project("flat", experiments = list(proteomics = inp))
  shiny::testServer(qc_view_server, args = list(current_project = shiny::reactiveVal(p)), {
    # A value left over from another layer is not used here.
    session$setInputs(missing_threshold = 0.5, outlier_method = "pca",
                      missing_filter = "any_group")
    expect_null(output$ui_missing_filter$html %||% NULL)
    expect_identical(last_bundle()$params$missing_filter, "global")
    expect_null(last_error())
  })
})

# ---- leave-one-out ----------------------------------------------------

test_that("a broken sample in a small layer is flagged and explained in plain words", {
  p <- tutorial_project()
  inp <- p$experiments$proteomics
  set.seed(4)
  inp$expr_mat[, 1] <- inp$expr_mat[, 1] +
    sample(c(-1.5, 1.5), nrow(inp$expr_mat), replace = TRUE)
  p$experiments <- list(proteomics = inp)
  s1 <- colnames(inp$expr_mat)[1]
  shiny::testServer(qc_view_server, args = list(current_project = shiny::reactiveVal(p)), {
    session$setInputs(missing_threshold = 0.5, outlier_method = "all", impute_method = "none")
    b <- last_bundle()
    expect_identical(b$params$outlier_method, c("pca", "connectivity", "iqr", "loo"))
    expect_true(s1 %in% b$results$qc_summary$outliers$by_method$loo$flagged_samples)
    html <- paste(unlist(output$notices), collapse = " ")
    expect_match(html, sprintf("Leave-one-out: %s correlates", s1), fixed = TRUE)
    if (ncol(inp$expr_mat) <= 10L) {
      expect_match(html, "The leave-one-out check can", fixed = TRUE)
      expect_false(grepl("z-score", html, fixed = TRUE))
    }
  })
})

# ---- the cleaning record ------------------------------------------------
# run_qc() no longer stores a copy of the cleaned input; the view reads the
# record instead, and results saved with the copy keep working.

test_that("the view reads kept samples, colours and missing cells from the record", {
  p <- tutorial_project()
  inp <- p$experiments$proteomics
  inp$expr_mat[1:3, 1:2] <- NA
  p$experiments <- list(proteomics = inp)
  shiny::testServer(qc_view_server, args = list(current_project = shiny::reactiveVal(p)), {
    session$setInputs(missing_threshold = 0.5, outlier_method = "pca", impute_method = "min")
    b <- last_bundle()
    expect_null(b$results$cleaned_input)
    expect_identical(qc_kept_samples(b), colnames(inp$expr_mat))
    expect_identical(qc_bundle_layer(current_project()$experiments, b), "proteomics")
    expect_true("group" %in% pca_color_choices())
    expect_identical(qc_missing_after(b), 0L)
    caption <- paste(unlist(output$missing_caption), collapse = " ")
    expect_match(caption, "after imputation for this view: 0.0%", fixed = TRUE)
    expect_s3_class(omicsCore::plot_qc(b, view = "pca", color_by = "group"), "ggplot")
  })
})

test_that("a QC result saved with its cleaned input still restores and draws", {
  p <- tutorial_project()
  inp <- p$experiments$proteomics
  b <- omicsCore::run_qc(inp, missing_threshold = 0.3, outlier_method = "iqr",
                         impute_method = "none")
  old <- b
  old$results <- list(qc_summary = b$results$qc_summary,
                      cleaned_input = omicsCore::qc_cleaned_input(b, inp))
  old$input_info$layer <- "proteomics"
  p$bundles <- list(qc = old)
  calls <- count_run_qc()
  testthat::local_mocked_bindings(qc_restore_controls = function(session, vals) NULL,
                                  .package = "omicsApp")
  shiny::testServer(qc_view_server, args = list(current_project = shiny::reactiveVal(p)), {
    session$flushReact()
    expect_identical(calls$n, 0L)
    expect_identical(last_bundle(), old)
    expect_identical(qc_kept_samples(old), colnames(inp$expr_mat))
    expect_true("group" %in% pca_color_choices())
    expect_match(paste(unlist(output$stats), collapse = " "), "Samples kept", fixed = TRUE)
    for (v in c("pca", "missing", "connectivity")) {
      expect_s3_class(omicsCore::plot_qc(old, view = v), "ggplot")
    }
  })
})

test_that("MinProb that cannot be estimated is said in plain words", {
  testthat::skip_if_not_installed("imputeLCMD")
  set.seed(5)
  x <- matrix(stats::rnorm(36, 20), 6, dimnames = list(paste0("p", 1:6), paste0("s", 1:6)))
  x[2:6, 1:3] <- NA
  inp <- omicsCore::omics_input(x, data.frame(group = rep(c("a", "b"), 3), row.names = colnames(x)),
                                data.frame(feature_id = rownames(x), row.names = rownames(x)),
                                omics_type = "proteomics", assay_type = "normalized_intensity")
  p <- omicsCore::omics_project("few", experiments = list(proteomics = inp))
  shiny::testServer(qc_view_server, args = list(current_project = shiny::reactiveVal(p)), {
    session$setInputs(missing_threshold = 1, outlier_method = "pca", impute_method = "MinProb")
    expect_null(last_error())
    b <- last_bundle()
    expect_identical(b$results$qc_summary$imputation$method, "MinDet")
    stats_html <- paste(unlist(output$stats), collapse = " ")
    expect_match(stats_html, "MinDet", fixed = TRUE)
    expect_match(stats_html, "in place of MinProb", fixed = TRUE)
    notes <- paste(unlist(output$notices), collapse = " ")
    expect_match(notes, "MinProb could not be used", fixed = TRUE)
  })
})
