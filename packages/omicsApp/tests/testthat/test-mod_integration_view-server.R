# testServer harness for the Integration view (slice 3E).
#
# Three scenarios:
#   1. NULL project + NULL diff_bundle → demo fixture; can_run() is FALSE.
#   2. Project with one experiment + a diff_bundle → can_run() is FALSE
#      (needs >=2 experiments) and the demo fixture is used.
#   3. Project with two experiments where the secondary's meta_df has
#      the same group_col / control / case levels → can_run() is TRUE
#      and run_integration(concordance) populates `integration_bundle()`.

test_that("integration view falls back to demo when prerequisites are absent", {
  current_project <- shiny::reactiveVal(NULL)
  diff_bundle <- shiny::reactiveVal(NULL)
  shiny::testServer(
    integration_view_server,
    args = list(current_project = current_project,
                diff_bundle     = diff_bundle),
    {
      session$setInputs(rerun = 0)
      expect_true(isTRUE(is_demo()))
      info <- can_run()
      expect_false(isTRUE(info$ok))
      expect_null(integration_bundle())
    }
  )
})

test_that("integration view stays on demo when project has only one experiment", {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("readxl")
  xlsx <- tempfile(fileext = ".xlsx")
  on.exit(unlink(xlsx), add = TRUE)
  write_tiny_omics_xlsx(xlsx)
  parsed <- omicsCore::read_omics(xlsx, omics_type = "proteomics")
  proj <- omicsCore::omics_project(
    name        = "single-layer",
    experiments = list(proteomics = parsed$input)
  )
  current_project <- shiny::reactiveVal(proj)
  diff_bundle <- shiny::reactiveVal(omicsApp:::example_diff_bundle())
  shiny::testServer(
    integration_view_server,
    args = list(current_project = current_project,
                diff_bundle     = diff_bundle),
    {
      session$setInputs(rerun = 0)
      expect_true(isTRUE(is_demo()))
      info <- can_run()
      expect_false(isTRUE(info$ok))
    }
  )
})

test_that("integration view runs concordance when two compatible layers exist", {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("readxl")
  xlsx <- tempfile(fileext = ".xlsx")
  on.exit(unlink(xlsx), add = TRUE)
  write_tiny_omics_xlsx(xlsx)
  parsed <- omicsCore::read_omics(xlsx, omics_type = "proteomics")
  inp_a <- parsed$input
  # Build a second compatible layer: same meta_df, fresh expression.
  inp_b <- inp_a
  inp_b$omics_type <- "rnaseq"
  inp_b$assay_type <- "log_expr"  # treat as log space so ttest is happy
  proj <- omicsCore::omics_project(
    name        = "dual",
    experiments = list(proteomics = inp_a, rnaseq = inp_b)
  )
  primary <- omicsCore::run_diff(
    input         = inp_a,
    method        = "ttest",
    analysis_type = "group",
    group_col     = "group",
    control_group = "G1",
    case_group    = "G2"
  )
  current_project <- shiny::reactiveVal(proj)
  diff_bundle <- shiny::reactiveVal(primary)
  shiny::testServer(
    integration_view_server,
    args = list(current_project = current_project,
                diff_bundle     = diff_bundle),
    {
      # Touch an input so the auto-run observer fires.
      session$setInputs(rerun = 0)
      info <- can_run()
      expect_true(isTRUE(info$ok))
      session$setInputs(rerun = 1)
      b <- integration_bundle()
      expect_s3_class(b, "analysis_bundle")
      expect_identical(b$analysis_name, "run_integration")
      expect_false(isTRUE(is_demo()))
    }
  )
})
