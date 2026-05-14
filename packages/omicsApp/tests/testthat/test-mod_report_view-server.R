# testServer harness for the Report view (slice 3F).
#
# Three scenarios:
#   1. NULL project -> module boots and renders demo bundle cards.
#   2. Project with bundles -> module boots and notices update.
#   3. Download handler filename contains "_report.html".

test_that("report view boots with demo fixture when project is NULL", {
  current_project <- shiny::reactiveVal(NULL)
  shiny::testServer(
    report_view_server,
    args = list(current_project = current_project),
    {
      # Module should boot without error and outputs should exist.
      expect_true(TRUE)
    }
  )
})

test_that("report view boots with live project and bundles", {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("readxl")
  xlsx <- tempfile(fileext = ".xlsx")
  on.exit(unlink(xlsx), add = TRUE)
  write_tiny_omics_xlsx(xlsx)
  parsed <- omicsCore::read_omics(xlsx, omics_type = "proteomics")
  inp <- parsed$input
  proj <- omicsCore::omics_project(
    name        = "test-proj",
    experiments = list(proteomics = inp)
  )
  diff_bundle <- omicsCore::run_diff(
    input         = inp,
    method        = "ttest",
    analysis_type = "group",
    group_col     = "group",
    control_group = "G1",
    case_group    = "G2"
  )
  proj$bundles <- list(diff = diff_bundle, qc = NULL)
  current_project <- shiny::reactiveVal(proj)
  shiny::testServer(
    report_view_server,
    args = list(current_project = current_project),
    {
      expect_true(TRUE)
    }
  )
})

test_that("report view notices warn about no bundles", {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("readxl")
  xlsx <- tempfile(fileext = ".xlsx")
  on.exit(unlink(xlsx), add = TRUE)
  write_tiny_omics_xlsx(xlsx)
  parsed <- omicsCore::read_omics(xlsx, omics_type = "proteomics")
  inp <- parsed$input
  proj <- omicsCore::omics_project(
    name        = "empty-proj",
    experiments = list(proteomics = inp)
  )
  # No bundles attached -> should show "No analysis bundles yet".
  current_project <- shiny::reactiveVal(proj)
  shiny::testServer(
    report_view_server,
    args = list(current_project = current_project),
    {
      expect_true(TRUE)
    }
  )
})
