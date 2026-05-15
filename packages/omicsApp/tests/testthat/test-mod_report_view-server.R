# testServer harness for the Report view (slice 4A).
#
# Three scenarios:
#   1. NULL project -> module boots without error (demo fixture).
#   2. Project with bundles -> module boots without error.
#   3. Project with no bundles -> module boots without error.
#
# downloadHandler outputs are not directly inspectable in testServer,
# but the boot tests catch signature changes and namespace errors.

test_that("report view boots with demo fixture when project is NULL", {
  current_project <- shiny::reactiveVal(NULL)
  expect_error(
    shiny::testServer(
      report_view_server,
      args = list(current_project = current_project),
      { }
    ),
    NA
  )
})

test_that("report view boots with live project and bundles attached", {
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
  proj$bundles <- list(diff = diff_bundle)
  current_project <- shiny::reactiveVal(proj)
  expect_error(
    shiny::testServer(
      report_view_server,
      args = list(current_project = current_project),
      { }
    ),
    NA
  )
})

test_that("report view boots when project has no bundles yet", {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("readxl")
  xlsx <- tempfile(fileext = ".xlsx")
  on.exit(unlink(xlsx), add = TRUE)
  write_tiny_omics_xlsx(xlsx)
  parsed <- omicsCore::read_omics(xlsx, omics_type = "proteomics")
  proj <- omicsCore::omics_project(
    name        = "empty-proj",
    experiments = list(proteomics = parsed$input)
  )
  current_project <- shiny::reactiveVal(proj)
  expect_error(
    shiny::testServer(
      report_view_server,
      args = list(current_project = current_project),
      { }
    ),
    NA
  )
})
