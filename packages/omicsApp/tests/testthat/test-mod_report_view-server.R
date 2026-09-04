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
  parsed <- omicsCore::read_omics(xlsx, omics_type = "proteomics", assay_type = "normalized_intensity")
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
  parsed <- omicsCore::read_omics(xlsx, omics_type = "proteomics", assay_type = "normalized_intensity")
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

# ---- reproducibility script ------------------------------------------

test_that("the code panel is empty until there is a project", {
  shiny::testServer(
    report_view_server,
    args = list(current_project = shiny::reactiveVal(NULL)),
    {
      expect_null(script_lines())
      html <- render_html(output$script_card)
      expect_match(html, "Import data and run an analysis", fixed = TRUE)
    }
  )
})

test_that("the code panel shows the calls that produced the project", {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("readxl")
  xlsx <- tempfile(fileext = ".xlsx")
  on.exit(unlink(xlsx), add = TRUE)
  write_tiny_omics_xlsx(xlsx)
  inp <- omicsCore::read_omics(xlsx, omics_type = "proteomics", assay_type = "normalized_intensity")$input
  inp$source_path <- xlsx
  proj <- omicsCore::omics_project("panel", experiments = list(proteomics = inp))
  proj$bundles <- list()

  shiny::testServer(
    report_view_server,
    args = list(current_project = shiny::reactiveVal(proj)),
    {
      lines <- script_lines()
      expect_true(any(grepl("library(omicsCore)", lines, fixed = TRUE)))
      expect_true(any(grepl("read_omics(", lines, fixed = TRUE)))
      html <- render_html(output$script_card)
      expect_match(html, "Analysis code", fixed = TRUE)
      expect_match(html, "Download .R", fixed = TRUE)
    }
  )
})

test_that("a project that cannot be rendered shows a comment, not a crash", {
  # export_script() rejects a non-project; the panel must degrade to a
  # message rather than taking the Report view down with it.
  shiny::testServer(
    report_view_server,
    args = list(current_project = shiny::reactiveVal(structure(
      list(name = "broken"), class = "omics_project"))),
    {
      expect_no_error(lines <- script_lines())
      expect_true(is.character(lines))
    }
  )
})
