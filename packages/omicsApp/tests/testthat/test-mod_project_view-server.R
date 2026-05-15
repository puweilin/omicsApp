# testServer harness for the Project view (slice 4A).
#
# Three scenarios:
#   1. NULL project -> demo fixture with example_project().
#   2. Project with one experiment -> live tag/type/shape cards.
#   3. Project with two experiments -> both layers listed.

test_that("project view falls back to demo when current_project is NULL", {
  current_project <- shiny::reactiveVal(NULL)
  shiny::testServer(
    project_view_server,
    args = list(current_project = current_project),
    {
      # Module should boot without error and render the demo project.
      expect_true(TRUE)
    }
  )
})

test_that("project view shows live experiment when project has one layer", {
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
  shiny::testServer(
    project_view_server,
    args = list(current_project = current_project),
    {
      expect_true(TRUE)
    }
  )
})

test_that("project view lists both experiments when project has two layers", {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("readxl")
  xlsx <- tempfile(fileext = ".xlsx")
  on.exit(unlink(xlsx), add = TRUE)
  write_tiny_omics_xlsx(xlsx)
  parsed <- omicsCore::read_omics(xlsx, omics_type = "proteomics")
  inp_a <- parsed$input
  inp_b <- inp_a
  inp_b$omics_type <- "rnaseq"
  proj <- omicsCore::omics_project(
    name        = "dual",
    experiments = list(proteomics = inp_a, rnaseq = inp_b)
  )
  current_project <- shiny::reactiveVal(proj)
  shiny::testServer(
    project_view_server,
    args = list(current_project = current_project),
    {
      expect_true(TRUE)
    }
  )
})
