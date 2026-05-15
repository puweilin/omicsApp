# Extended edge-case tests for view modules.

# ---- import view tests -------------------------------------------------

test_that("import view boots without error", {
  expect_error(
    shiny::testServer(
      import_view_server,
      args = list(id = "test_import"),
      {
        session$setInputs(file = NULL)
        expect_true(TRUE)
      }
    ), NA)
})

test_that("import view can toggle format", {
  expect_error(
    shiny::testServer(
      import_view_server,
      args = list(id = "test_import"),
      {
        session$setInputs(file = NULL, confirm = 0)
        expect_true(TRUE)
      }
    ), NA)
})

# ---- qc view tests ----------------------------------------------------

test_that("qc view boots without error with NULL project", {
  current_project <- shiny::reactiveVal(NULL)
  expect_error(
    shiny::testServer(
      qc_view_server,
      args = list(current_project = current_project),
      {
        session$setInputs(missing_threshold = 0.5, outlier_method = "pca")
        expect_true(TRUE)
      }
    ), NA)
})

test_that("qc view handles project with single experiment", {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("readxl")

  xlsx <- tempfile(fileext = ".xlsx")
  on.exit(unlink(xlsx), add = TRUE)
  write_tiny_omics_xlsx(xlsx)
  parsed <- omicsCore::read_omics(xlsx, omics_type = "proteomics")
  proj <- omicsCore::omics_project(
    name = "test",
    experiments = list(proteomics = parsed$input)
  )
  current_project <- shiny::reactiveVal(proj)

  expect_error(
    shiny::testServer(
      qc_view_server,
      args = list(current_project = current_project),
      {
        session$setInputs(
          missing_threshold = 0.5,
          outlier_method = "pca",
          rerun = 1
        )
        expect_true(TRUE)
      }
    ), NA)
})

test_that("qc view can switch outlier methods without crashing", {
  current_project <- shiny::reactiveVal(NULL)
  expect_error(
    shiny::testServer(
      qc_view_server,
      args = list(current_project = current_project),
      {
        for (method in c("pca", "connectivity", "iqr")) {
          session$setInputs(outlier_method = method)
        }
        expect_true(TRUE)
      }
    ), NA)
})

# ---- project view tests ------------------------------------------------

test_that("project view boots with NULL project", {
  current_project <- shiny::reactiveVal(NULL)
  expect_error(
    shiny::testServer(
      project_view_server,
      args = list(current_project = current_project),
      { expect_true(TRUE) }
    ), NA)
})

test_that("project view boots with one experiment", {
  inp <- example_input("proteomics")
  proj <- omicsCore::omics_project(
    name = "test_proj",
    experiments = list(proteomics = inp)
  )
  current_project <- shiny::reactiveVal(proj)
  expect_error(
    shiny::testServer(
      project_view_server,
      args = list(current_project = current_project),
      { expect_true(TRUE) }
    ), NA)
})

test_that("project view boots with two experiments", {
  inp <- example_input("proteomics")
  proj <- omicsCore::omics_project(
    name = "dual",
    experiments = list(prot = inp, rna = example_input("rnaseq"))
  )
  current_project <- shiny::reactiveVal(proj)
  expect_error(
    shiny::testServer(
      project_view_server,
      args = list(current_project = current_project),
      { expect_true(TRUE) }
    ), NA)
})

# ---- report view tests -------------------------------------------------

test_that("report view boots without crashing with NULL project", {
  current_project <- shiny::reactiveVal(NULL)
  expect_error(
    shiny::testServer(
      report_view_server,
      args = list(current_project = current_project),
      { expect_true(TRUE) }
    ), NA)
})

test_that("report view boots with project but no bundles", {
  inp <- example_input("proteomics")
  proj <- omicsCore::omics_project(
    name = "test",
    experiments = list(prot = inp)
  )
  current_project <- shiny::reactiveVal(proj)
  expect_error(
    shiny::testServer(
      report_view_server,
      args = list(current_project = current_project),
      { expect_true(TRUE) }
    ), NA)
})

test_that("report view boots with project and bundles", {
  inp <- example_input("proteomics")
  proj <- omicsCore::omics_project(
    name = "test",
    experiments = list(prot = inp)
  )
  bundle <- run_diff(inp, method = "ttest", analysis_type = "group",
                     group_col = "group",
                     control_group = "G1", case_group = "G2")
  proj$bundles <- list(diff = bundle)
  current_project <- shiny::reactiveVal(proj)
  expect_error(
    shiny::testServer(
      report_view_server,
      args = list(current_project = current_project),
      { expect_true(TRUE) }
    ), NA)
})
