test_that("applying QC updates downstream data without losing the original layer", {
  inp <- example_proteomics_input()
  inp$expr_mat[1, ] <- NA_real_
  current_project <- shiny::reactiveVal(omicsCore::omics_project("desktop", list(protein = inp)))
  shiny::testServer(qc_view_server, args = list(current_project = current_project), {
    session$setInputs(missing_threshold = 0.5, outlier_method = "none", impute_method = "min")
    expect_null(last_error())
    session$setInputs(apply_to_downstream = 1)
    p <- current_project()
    expect_equal(nrow(p$experiments$protein$expr_mat), nrow(inp$expr_mat) - 1L)
    expect_identical(active()$input$expr_mat, inp$expr_mat)
    expect_false(anyNA(p$experiments$protein$expr_mat))
  })
})

test_that("an explicit startup project wins over the previous autosave", {
  store <- withr::local_tempdir()
  withr::local_envvar(OMICSAPP_DATA_DIR = store)
  chosen <- omicsCore::omics_project("chosen", list(protein = example_proteomics_input()))
  prior <- chosen
  prior$name <- "prior autosave"
  store_autosave(prior)
  path <- withr::local_tempfile(fileext = ".omp")
  omicsCore::save_project(chosen, path)
  withr::local_options(omicsApp.project = path)
  p <- shiny::reactiveVal(NULL)
  shiny::testServer(project_view_server, args = list(current_project = p), {
    session$flushReact()
    expect_equal(p()$name, "chosen")
  })
})

test_that("restored differential settings populate controls then allow user changes", {
  inp <- example_proteomics_input()
  p <- omicsCore::omics_project("saved settings", list(protein = inp))
  p$bundles <- list(diff = omicsCore::run_diff(inp, method = "limma", group_col = "group",
                                             control_group = "G1", case_group = "G2"))
  p$bundles$diff$input_info$layer <- "protein"
  current_project <- shiny::reactiveVal(NULL)
  shiny::testServer(diff_view_server, args = list(current_project = current_project), {
    session$setInputs(method = "auto")
    current_project(p)
    session$flushReact()
    expect_equal(pick("method", input$method), "limma")
    expect_match(output$ui_method$html, 'value="limma" selected', fixed = TRUE)
    session$setInputs(method = "limma")
    session$setInputs(method = "auto")
    expect_equal(pick("method", input$method), "auto")
  })
})
