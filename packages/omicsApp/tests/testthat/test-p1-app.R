# The remaining P1 items on the app side.

test_that("a result whose settings were changed afterwards says so until re-run", {
  proj <- shiny::reactiveVal(tutorial_project())
  shiny::testServer(diff_view_server, args = list(current_project = proj), {
    stale_note <- function() grepl("settings have changed",
                                   paste(as.character(output$notices), collapse = ""), fixed = TRUE)
    session$setInputs(layer = "proteomics", group_col = "group", control = "Control",
                      case = "TreatA", method = "limma", rerun = 1)
    expect_false(stale_note())
    session$setInputs(case = c("TreatA", "TreatB"))
    expect_true(stale_note())
    # Back to what was run: no longer out of date.
    session$setInputs(case = "TreatA")
    expect_false(stale_note())
    session$setInputs(covariates = "batch")
    expect_true(stale_note())
    session$setInputs(rerun = 2)
    expect_false(stale_note())
  })
})

test_that("a restored result is not flagged until a control is changed", {
  p <- tutorial_project()
  p$bundles <- list(diff = omicsCore::run_diff(p$experiments$proteomics, method = "limma",
                                               group_col = "group", control_group = "Control",
                                               case_group = "TreatA"))
  proj <- shiny::reactiveVal(NULL)
  shiny::testServer(diff_view_server, args = list(current_project = proj), {
    proj(p)
    session$flushReact()
    expect_false(is.null(diff_bundle()))
    html <- paste(as.character(output$notices), collapse = "")
    expect_false(grepl("settings have changed", html, fixed = TRUE))
    session$setInputs(method = "ttest")
    expect_match(paste(as.character(output$notices), collapse = ""), "settings have changed",
                 fixed = TRUE)
  })
})
