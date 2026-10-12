# The app side of the items carried over from the round-4 review.

test_that("the page declares its own icon, so nothing asks for /favicon.ico", {
  head <- as.character(htmltools::renderTags(app_ui())$head)
  expect_match(head, 'rel="icon"', fixed = TRUE)
  expect_match(head, "data:image/svg+xml,", fixed = TRUE)
})

test_that("static plots are drawn at 96 dpi and their rows grow with wrapped names", {
  expect_identical(PLOT_RES, 96)
  short <- label_rows_px(c("B_vs_A", "C_vs_A"), 14L)
  long <- label_rows_px(c("A treatment with a long descriptive name_vs_A",
                          "Another treatment with a long name_vs_A"), 14L)
  expect_gt(long, short)
})

test_that("a worker's steps reach the notification through its step file", {
  f <- withr::local_tempfile(fileext = ".txt")
  run <- with_step_file(detached_call(function() {
    omicsCore:::report_progress("Estimating dispersions (2 of 4)")
    42
  }), f)
  expect_identical(run(), 42)
  expect_identical(read_step(f), "Estimating dispersions (2 of 4)")
  expect_null(read_step(tempfile()))
  ui <- as.character(async_progress_ui("Running DESeq2...", 12, "t1", step = read_step(f)))
  expect_match(ui, "Estimating dispersions (2 of 4)", fixed = TRUE)
  expect_false(grepl("async-progress-step",
                     as.character(async_progress_ui("Running...", 1, "t2")), fixed = TRUE))
})

test_that("main-process work runs unchanged under the step panel", {
  # testServer has no client for a progress panel; the work still runs.
  val <- NULL
  shiny::testServer(function(input, output, session) {
    val <<- with_step_progress("Running quality control", {
      omicsCore:::report_progress("Counting missing values")
      7
    })
  }, {})
  expect_identical(val, 7)
})

test_that("the global test and enrichment buttons come back disabled mid-run", {
  proj <- shiny::reactiveVal(tutorial_project())
  shiny::testServer(diff_view_server, args = list(current_project = proj), {
    session$setInputs(layer = "proteomics", group_col = "group", control = "Control",
                      case = c("TreatA", "TreatB"), method = "limma", rerun = 1)
    # Re-rendered while a run is in flight (anything the card reads can
    # change meanwhile): it must not come back pressable.
    anova_running(TRUE)
    anova_error("an earlier failure")
    session$flushReact()
    html <- paste(as.character(output$anova_card), collapse = "")
    btn <- regmatches(html, regexpr('<button[^>]*run_anova[^>]*>', html))
    expect_match(btn, "disabled", fixed = TRUE)
    anova_running(FALSE)
    anova_error(NULL)
    session$setInputs(run_anova = 1)
    expect_false(anova_running())
    html <- paste(as.character(output$anova_card), collapse = "")
    btn <- regmatches(html, regexpr('<button[^>]*run_anova[^>]*>', html))
    expect_false(grepl("disabled", btn, fixed = TRUE))
  })
})

test_that("the Top hits card says what to do before the first run", {
  proj <- shiny::reactiveVal(tutorial_project())
  shiny::testServer(diff_view_server, args = list(current_project = proj), {
    session$flushReact()
    expect_match(paste(as.character(output$hits_empty), collapse = ""),
                 "Run the analysis", fixed = TRUE)
    session$setInputs(layer = "proteomics", group_col = "group", control = "Control",
                      case = "TreatA", method = "limma", rerun = 1)
    expect_identical(paste(as.character(output$hits_empty), collapse = ""), "")
  })
})

test_that("cards and tables say log2FC and adjusted p, not effect and adj.P", {
  proj <- shiny::reactiveVal(tutorial_project())
  shiny::testServer(diff_view_server, args = list(current_project = proj), {
    session$setInputs(layer = "proteomics", group_col = "group", control = "Control",
                      case = "TreatA", method = "limma", rerun = 1)
    html <- paste(as.character(output$stats), collapse = "")
    html <- gsub("&gt;", ">", gsub("&lt;", "<", html, fixed = TRUE), fixed = TRUE)
    expect_match(html, "log2FC > 0.26", fixed = TRUE)
    expect_match(html, "adjusted p <", fixed = TRUE)
    expect_false(grepl("effect >", html, fixed = TRUE))
  })
})

test_that("in a browser: no favicon request, and the run button stays on screen", {
  skip_on_cran()
  skip_if_not_installed("shinytest2")
  skip_if_not_installed("chromote")
  skip_if(is.null(tryCatch(chromote::find_chrome(), error = function(e) NULL)),
          "No Chrome/Chromium available for chromote")
  where <- smoke_app_dir()
  if (!nzchar(where$dir)) skip("omicsApp is not installed and no source tree is present.")
  withr::local_envvar(OMICSAPP_DEV_ROOT = where$dev_root,
                      OMICSAPP_DATA_DIR = withr::local_tempdir())
  app <- tryCatch(
    shinytest2::AppDriver$new(where$dir, name = "omicsApp-backlog", width = 1440,
                              height = 900, load_timeout = 60000, timeout = 60000),
    error = function(e) skip(sprintf("AppDriver launch failed: %s", conditionMessage(e))))
  on.exit(app$stop(), add = TRUE)
  app$click("nav_diff")
  app$wait_for_idle(timeout = 60000)
  sticky <- app$get_js(paste0(
    "(() => { const b = document.querySelector('.run-sticky');",
    " return [getComputedStyle(b).position, b.getBoundingClientRect().bottom <= window.innerHeight]; })()"))
  expect_identical(sticky[[1]], "sticky")
  expect_true(sticky[[2]])
  logs <- paste(utils::capture.output(print(app$get_logs())), collapse = "\n")
  expect_false(grepl("favicon.ico", logs, fixed = TRUE))
})

test_that("a burst of project changes is saved once, and a pending one on session end", {
  written <- character(0)
  writer <- function(p) written <<- c(written, p$name)
  mk <- function(nm) { p <- tutorial_project(); p$name <- nm; p }
  shiny::testServer(
    function(input, output, session) {
      current_project <- shiny::reactiveVal(NULL)
      wire_autosave(current_project, writer = writer, delay_ms = 1000)
      shiny::observeEvent(input$burst, {
        current_project(mk("diff"))
      })
      shiny::observeEvent(input$more, current_project(mk("enrich")))
      shiny::observeEvent(input$last, current_project(mk("integration")))
    },
    {
      session$setInputs(burst = 1)
      session$elapse(300)
      session$setInputs(more = 1)
      session$elapse(300)
      expect_length(written, 0L)
      session$elapse(1200)
      expect_identical(written, "enrich")
      session$setInputs(last = 1)
      session$close()
      expect_identical(written, c("enrich", "integration"))
    }
  )
})
