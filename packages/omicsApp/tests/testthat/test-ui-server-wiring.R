# What the server listens for is what the UI draws.
#
# A control whose id the server never reads does nothing when clicked,
# and an `input$` the UI never declares waits forever; neither fails a
# test that drives the server function directly, because testServer()
# sets whatever inputs the test names. Read from the source instead:
# every `input$x` in a module's server must be an `ns("x")` in the same
# file, every `output$y` must have its `ns("y")` placeholder, and every
# control the UI draws must be read back.
#
# `input` is also the name modules give an omics_input in helper
# functions, so `input$meta_df` inside `function(input)` is not a
# Shiny input, and `active()$input$omics_type` is not one either. The
# walk below follows the parse tree and keeps track of both: a
# function whose formals are `input, output, session` is a server, and
# any other function with an `input` formal hides the Shiny one.

module_files <- function() {
  r_dir <- file.path("..", "..", "R")
  skip_if_not(dir.exists(r_dir), "package source is not beside the tests")
  list.files(r_dir, pattern = "^mod_.*\\.R$", full.names = TRUE)
}

# One row per reference: kind = input / output / ns, the id, and the
# name of the call the reference sits in.
collect_ids <- function(file) {
  rows <- list()
  note <- function(kind, id, caller) {
    rows[[length(rows) + 1L]] <<- data.frame(kind = kind, id = id, caller = caller,
                                             stringsAsFactors = FALSE)
  }
  call_name <- function(fn) {
    if (is.symbol(fn)) return(as.character(fn))
    if (is.call(fn) && identical(fn[[1L]], as.name("::"))) return(as.character(fn[[3L]]))
    ""
  }
  walk <- function(expr, shadowed, caller) {
    if (!is.call(expr)) return(invisible())
    fn <- expr[[1L]]
    if (identical(fn, as.name("function"))) {
      formals_here <- names(expr[[2L]])
      is_server <- all(c("input", "output", "session") %in% formals_here)
      hides <- "input" %in% formals_here && !is_server
      walk(expr[[3L]], if (is_server) FALSE else shadowed || hides, caller)
      return(invisible())
    }
    if (identical(fn, as.name("$")) && length(expr) == 3L && is.symbol(expr[[2L]])) {
      who <- as.character(expr[[2L]])
      if (who == "input" && !shadowed) note("input", as.character(expr[[3L]]), caller)
      if (who == "output") note("output", as.character(expr[[3L]]), caller)
    }
    if (identical(fn, as.name("[[")) && length(expr) == 3L && is.symbol(expr[[2L]]) &&
        is.character(expr[[3L]])) {
      who <- as.character(expr[[2L]])
      if (who == "input" && !shadowed) note("input", expr[[3L]], caller)
      if (who == "output") note("output", expr[[3L]], caller)
    }
    is_ns <- identical(fn, as.name("ns")) ||
      (is.call(fn) && identical(fn[[1L]], as.name("$")) &&
         identical(fn[[2L]], as.name("session")) && identical(fn[[3L]], as.name("ns")))
    if (is_ns && length(expr) == 2L && is.character(expr[[2L]])) {
      note("ns", expr[[2L]], caller)
    }
    here <- call_name(fn)
    for (i in seq_along(expr)[-1L]) {
      if (!is.null(expr[[i]])) walk(expr[[i]], shadowed, if (nzchar(here)) here else caller)
    }
    invisible()
  }
  # Parsed without a locale warning about the comments' typography.
  for (e in suppressWarnings(parse(file, keep.source = FALSE))) walk(e, FALSE, "")
  if (length(rows) == 0L) {
    return(data.frame(kind = character(0), id = character(0), caller = character(0)))
  }
  unique(do.call(rbind, rows))
}

# Ids built at run time from a pattern, on both sides, so there is no
# literal to match. Named here so the check stays exact everywhere else.
dynamic_ids <- list(
  # The report's download buttons come from a helper that takes the id.
  mod_report_view.R = c("download_html", "download_pdf")
)

control_re <- "(Input|Button|Buttons|Link|Switch|Slider|Checkbox|Picker)$"
download_re <- "^download(Button|Link)$"
output_re <- "Output$"

test_that("every input a module reads is a control it draws", {
  for (file in module_files()) {
    ids <- collect_ids(file)
    declared <- c(ids$id[ids$kind == "ns"], dynamic_ids[[basename(file)]])
    read <- unique(ids$id[ids$kind == "input"])
    missing <- setdiff(read, declared)
    expect_identical(missing, character(0),
                     label = sprintf("%s reads inputs it never declares", basename(file)))
  }
})

test_that("every output a module assigns has a placeholder", {
  for (file in module_files()) {
    ids <- collect_ids(file)
    declared <- c(ids$id[ids$kind == "ns"], dynamic_ids[[basename(file)]])
    assigned <- unique(ids$id[ids$kind == "output"])
    missing <- setdiff(assigned, declared)
    expect_identical(missing, character(0),
                     label = sprintf("%s assigns outputs it never places", basename(file)))
  }
})

test_that("every control a module draws is read back, and every placeholder filled", {
  for (file in module_files()) {
    ids <- collect_ids(file)
    ns_ids <- ids[ids$kind == "ns", , drop = FALSE]
    is_download <- grepl(download_re, ns_ids$caller)
    controls <- ns_ids$id[grepl(control_re, ns_ids$caller) & !is_download]
    placeholders <- ns_ids$id[grepl(output_re, ns_ids$caller) | is_download]
    read <- ids$id[ids$kind == "input"]
    assigned <- ids$id[ids$kind == "output"]
    expect_identical(setdiff(controls, read), character(0),
                     label = sprintf("%s draws controls nothing reads", basename(file)))
    expect_identical(setdiff(placeholders, assigned), character(0),
                     label = sprintf("%s places outputs nothing fills", basename(file)))
  }
})

test_that("the top level wires its nav links and mounts each module under one id", {
  r_dir <- file.path("..", "..", "R")
  skip_if_not(dir.exists(r_dir), "package source is not beside the tests")
  ui <- readLines(file.path(r_dir, "app_ui.R"), warn = FALSE)
  server <- readLines(file.path(r_dir, "app_server.R"), warn = FALSE)
  nav_ids <- regmatches(ui, regexpr('nav_item\\("[a-z_]+"', ui))
  nav_ids <- sub('nav_item\\("', "", nav_ids)
  nav_ids <- sub('"$', "", nav_ids)
  read <- unique(unlist(regmatches(server, gregexpr("input\\$nav_[a-z_]+", server))))
  read <- sub("input\\$", "", read)
  expect_setequal(read, nav_ids)
  expect_true(any(grepl('uiOutput\\("project_picker"\\)', ui)))
  expect_true(any(grepl("output\\$project_picker", server)))
  # `xxx_view_ui("id")` in the UI and `xxx_view_server("id", ...)` in the
  # server, module by module.
  ui_mounts <- regmatches(ui, regexpr('[a-z]+_view_ui\\("[a-z]+"', ui))
  server_mounts <- regmatches(server, regexpr('[a-z]+_view_server\\("[a-z]+"', server))
  expect_length(ui_mounts, 7L)
  expect_setequal(sub("_ui", "", ui_mounts), sub("_server", "", server_mounts))
  # The tab that shows a view has the view's id.
  tabs <- regmatches(ui, regexpr('tabPanelBody\\("[a-z]+"', ui))
  expect_setequal(sub('tabPanelBody\\("', "", tabs), sub('.*\\("', "", ui_mounts))
})

test_that("the walker tells a Shiny input from an omics_input", {
  tmp <- withr::local_tempfile(fileext = ".R")
  writeLines(c(
    "srv <- function(id) {",
    "  shiny::moduleServer(id, function(input, output, session) {",
    "    ns <- session$ns",
    "    output$x <- shiny::renderUI(shiny::textInput(ns('a'), 'A'))",
    "    shiny::observe(input$a)",
    "    shiny::observe(input[['b']])",
    "    pick <- function(input) input$meta_df",
    "    shiny::observe(active()$input$omics_type)",
    "    output$y <- shiny::renderUI(session$ns('c'))",
    "  })",
    "}"), tmp)
  ids <- collect_ids(tmp)
  expect_setequal(ids$id[ids$kind == "input"], c("a", "b"))
  expect_setequal(ids$id[ids$kind == "output"], c("x", "y"))
  expect_setequal(ids$id[ids$kind == "ns"], c("a", "c"))
  expect_identical(ids$caller[ids$kind == "ns" & ids$id == "a"], "textInput")
})
