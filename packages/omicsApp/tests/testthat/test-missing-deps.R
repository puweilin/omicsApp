# What each view says when an optional package is not there.
#
# These branches exist for the collaborator on a locked-down machine who
# installed the light package and nothing else. They could not be tested
# before: every developer machine has everything, and the one test that
# tried was skipped with "clusterProfiler is installed" on every run.
# The views now ask has_pkg(), which a test can answer.

absent <- function(...) {
  gone <- c(...)
  testthat::local_mocked_bindings(
    has_pkg = function(pkg) {
      !(pkg %in% gone) && requireNamespace(pkg, quietly = TRUE)
    },
    .package = "omicsApp",
    .env = parent.frame()
  )
}

test_that("the enrichment view points at the missing package instead of computing", {
  absent("clusterProfiler")
  diff_bundle <- shiny::reactiveVal(example_diff_bundle())
  shiny::testServer(
    enrich_view_server,
    args = list(diff_bundle = diff_bundle),
    {
      session$setInputs(type = "ora", database = "hallmark",
                        direction = "both", rerun = 1)
      expect_true(isTRUE(is_demo()))
      expect_null(enrich_bundle())
      expect_match(enrich_error(), "clusterProfiler", fixed = TRUE)
      expect_match(enrich_error(), "install", ignore.case = TRUE)
    }
  )
})

test_that("the report view offers no download without rmarkdown", {
  absent("rmarkdown")
  shiny::testServer(
    report_view_server,
    args = list(current_project = shiny::reactiveVal(NULL)),
    {
      html <- render_html(output$notices)
      expect_match(html, "rmarkdown unavailable", fixed = TRUE)
      expect_match(html, "install.packages", fixed = TRUE)
      # And the buttons are there but disabled, not gone
      head_html <- render_html(output$header)
      expect_match(head_html, "disabled", fixed = TRUE)
    }
  )
})

test_that("the differential view names the engines it cannot offer", {
  absent("DESeq2", "edgeR")
  expect_setequal(diff_missing_engines(), c("DESeq2", "edgeR"))

  absent("DESeq2", "edgeR", "limma")
  expect_setequal(diff_missing_engines(), c("limma", "DESeq2", "edgeR"))
})

test_that("with nothing missing the differential view reports nothing missing", {
  skip_if_not_installed("limma")
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("edgeR")
  expect_length(diff_missing_engines(), 0L)
})

test_that("the template writer refuses without openxlsx", {
  absent("openxlsx")
  expect_error(write_import_template(tempfile(fileext = ".xlsx")),
               "openxlsx")
})

test_that("every optional-package question in the app goes through has_pkg()", {
  # The mock above only reaches code that asks has_pkg(). A view that
  # goes back to requireNamespace() directly silently leaves this file.
  src_dir <- system.file("R", package = "omicsApp") |> dirname() |>
    file.path("R")
  src <- list.files(src_dir, pattern = "\\.R$", full.names = TRUE)
  skip_if(length(src) == 0L, "package sources not available")
  hits <- unlist(lapply(src, function(f) {
    lines <- readLines(f, warn = FALSE)
    code <- lines[!grepl("^\\s*#", lines)]
    if (any(grepl("requireNamespace(", code, fixed = TRUE))) basename(f)
  }))
  expect_identical(hits, "deps.R")
})
