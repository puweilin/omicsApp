# One user's afternoon, in a real browser.
#
# testServer() drives the server function with inputs the test sets
# itself; it cannot see whether the button a user clicks is wired to the
# input the server listens for, whether a modal's Replace button exists,
# or whether a download link produces a file. This is the acceptance
# table again, through Chrome: upload, confirm, run, enrich, download
# the script, then replace the data and watch the results clear.
#
# Same gates as the smoke test: shinytest2, chromote, a Chrome, and the
# source tree (see smoke_app_dir() in test-app-smoke.R).

journey_text <- function(value) {
  if (is.null(value)) return("")
  paste(unlist(value), collapse = " ")
}

test_that("upload, run, enrich, download, replace: the wiring holds end to end", {
  skip_on_cran()
  skip_if_not_installed("shinytest2")
  skip_if_not_installed("chromote")
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("withr")
  skip_if_not(!is.null(tryCatch(chromote::find_chrome(), error = function(e) NULL)),
              "No Chrome/Chromium available for chromote")

  where <- smoke_app_dir()
  skip_if(!nzchar(where$dir), "no app directory to launch")
  store <- file.path(tempfile("journey-"), "store")
  dir.create(store, recursive = TRUE)
  withr::local_envvar(OMICSAPP_DEV_ROOT = where$dev_root,
                      OMICSAPP_DATA_DIR = store)
  xlsx_a <- tempfile(fileext = ".xlsx")
  xlsx_b <- tempfile(fileext = ".xlsx")
  on.exit(unlink(c(xlsx_a, xlsx_b, dirname(store)), recursive = TRUE), add = TRUE)
  write_tiny_omics_xlsx(xlsx_a, n_features = 40L, n_samples = 8L, seed = 1)
  write_tiny_omics_xlsx(xlsx_b, n_features = 40L, n_samples = 8L, seed = 2)

  app <- tryCatch(
    shinytest2::AppDriver$new(where$dir, name = "omicsApp-journey",
                              load_timeout = 30000, seed = 1),
    error = function(e) skip(sprintf("AppDriver launch failed: %s", conditionMessage(e)))
  )
  on.exit(app$stop(), add = TRUE)

  # ---- upload and confirm --------------------------------------------
  app$click(selector = "#nav_import")
  app$wait_for_idle(timeout = 5000)
  app$upload_file(`import-file` = xlsx_a)
  app$wait_for_idle(timeout = 15000)
  expect_match(journey_text(app$get_value(output = "project_picker")),
               "no project loaded", fixed = TRUE)
  app$click("import-confirm")
  app$wait_for_idle(timeout = 15000)
  expect_match(journey_text(app$get_value(output = "project_picker")),
               "User project", fixed = TRUE)
  # The snapshot and the archived upload are on disk already
  expect_true(file.exists(file.path(store, "_autosave.omp")))
  expect_length(list.files(file.path(store, "raw")), 1L)

  # ---- differential --------------------------------------------------
  app$click(selector = "#nav_diff")
  # Entering the view asks has_pkg() about every engine, and the first
  # answer loads DESeq2 and its sixty-odd namespaces into the app
  # process. On a CI runner that alone is longer than the 5 s the other
  # views need.
  app$wait_for_idle(timeout = 30000)
  app$set_inputs(`diff-method` = "ttest")
  app$click("diff-rerun")
  app$wait_for_idle(timeout = 30000)
  stats <- journey_text(app$get_value(output = "diff-stats"))
  expect_match(stats, "Tested features", fixed = TRUE)
  expect_match(stats, "40", fixed = TRUE)

  # ---- enrichment reads the same thresholds --------------------------
  app$click(selector = "#nav_enrich")
  app$wait_for_idle(timeout = 5000)
  summary <- journey_text(app$get_value(output = "enrich-input_summary"))
  expect_match(summary, "40", fixed = TRUE)

  # ---- the script comes out of the browser as a file -----------------
  app$click(selector = "#nav_report")
  app$wait_for_idle(timeout = 5000)
  script_path <- app$get_download("report-download_script")
  expect_true(file.exists(script_path))
  script <- readLines(script_path, warn = FALSE)
  expect_no_error(parse(text = script))
  expect_true(any(grepl("run_diff(", script, fixed = TRUE)))

  # ---- a different file asks before clearing --------------------------
  app$click(selector = "#nav_import")
  app$wait_for_idle(timeout = 5000)
  app$upload_file(`import-file` = xlsx_b)
  app$wait_for_idle(timeout = 15000)
  app$click("import-confirm")
  app$wait_for_idle(timeout = 5000)
  # The modal is in the page, with its Replace button
  expect_true(any(grepl("confirm_replace", app$get_html("body"), fixed = TRUE)))
  app$click("import-confirm_replace")
  app$wait_for_idle(timeout = 15000)
  app$click(selector = "#nav_diff")
  app$wait_for_idle(timeout = 30000)   # same headroom as the first entry
  stats_after <- journey_text(app$get_value(output = "diff-stats"))
  expect_false(grepl("Tested features", stats_after, fixed = TRUE))
  # And the snapshot on disk no longer carries the old analysis
  snap <- omicsCore::load_project(file.path(store, "_autosave.omp"))
  expect_false("diff" %in% names(snap$bundles))
})
