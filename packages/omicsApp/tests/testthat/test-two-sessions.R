# Two people on one process.
#
# Outside ShinyProxy -- on a laptop, or a plain Shiny Server -- one R
# process serves every browser tab, and any state a module keeps
# outside its session (a package-level environment, a global) is
# shared between them. testServer() runs one session at a time and
# cannot see that. Two browsers on one app can: what one uploads must
# not appear in the other.
#
# Same gates as the journey test: shinytest2, chromote, a Chrome, and
# the source tree.

test_that("what one session imports, the other does not see", {
  skip_on_cran()
  skip_if_not_installed("shinytest2")
  skip_if_not_installed("chromote")
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("withr")
  skip_if_not(!is.null(tryCatch(chromote::find_chrome(), error = function(e) NULL)),
              "No Chrome/Chromium available for chromote")

  where <- smoke_app_dir()
  skip_if(!nzchar(where$dir), "no app directory to launch")
  store <- file.path(tempfile("sessions-"), "store")
  dir.create(store, recursive = TRUE)
  withr::local_envvar(OMICSAPP_DEV_ROOT = where$dev_root, OMICSAPP_DATA_DIR = store)
  xlsx <- tempfile(fileext = ".xlsx")
  on.exit(unlink(c(xlsx, dirname(store)), recursive = TRUE), add = TRUE)
  write_tiny_omics_xlsx(xlsx, n_features = 40L, n_samples = 8L, seed = 1)

  first <- tryCatch(
    shinytest2::AppDriver$new(where$dir, name = "omicsApp-session-a",
                              load_timeout = 30000, seed = 1),
    error = function(e) skip(sprintf("AppDriver launch failed: %s", conditionMessage(e))))
  on.exit(first$stop(), add = TRUE)
  # A second browser on the same process, not a second process.
  second <- tryCatch(
    shinytest2::AppDriver$new(first$get_url(), name = "omicsApp-session-b",
                              load_timeout = 30000, seed = 2),
    error = function(e) skip(sprintf("second AppDriver failed: %s", conditionMessage(e))))
  on.exit(second$stop(), add = TRUE)

  picker <- function(app) paste(unlist(app$get_value(output = "project_picker")), collapse = " ")

  first$click(selector = "#nav_import")
  first$wait_for_idle(timeout = 5000)
  first$upload_file(`import-file` = xlsx)
  first$wait_for_idle(timeout = 15000)
  first$click("import-confirm")
  first$wait_for_idle(timeout = 15000)
  expect_match(picker(first), "User project", fixed = TRUE)

  # The other browser still has nothing, and its views still show the demo.
  expect_match(picker(second), "no project loaded", fixed = TRUE)
  second$click(selector = "#nav_diff")
  second$wait_for_idle(timeout = 5000)
  expect_match(picker(second), "no project loaded", fixed = TRUE)

  # And running the demo there changes nothing here. The demo result is
  # what that view already shows, so no output need change: wait for
  # the server to go idle, not for an update.
  second$click("diff-rerun", wait_ = FALSE)
  second$wait_for_idle(timeout = 15000)
  expect_match(picker(first), "User project", fixed = TRUE)
  expect_match(picker(second), "no project loaded", fixed = TRUE)
})
