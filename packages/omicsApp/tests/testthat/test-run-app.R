# `launch()` configures the process before handing off to Shiny, and the
# whole point of Phase 0 is that those side effects actually happen.
# `shiny::runApp()` is stubbed so we can inspect the state that is in
# force at the moment the server would have started, without binding a
# port.

launch_state <- function(...) {
  seen <- NULL
  testthat::local_mocked_bindings(
    runApp = function(appDir, port, host, launch.browser, ...) {
      seen <<- list(
        app_dir    = appDir,
        host       = host,
        port       = port,
        max_bytes  = getOption("shiny.maxRequestSize"),
        plan_class = class(future::plan())
      )
      invisible(NULL)
    },
    .package = "shiny"
  )
  omicsApp::launch(launch.browser = FALSE, ...)
  seen
}

test_that("launch passes host and port through to Shiny", {
  seen <- launch_state(host = "0.0.0.0", port = 3838, workers = 0)
  expect_equal(seen$host, "0.0.0.0")
  expect_equal(seen$port, 3838)
  expect_true(nzchar(seen$app_dir))
})

test_that("launch raises the upload ceiling above Shiny's 5 MB default", {
  seen <- launch_state(max_upload_mb = 500, workers = 0)
  expect_equal(seen$max_bytes, 500 * 1024^2)
})

test_that("launch restores the upload option when the app stops", {
  before <- getOption("shiny.maxRequestSize")
  launch_state(max_upload_mb = 123, workers = 0)
  expect_identical(getOption("shiny.maxRequestSize"), before)
})

test_that("workers > 0 installs a parallel future plan", {
  seen <- launch_state(workers = 2)
  # multicore where forking is available (Linux containers), otherwise
  # multisession. Either way it must no longer be sequential, or
  # run_async() silently blocks the event loop.
  expect_false("sequential" %in% seen$plan_class)
  expected <- if (future::supportsMulticore()) "multicore" else "multisession"
  expect_true(expected %in% seen$plan_class)
})

test_that("workers = 0 leaves the caller's plan untouched", {
  old <- future::plan()
  on.exit(future::plan(old), add = TRUE)
  future::plan(future::sequential)
  seen <- launch_state(workers = 0)
  expect_true("sequential" %in% seen$plan_class)
})

test_that("launch restores the previous future plan when the app stops", {
  old <- future::plan()
  on.exit(future::plan(old), add = TRUE)
  future::plan(future::sequential)
  launch_state(workers = 2)
  expect_true("sequential" %in% class(future::plan()))
})

test_that("launch validates its numeric arguments", {
  expect_error(omicsApp::launch(launch.browser = FALSE, workers = -1),
               "non-negative")
  expect_error(omicsApp::launch(launch.browser = FALSE, max_upload_mb = 0),
               "positive")
  expect_error(omicsApp::launch(launch.browser = FALSE, max_upload_mb = NA),
               "positive")
})

# --- launch(project = ) ---------------------------------------------------

# The argument was accepted and did nothing: run_app() wrote the path into
# an option, and nothing in the package read it. `app_server()` started from
# NULL and then restored the autosave, so a caller who named a file got
# somebody else's last session -- and no indication that their argument had
# been ignored.
test_that("a project named at launch is the one that opens", {
  skip_if_not_installed("qs2")
  skip_if_not_installed("withr")
  store <- withr::local_tempdir()
  withr::local_envvar(OMICSAPP_DATA_DIR = store)

  # An autosave exists, so this also shows which of the two wins.
  store_autosave(omicsCore::omics_project("the autosave"))
  named <- withr::local_tempfile(fileext = ".omp")
  omicsCore::save_project(omicsCore::omics_project("the named one"), named)
  withr::local_options(omicsApp.project = named)

  shiny::testServer(app_server, {
    session$flushReact()
    expect_equal(current_project()$name, "the named one")
    # Opening it counts as an activation, like Open and Restore do.
    expect_gt(activation(), 0L)
  })
})

test_that("without the option the autosave is still restored", {
  skip_if_not_installed("qs2")
  skip_if_not_installed("withr")
  store <- withr::local_tempdir()
  withr::local_envvar(OMICSAPP_DATA_DIR = store)
  store_autosave(omicsCore::omics_project("the autosave"))
  withr::local_options(omicsApp.project = NULL)

  shiny::testServer(app_server, {
    session$flushReact()
    expect_equal(current_project()$name, "the autosave")
  })
})

test_that("a named project that cannot be opened is reported, not skipped", {
  skip_if_not_installed("qs2")
  skip_if_not_installed("withr")
  store <- withr::local_tempdir()
  withr::local_envvar(OMICSAPP_DATA_DIR = store)
  store_autosave(omicsCore::omics_project("the autosave"))
  withr::local_options(omicsApp.project = file.path(store, "absent.omp"))

  shiny::testServer(app_server, {
    session$flushReact()
    # Falling through to the autosave without a word would look like the
    # argument had worked.
    expect_false(identical(current_project()$name, "absent"))
  })
})

test_that("launch() does not leave the option set for the next call", {
  skip_if_not_installed("withr")
  f <- withr::local_tempfile(fileext = ".omp")
  omicsCore::save_project(omicsCore::omics_project("x"), f)
  before <- getOption("omicsApp.project")
  # runApp() is what launch() ends on; stub it so the body runs and exits.
  testthat::local_mocked_bindings(runApp = function(...) invisible(NULL),
                                  .package = "shiny")
  launch(project = f, launch.browser = FALSE, workers = 0L)
  expect_identical(getOption("omicsApp.project"), before)
})
