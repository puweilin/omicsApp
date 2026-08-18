#' Launch the omicsApp Shiny application
#'
#' Starts the multi-omics analysis web interface.
#'
#' # Running as a multi-user server
#'
#' Two settings matter when omicsApp is served rather than run on a
#' laptop, and both default to laptop-friendly values:
#'
#' * `host = "0.0.0.0"` binds all interfaces so a reverse proxy (or a
#'   container port mapping) can reach the app. Keep the default
#'   `"127.0.0.1"` for local use.
#' * `workers` controls how many background R processes long analyses
#'   are offloaded to. Without this, `future::future()` in
#'   `run_async()` falls back to the `sequential` plan and every
#'   computation blocks the Shiny event loop — one user running
#'   DESeq2 freezes every other session sharing the process.
#'
#' On Linux (including inside a container) `workers > 0` selects
#' `future::multicore`, which forks and therefore shares the parent's
#' expression matrices copy-on-write. Where forking is unavailable or
#' unsafe — Windows, RStudio — it falls back to `future::multisession`,
#' which serialises captured data to each worker instead.
#'
#' @param port Integer port number. If `NULL` (default), Shiny chooses a free port.
#' @param host Host address to bind to. Default `"127.0.0.1"` (local only).
#'   Use `"0.0.0.0"` when running inside a container behind a proxy.
#' @param project Optional path to an existing `.omp` project file to open on startup.
#' @param launch.browser Logical. If `TRUE` (default), open the app in the user's
#'   default browser. Set to `FALSE` for headless server environments.
#' @param workers Integer number of background workers for asynchronous
#'   analyses. `0` keeps the current `future` plan untouched (useful in
#'   tests, where `sequential` is what we want).
#' @param max_upload_mb Maximum accepted upload size in megabytes. Shiny's
#'   built-in default is 5 MB, which is below the size of a typical omics
#'   expression workbook.
#' @param ... Additional arguments passed to [shiny::runApp()].
#'
#' @return Invisible `NULL`. Called for its side effect (starts the Shiny server).
#' @export
launch <- function(port = NULL,
                   host = "127.0.0.1",
                   project = NULL,
                   launch.browser = interactive(),
                   workers = 2L,
                   max_upload_mb = 500,
                   ...) {
  app_dir <- system.file("app", package = "omicsApp")
  if (!nzchar(app_dir)) {
    stop("Could not find Shiny app directory inside the installed omicsApp package.")
  }

  # Pass project path via options (read by inst/app/app.R)
  if (!is.null(project)) {
    options(omicsApp.project = normalizePath(project, mustWork = TRUE))
  }

  if (!is.numeric(max_upload_mb) || length(max_upload_mb) != 1L ||
      is.na(max_upload_mb) || max_upload_mb <= 0) {
    stop("`max_upload_mb` must be a positive number.")
  }
  old_opts <- options(shiny.maxRequestSize = max_upload_mb * 1024^2)
  on.exit(options(old_opts), add = TRUE)

  if (!is.numeric(workers) || length(workers) != 1L || is.na(workers) ||
      workers < 0) {
    stop("`workers` must be a non-negative number.")
  }
  workers <- as.integer(workers)
  if (workers > 0L) {
    # `future::plan()` returns the previous plan, so we can hand it back
    # when the app stops rather than leaving the caller's session
    # reconfigured.
    strategy <- if (future::supportsMulticore()) {
      future::multicore
    } else {
      future::multisession
    }
    previous_plan <- future::plan(strategy, workers = workers)
    on.exit(future::plan(previous_plan), add = TRUE)
  }

  shiny::runApp(
    appDir = app_dir,
    port = port,
    host = host,
    launch.browser = launch.browser,
    ...
  )
}
