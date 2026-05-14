#' Launch the omicsApp Shiny application
#'
#' Starts the multi-omics analysis web interface.
#'
#' @param port Integer port number. If `NULL` (default), Shiny chooses a free port.
#' @param host Host address to bind to. Default `"127.0.0.1"` (local only).
#' @param project Optional path to an existing `.omp` project file to open on startup.
#' @param launch.browser Logical. If `TRUE` (default), open the app in the user's
#'   default browser. Set to `FALSE` for headless server environments.
#' @param ... Additional arguments passed to [shiny::runApp()].
#'
#' @return Invisible `NULL`. Called for its side effect (starts the Shiny server).
#' @export
launch <- function(port = NULL,
                   host = "127.0.0.1",
                   project = NULL,
                   launch.browser = interactive(),
                   ...) {
  app_dir <- system.file("app", package = "omicsApp")
  if (!nzchar(app_dir)) {
    stop("Could not find Shiny app directory inside the installed omicsApp package.")
  }

  # Pass project path via options (read by inst/app/app.R)
  if (!is.null(project)) {
    options(omicsApp.project = normalizePath(project, mustWork = TRUE))
  }

  shiny::runApp(
    appDir = app_dir,
    port = port,
    host = host,
    launch.browser = launch.browser,
    ...
  )
}
