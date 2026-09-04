# Shared by the shinytest2 files (test-app-smoke.R, test-app-journey.R).
#
# Where the app under test comes from. In the source tree the app
# directory is `inst/app` and both packages are loaded from source via
# OMICSAPP_DEV_ROOT (see inst/app/app.R). The installed copy is the
# fallback for R CMD check, and *only* there: it used to be the only
# option, and the installed omicsApp was four months older than the
# code being tested. A browser test of stale code is worse than none,
# because it says "all views render" about a version nobody is editing.
smoke_app_dir <- function() {
  src_app <- file.path("..", "..", "inst", "app")
  packages_dir <- normalizePath(file.path("..", "..", ".."), mustWork = FALSE)
  if (file.exists(file.path(src_app, "app.R")) &&
      file.exists(file.path(packages_dir, "omicsCore", "DESCRIPTION")) &&
      requireNamespace("pkgload", quietly = TRUE)) {
    return(list(dir = normalizePath(src_app), dev_root = packages_dir))
  }
  list(dir = system.file("app", package = "omicsApp"), dev_root = "")
}

