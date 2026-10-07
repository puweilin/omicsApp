# Shared by the shinytest2 files (test-app-smoke.R, test-app-journey.R).
#
# Where the app under test comes from: the same packages the rest of the
# suite tests (see setup.R). When omicsApp runs from source the app
# directory is `inst/app` and the app process loads omicsApp from source
# via OMICSAPP_DEV_ROOT (see inst/app/app.R) -- and omicsCore too, unless
# setup.R settled on the installed one (OMICSAPP_TEST_CORE=installed).
# The installed copy is for a run against the installed omicsApp (R CMD
# check), and *only* there: it used to be the only option, and the
# installed omicsApp was four months older than the code being tested. A
# browser test of stale code is worse than none, because it says "all
# views render" about a version nobody is editing.
smoke_app_dir <- function() {
  src_app <- file.path("..", "..", "inst", "app")
  if (file.exists(file.path(src_app, "app.R")) &&
      requireNamespace("pkgload", quietly = TRUE) &&
      pkgload::is_dev_package("omicsApp")) {
    packages_dir <- normalizePath(file.path("..", "..", ".."), mustWork = FALSE)
    return(list(dir = normalizePath(src_app), dev_root = packages_dir))
  }
  list(dir = system.file("app", package = "omicsApp"), dev_root = "")
}

