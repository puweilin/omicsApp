# Shiny entry point for omicsApp.
#
# `omicsApp::launch()` runs `shiny::runApp()` against this directory.
# Everything that defines the UI / server lives in the package's
# `R/` folder so it benefits from roxygen, R CMD check, and reuse
# from tests. Keep this file thin: it returns `omicsApp::shiny_app()`.

# Development hook. When OMICSAPP_DEV_ROOT names the monorepo's
# `packages/` directory, both packages are loaded from source instead of
# from the library. The shinytest2 harness sets it, because otherwise a
# browser test exercises whatever was last installed -- which was four
# months old the day this was added -- and its verdict says nothing
# about the code being changed. omicsCore stays the installed one when
# the test suite settled on that (OMICSAPP_TEST_CORE=installed, see
# tests/testthat/setup.R), so the browser tests and the rest of the
# suite test the same omicsCore.
dev_root <- Sys.getenv("OMICSAPP_DEV_ROOT", "")
if (nzchar(dev_root)) {
  if (!identical(Sys.getenv("OMICSAPP_TEST_CORE"), "installed")) {
    pkgload::load_all(file.path(dev_root, "omicsCore"), quiet = TRUE)
  }
  pkgload::load_all(file.path(dev_root, "omicsApp"), quiet = TRUE)
}

# Through the exported constructor rather than `omicsApp:::app_ui`:
# reaching into the namespace ties this file to internal names that
# R CMD check cannot see being used.
omicsApp::shiny_app()
