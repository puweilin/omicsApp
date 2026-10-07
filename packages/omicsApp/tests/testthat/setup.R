# Test fixtures and global setup for omicsApp.
#
# Kept intentionally small — most setup belongs inside the individual
# test files. We use this file only to silence package-startup chatter
# so test output stays scannable.

# Which omicsCore the tests run against.
#
#   source     omicsApp itself is running from source (devtools::test(),
#              testthat::test_local(), load_all(), a bare test_dir() in
#              the source tree) and the omicsCore source tree is beside
#              it: that tree, loaded with pkgload. `library(omicsCore)`
#              alone loads whatever happens to be installed, and that
#              copy is as old as the last time someone ran
#              install_local(): a day behind and nine tests failed here
#              for reasons that had nothing to do with omicsApp.
#   installed  everything else -- R CMD check, test_check() or
#              test_package() on an installed omicsApp: the installed
#              packages only, and no source is ever loaded. Testing an
#              installed omicsApp against whatever omicsCore happens to
#              sit next to the check directory would test a combination
#              nobody installs. CI installs omicsCore from the same
#              checkout first (.github/workflows/R-CMD-check.yaml).
#
# OMICSAPP_TEST_CORE=installed or =source forces one or the other. The
# answer is written back to OMICSAPP_TEST_CORE so that the browser tests'
# app process (inst/app/app.R) runs the same omicsCore.
core_src <- file.path("..", "..", "..", "omicsCore")
app_src <- file.path("..", "..")
has_pkgload <- requireNamespace("pkgload", quietly = TRUE)
app_from_source <- has_pkgload &&
  if ("omicsApp" %in% loadedNamespaces()) {
    pkgload::is_dev_package("omicsApp")
  } else {
    # A bare test_dir(): omicsApp is loaded below, from the source tree
    # when the tests sit in one.
    file.exists(file.path(app_src, "DESCRIPTION"))
  }
test_core <- tolower(Sys.getenv("OMICSAPP_TEST_CORE", ""))
if (!test_core %in% c("", "installed", "source")) {
  stop("OMICSAPP_TEST_CORE must be 'installed' or 'source', not '", test_core, "'.",
       call. = FALSE)
}
if (!nzchar(test_core)) {
  test_core <- if (app_from_source && file.exists(file.path(core_src, "DESCRIPTION"))) {
    "source"
  } else {
    "installed"
  }
}
if (identical(test_core, "source") &&
    !(has_pkgload && file.exists(file.path(core_src, "DESCRIPTION")))) {
  stop("OMICSAPP_TEST_CORE=source, but there is no omicsCore source tree at ",
       normalizePath(core_src, mustWork = FALSE), " (or pkgload is not installed).",
       call. = FALSE)
}
withr::local_envvar(OMICSAPP_TEST_CORE = test_core, .local_envir = testthat::teardown_env())

suppressPackageStartupMessages({
  library(shiny)
  library(htmltools)
  # Safe to swap in after omicsApp is loaded because omicsApp imports
  # nothing by name -- every call is `omicsCore::`, which resolves
  # against the registered namespace at call time.
  if (identical(test_core, "source") && !pkgload::is_dev_package("omicsCore")) {
    pkgload::load_all(core_src, quiet = TRUE)
  }
  library(omicsCore)
  # omicsApp itself, only when nothing has loaded it yet (a bare
  # test_dir() call). Under devtools::test() it is already loaded, and
  # loading it *again* here would leave every test bound to the old,
  # unregistered namespace: testthat parents the test environments on
  # asNamespace("omicsApp") before it sources this file, so a second
  # load_all() puts local_mocked_bindings() and the code under test in
  # two different copies of the package. The mocks then apply to
  # nothing, silently.
  if (!"omicsApp" %in% loadedNamespaces() && app_from_source) {
    pkgload::load_all(app_src, quiet = TRUE)
  }
  library(omicsApp)
})
rm(core_src, app_src, has_pkgload, app_from_source, test_core)

# Helper: render a tag (or tagList) to HTML for snapshot-style asserts.
render_html <- function(tag) {
  htmltools::renderTags(tag)$html
}

# Server-side tests set an input and read the result in the same flush;
# the QC slider's debounce would defer it.
options(omicsApp.qc_debounce_ms = 0, omicsApp.autosave_debounce_ms = 0)
