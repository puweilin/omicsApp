# Test fixtures and global setup for omicsApp.
#
# Kept intentionally small — most setup belongs inside the individual
# test files. We use this file only to silence package-startup chatter
# so test output stays scannable.

suppressPackageStartupMessages({
  library(shiny)
  library(htmltools)
  # omicsCore from the source tree beside this package, when there is one.
  # `library(omicsCore)` alone loads whatever happens to be installed, and
  # that copy is as old as the last time someone ran install_local(): a
  # day behind and nine tests failed here for reasons that had nothing to
  # do with omicsApp. Safe to swap in after omicsApp is loaded because
  # omicsApp imports nothing by name -- every call is `omicsCore::`, which
  # resolves against the registered namespace at call time.
  core_src <- file.path("..", "..", "..", "omicsCore")
  if (requireNamespace("pkgload", quietly = TRUE) &&
      file.exists(file.path(core_src, "DESCRIPTION")) &&
      !pkgload::is_dev_package("omicsCore")) {
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
  if (requireNamespace("pkgload", quietly = TRUE) &&
      file.exists("../../DESCRIPTION") &&
      !pkgload::is_dev_package("omicsApp")) {
    pkgload::load_all("../..", quiet = TRUE)
  }
})

# Helper: render a tag (or tagList) to HTML for snapshot-style asserts.
render_html <- function(tag) {
  htmltools::renderTags(tag)$html
}
