# Test fixtures and global setup for omicsApp.
#
# Kept intentionally small — most setup belongs inside the individual
# test files. We use this file only to silence package-startup chatter
# so test output stays scannable.

suppressPackageStartupMessages({
  library(shiny)
  library(htmltools)
  library(omicsCore)
  if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all("../..", quiet = TRUE)
  }
})

# Helper: render a tag (or tagList) to HTML for snapshot-style asserts.
render_html <- function(tag) {
  htmltools::renderTags(tag)$html
}
