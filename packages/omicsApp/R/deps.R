# One seam for "is this optional package here?"
#
# The Report view checks for rmarkdown, the Enrichment view for
# clusterProfiler, the Differential view for whichever engine a layer
# needs. Each asked `requireNamespace()` directly, and on a machine where
# everything is installed -- every developer laptop, every CI runner --
# the branch that tells the user what to install could not be reached,
# so it was never run. Routing the question through here lets a test
# answer it.
#
# @param pkg Package name.
# @return `TRUE` when the namespace can be loaded.
# @keywords internal
# @noRd
has_pkg <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}
