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
# Asks whether the package is installed, not whether it loads:
# requireNamespace() loaded it, and drawing the Differential view's
# method list loaded DESeq2 -- six seconds -- before anything was asked
# of it. A package that is installed and broken still fails, later, with
# its own message at the point of use.
#
# @param pkg Package name.
# @return `TRUE` when the package is installed.
# @keywords internal
# @noRd
has_pkg <- function(pkg) {
  nzchar(system.file(package = pkg))
}
