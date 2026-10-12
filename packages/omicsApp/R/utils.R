# Small internal helpers shared across the package.

# Null-coalescing: `a` unless it is NULL. Base R has had `%||%` since 4.4,
# but the package supports R >= 4.2, so it keeps its own. One definition,
# here; the server and each module used to carry a copy.
`%||%` <- function(a, b) if (is.null(a)) b else a

# An omics type as the page names it.
omics_display <- function(t) {
  switch(t %||% "",
         proteomics = "Proteomics",
         rnaseq     = "RNA-seq",
         "\u2014")
}
