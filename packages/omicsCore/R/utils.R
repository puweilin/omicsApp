# Small internal helpers shared across the package.

# Null-coalescing: `a` unless it is NULL. Base R has had `%||%` since 4.4,
# but the package supports R >= 4.2, and rlang's would be one more import
# to keep in step. One definition, here; two files used to carry a copy.
`%||%` <- function(a, b) if (is.null(a)) b else a
