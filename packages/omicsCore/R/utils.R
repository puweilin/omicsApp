# Small internal helpers shared across the package.

# Null-coalescing: `a` unless it is NULL. Base R has had `%||%` since 4.4,
# but the package supports R >= 4.2, and rlang's would be one more import
# to keep in step. One definition, here; two files used to carry a copy.
`%||%` <- function(a, b) if (is.null(a)) b else a

# A plot with nothing to show but a sentence saying why: a result with
# no rows draws this rather than an empty panel or an error.
empty_plot <- function(label) {
  ggplot2::ggplot() +
    ggplot2::theme_void() +
    ggplot2::annotate("text", x = 0.5, y = 0.5, label = label,
                      color = "#4D4D4D", size = 4) +
    ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1)
}
