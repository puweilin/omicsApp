# Plot colour palette.
#
# This lived in omicsApp, next to the SCSS it mirrors, back when the
# view modules drew their own figures. Now that the figures come from
# `plot_*()`, the palette has to live with them -- otherwise the app and
# the exported report tint the same data differently, and nobody notices
# until a reviewer compares a screenshot with a PDF.
#
# The values match `omicsApp/inst/app/www/styles.scss`. Change them
# together.

#' Plot colour palette
#'
#' Named colours used by every `plot_*()` function, mirroring the design
#' tokens in the Shiny front end so a figure on screen and the same
#' figure in an exported report are tinted identically.
#'
#' @format A named list of hex colour strings.
#' \describe{
#'   \item{up, down, ns}{Direction of change, and features that reach no
#'     threshold.}
#'   \item{fg_dark, border}{Axis text and rule colours.}
#'   \item{scale_low, scale_high}{Endpoints of continuous colour scales,
#'     low to high significance.}
#'   \item{conc_up_up, conc_down_down, conc_up_down, conc_down_up}{The
#'     four concordance quadrants of a two-omics comparison.}
#'   \item{shared, unique_}{Pathways found in both layers versus one.}
#' }
#' @export
#' @family plot
#' @examples
#' omics_colors$up
omics_colors <- list(
  up        = "#C0392B",
  down      = "#2C3E99",
  ns        = "#9AA3AE",
  fg_dark   = "#1A2541",
  border    = "#E5E7EB",

  # Concordance quadrants. Deliberately not a rainbow: agreement
  # (up_up / down_down) reads as the same warm/cool pair used for
  # direction elsewhere, and disagreement gets the two colours that
  # belong to no other meaning in the app.
  conc_up_up     = "#C0392B",
  conc_down_down = "#1F4E96",
  conc_up_down   = "#E0A030",
  conc_down_up   = "#7A4FA0",

  # Continuous scales run from "not interesting" to "up", so a reader
  # who has learnt the volcano colours reads a dot plot the same way.
  scale_low  = "#9AA3AE",
  scale_high = "#C0392B",

  shared  = "#1F4E96",
  unique_ = "#9AA3AE"
)

#' Colours for the concordance quadrants
#'
#' @return A named character vector covering the four quadrants plus the
#'   fallback used when a feature reaches no threshold in either layer.
#' @keywords internal
#' @noRd
quadrant_palette <- function() {
  c(up_up     = omics_colors$conc_up_up,
    down_down = omics_colors$conc_down_down,
    up_down   = omics_colors$conc_up_down,
    down_up   = omics_colors$conc_down_up,
    ns        = omics_colors$ns,
    `n/a`     = omics_colors$ns)
}
