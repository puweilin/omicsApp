#' Plot colour palette
#'
#' Mirrors the SCSS design tokens in `inst/app/www/styles.scss` so
#' ggplot2 / plotly calls don't drift away from the rest of the UI
#' if a designer revises the palette. Use these constants instead
#' of raw hex literals in any new module.
#'
#' @keywords internal
#' @noRd
omics_colors <- list(
  up        = "#C0392B",  # $omics-up
  down      = "#2C3E99",  # $omics-down
  ns        = "#9AA3AE",  # $omics-ns
  fg_dark   = "#1A2541",  # axis text on plots
  border    = "#E5E7EB",  # $border-color
  # Integration quadrant palette
  conc_up_up     = "#C0392B",
  conc_down_down = "#1F4E96",
  conc_up_down   = "#E0A030",
  conc_down_up   = "#7A4FA0",
  # Enrichment dotplot scale endpoints
  scale_low  = "#9AA3AE",
  scale_high = "#C0392B",
  # Integration shared-pathway dotplot
  shared  = "#1F4E96",
  unique_ = "#9AA3AE"
)
