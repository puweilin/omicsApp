#' Plot colour palette
#'
#' Re-exported from omicsCore, which is where the palette moved when the
#' view modules stopped drawing their own figures. Keeping the name
#' bound here means the app's remaining non-plot uses -- legend
#' swatches, inline chips -- read the same as before.
#'
#' The SCSS design tokens in `inst/app/www/styles.scss` mirror these
#' values; change them together.
#'
#' @keywords internal
#' @noRd
omics_colors <- omicsCore::omics_colors
