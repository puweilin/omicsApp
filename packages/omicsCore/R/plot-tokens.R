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
#'   \item{layer_a, layer_b}{The two layers of an integration, where a
#'     figure shows each layer's value side by side.}
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
  unique_ = "#9AA3AE",

  # The two layers side by side (the top-hits plot). Neither is a
  # direction colour: a blue dot there means "this layer", not "down".
  # Checked as a pair for colour-blind separation; a shape tells them
  # apart as well. Plot-only, so there is no SCSS twin.
  layer_a = "#2A78D6",
  layer_b = "#EB6834"
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

#' Wrap long labels for a plot axis or legend
#'
#' Breaks a label over several lines at spaces, and inside a long word at
#' `_`, `-`, `.` or `/` (a group called `Treatment_high_dose_week_12`
#' has no spaces to break at), cutting a word with none of those where
#' it has to. Past `max_lines` the label ends in an ellipsis. Widths are
#' counted in display columns, so a Chinese label wraps at about half the
#' characters of an English one.
#'
#' @param x Character vector of labels.
#' @param width Widest line, in display columns.
#' @param max_lines Most lines kept per label.
#' @return Character vector the length of `x`, lines joined by `"\n"`.
#' @export
#' @family plot
#' @examples
#' wrap_label("Treatment with a long descriptive name vs Control", width = 20)
wrap_label <- function(x, width = 24L, max_lines = 3L) {
  assert_count(width, "width")
  assert_count(max_lines, "max_lines")
  x <- as.character(x)
  vapply(x, wrap_one_label, character(1), width = width, max_lines = max_lines,
         USE.NAMES = FALSE)
}

# A long label in a fixed-width panel squeezes the panel rather than
# itself: one 60-character group name left the "Hits per comparison"
# bars a 40 px sliver, and pushed the overlap plot off the right edge.
wrap_one_label <- function(s, width, max_lines) {
  if (is.na(s) || nchar(s, type = "width") <= width) return(s)
  cols <- function(t) nchar(t, type = "width")
  # Pieces with the separator that precedes them: " " between words,
  # "" between the parts of a word split at _ - . / or by length.
  pieces <- character(0)
  seps <- character(0)
  for (w in strsplit(s, " ", fixed = TRUE)[[1L]]) {
    if (!nzchar(w)) next
    parts <- if (cols(w) <= width) w else {
      p <- regmatches(w, gregexpr("[^_./-]*[_./-]?", w))[[1L]]
      p[nzchar(p)]
    }
    parts <- unlist(lapply(parts, function(p) {
      if (cols(p) <= width) return(p)
      chars <- strsplit(p, "")[[1L]]
      grp <- cumsum(c(0, utils::head(cumsum(nchar(chars, type = "width")), -1L)) %/% width)
      vapply(split(chars, grp), paste, character(1), collapse = "")
    }), use.names = FALSE)
    pieces <- c(pieces, parts)
    seps <- c(seps, if (length(seps)) " " else "", rep("", length(parts) - 1L))
  }
  lines <- character(0)
  cur <- ""
  for (i in seq_along(pieces)) {
    cand <- paste0(cur, seps[i], pieces[i])
    if (!nzchar(cur) || cols(cand) <= width) {
      cur <- if (nzchar(cur)) cand else pieces[i]
    } else {
      lines <- c(lines, cur)
      cur <- pieces[i]
    }
  }
  lines <- c(lines, cur)
  if (length(lines) > max_lines) {
    lines <- lines[seq_len(max_lines)]
    lines[max_lines] <- paste0(lines[max_lines], "…")
  }
  paste(lines, collapse = "\n")
}

# Scales for a sample-group legend: labels wrapped so a long group name
# does not take the panel's width (at phone width a 60-character name
# left a PCA panel narrower than its legend), and, for up to six groups,
# the same groups drawn as shapes as well as hues -- six default hues
# include a red, a pink and a gold that are hard to tell apart, and are
# the same grey to a reader who cannot see colour.
GROUP_SHAPES <- c(16L, 17L, 15L, 18L, 4L, 8L)

group_legend_scales <- function(values, redundant_shape = TRUE) {
  if (is.null(values) || is.numeric(values)) return(list())
  n <- length(unique(stats::na.omit(as.character(values))))
  lab <- function(x) wrap_label(x, width = 18L, max_lines = 3L)
  # One entry per line wherever the legend is placed: across the bottom
  # of a narrow panel, a row of long group names ran off the edge.
  guide <- ggplot2::guide_legend(ncol = 1)
  out <- list(ggplot2::scale_colour_discrete(labels = lab, guide = guide))
  if (redundant_shape && n <= length(GROUP_SHAPES)) {
    out <- c(out, list(ggplot2::scale_shape_manual(values = GROUP_SHAPES, labels = lab,
                                                   guide = guide)))
  }
  out
}

use_group_shape <- function(values) {
  !is.null(values) && !is.numeric(values) &&
    length(unique(stats::na.omit(as.character(values)))) <= length(GROUP_SHAPES)
}
