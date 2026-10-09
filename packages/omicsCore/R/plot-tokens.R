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

# ---- sample groups -----------------------------------------------------
#
# ggplot's default hues coloured groups until this: the first group --
# usually the control -- came out salmon red, next to the red that means
# "up", and the third a teal-blue next to "down", so a reader could take
# a group for a direction. Six of them were also too close for a reader
# with deuteranopia to tell apart.
#
# These eight were chosen against the Machado (2009) simulation of
# protanopia and deuteranopia (distances in OKLab x 100, the dataviz
# palette rules):
#   - the first six are at least 8.8 apart from one another, every pair,
#     with and without colour vision deficiency (17 without), so any of
#     them can sit beside any other in a PCA;
#   - seven and eight are at least 8.8 from their neighbours in the order
#     (a heatmap's group bar puts neighbours side by side; eight's
#     neighbour past eight groups is one), but closer to some of the
#     first six; past six groups the PCA's shapes carry the difference
#     too;
#   - every one is at least 21 from the up red and the down blue (11
#     under the simulations), so none of them reads as a direction.
#     There is no red, and no blue as dark as the down blue: under
#     deuteranopia the up red turns a dark olive and a mid green turns
#     the same olive, so there is no mid green either.
# They sit at OKLCH lightness 0.59-0.76, lighter than most plot colours,
# because the dark end is where up and down live. Lighter colours have
# less contrast with the white panel (2.1:1 at worst), so a group is
# never told by colour alone: the PCA's legend names it and gives it a
# shape as well, the boxplot's axis names it, and so does the heatmap's
# legend.
#
# Plot-only, so there is no SCSS twin. Checked in test-plot-tokens.R.
GROUP_PALETTE <- c(
  teal       = "#0A9282",
  orange     = "#F59A0B",
  violet     = "#9B59D0",
  sky        = "#3AB0F5",
  pink       = "#EE7FC0",
  green      = "#5DC77E",
  periwinkle = "#8387E3",
  olive      = "#9F9615"
)

#' Colours for sample groups
#'
#' The colours every `plot_*()` function gives sample groups: the PCA's
#' points, the selected feature's boxplot, the heatmap's group bar. The
#' eight colours stay apart for readers with red-green colour blindness,
#' and none of them is a red or a dark blue, so a group never looks like
#' a direction of change (`omics_colors$up` and `omics_colors$down`).
#'
#' The figures deal them out in the study design's group order -- the
#' reference group first, then the others -- over every group of the
#' column, not only the ones a figure shows, so a group has the same
#' colour in every figure. Past eight groups the colours repeat, and the
#' PCA tells the ninth group from the first by its shape.
#'
#' @param n Number of groups.
#' @return Character vector of `n` hex colours.
#' @export
#' @family plot
#' @examples
#' group_palette(4)
group_palette <- function(n) {
  assert_count(n, "n")
  unname(GROUP_PALETTE[(seq_len(n) - 1L) %% length(GROUP_PALETTE) + 1L])
}

# The groups of a column in the order their colours are dealt: the
# reference first, then a factor's levels or else sorted. Taken over the
# whole column, so a figure that shows two of four groups, or shows the
# comparison's reference first, colours each group as the PCA does.
group_order <- function(values, reference = NULL) {
  lv <- if (is.factor(values)) levels(droplevels(values))
        else sort(unique(as.character(stats::na.omit(values))))
  ref <- as.character(reference)
  if (length(ref) == 1L && !is.na(ref) && ref %in% lv) lv <- c(ref, setdiff(lv, ref))
  lv
}

# One colour per group, named by group.
group_colours <- function(values, reference = NULL) {
  lv <- group_order(values, reference)
  stats::setNames(group_palette(length(lv)), lv)
}

# The reference group recorded on a layer, when `column` is the layer's
# group column; NULL otherwise (a figure coloured by batch has none).
design_reference <- function(input, column) {
  d <- if (is_omics_input(input) && !is.null(column)) {
    tryCatch(study_design(input), error = function(e) NULL)
  }
  if (!is.null(d) && identical(d$group_col, column)) d$reference
}

# Shapes for a sample-group legend, the second channel beside colour: up
# to eight groups each have their own. Past eight, where the colours
# repeat, the first eight groups are circles, the next eight triangles,
# and so on, so no two of the first 64 groups share colour and shape.
GROUP_SHAPES <- c(16L, 17L, 15L, 18L, 4L, 8L, 1L, 2L)

group_shapes <- function(n) {
  if (n <= length(GROUP_SHAPES)) return(GROUP_SHAPES[seq_len(n)])
  GROUP_SHAPES[((seq_len(n) - 1L) %/% length(GROUP_PALETTE)) %% length(GROUP_SHAPES) + 1L]
}

# Scales for a sample-group legend: the group colours, labels wrapped so
# a long group name does not take the panel's width (at phone width a
# 60-character name left a PCA panel narrower than its legend), and with
# `redundant_shape` the same groups drawn as shapes as well. `values` is
# the whole column (see group_order()); the legend lists the groups
# drawn, in the colours' order.
group_legend_scales <- function(values, redundant_shape = TRUE, reference = NULL) {
  if (is.null(values) || is.numeric(values)) return(list())
  pal <- group_colours(values, reference)
  lab <- function(x) wrap_label(x, width = 18L, max_lines = 3L)
  # One entry per line wherever the legend is placed: across the bottom
  # of a narrow panel, a row of long group names ran off the edge.
  guide <- ggplot2::guide_legend(ncol = 1)
  out <- list(ggplot2::scale_colour_manual(values = pal, breaks = names(pal),
                                           labels = lab, guide = guide))
  if (redundant_shape) {
    shapes <- stats::setNames(group_shapes(length(pal)), names(pal))
    out <- c(out, list(ggplot2::scale_shape_manual(values = shapes, breaks = names(pal),
                                                   labels = lab, guide = guide)))
  }
  out
}

# Every discrete grouping is drawn with shapes as well: up to eight
# groups as a second channel beside colour, past eight to tell apart the
# groups whose colours repeat.
use_group_shape <- function(values) {
  !is.null(values) && !is.numeric(values) &&
    length(stats::na.omit(values)) > 0L
}
