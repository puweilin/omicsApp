#' Enrichment summary plot
#'
#' Visualises the standardized enrichment table produced by
#' [run_enrichment()]. The `"dot"` view draws a dotplot of the top `top_n`
#' pathways per database: for ORA, x is the share of each pathway's genes
#' found in the gene list and point size the number found; for GSEA, x is
#' the NES and point size the gene set size. Colour is -log10 of the
#' chosen p-value, capped so that one extreme pathway does not wash out
#' the rest: values above the cap take the top colour, and the top
#' legend label is marked as "at least" that value. The
#' `"bar"` view draws a horizontal bar chart of the same selection. For
#' GSEA bundles a `"gsea_dot"` view is available that splits pathways by
#' direction (up / down).
#'
#' An ORA run with up- and down-regulated genes tested separately (the
#' default of [run_enrichment()]) is drawn with one panel per gene list,
#' headed "Up-regulated genes" and "Down-regulated genes", and `top_n`
#' is shared between the two lists of a database (alternating best of
#' each), so neither direction crowds the other out.
#'
#' @param bundle An [`analysis_bundle`][is_analysis_bundle()] produced by
#'   [run_enrichment()].
#' @param top_n Number of pathways to display per database.
#' @param view One of `"dot"` (default), `"bar"`, or `"gsea_dot"`.
#' @param p_preference `"adjusted"` (default), `"raw"`, or `"qvalue"`.
#' @param p_cutoff Optional significance cutoff. If `NULL`, all rows are
#'   shown (subject to `top_n`).
#' @param database Optional vector of databases to restrict to.
#'
#' @return A `ggplot` object.
#' @export
#' @family enrich
plot_enrichment <- function(
  bundle,
  top_n = 20L,
  view = c("dot", "bar", "gsea_dot"),
  p_preference = c("adjusted", "raw", "qvalue"),
  p_cutoff = NULL,
  database = NULL
) {
  assert_count(top_n, "top_n")
  assert_number(p_cutoff, "p_cutoff", lower = 0, upper = 1, allow_null = TRUE)
  assert_character(database, "database", allow_null = TRUE)
  view <- match.arg(view)
  p_preference <- match.arg(p_preference)
  df <- enrich_result_from_bundle(bundle)

  if (!is.null(database)) {
    keep_db <- vapply(database, normalize_enrich_database, character(1L))
    df <- df[df$database %in% keep_db, , drop = FALSE]
  }
  if (!is.null(p_cutoff)) {
    df <- filter_enrich_results(df, p_cutoff = p_cutoff, p_preference = p_preference)
  }
  if (nrow(df) == 0L) {
    return(empty_plot("No pathways to plot."))
  }

  p_col <- resolve_enrich_p_col(df, p_preference)
  df <- pick_top_per_database(df, top_n, p_col)

  switch(view,
    dot      = plot_enrich_dot(df, p_col),
    bar      = plot_enrich_bar(df, p_col),
    gsea_dot = plot_enrich_gsea_dot(df, p_col)
  )
}

#' GSEA running-score plot
#'
#' Plots the `clusterProfiler::gseaplot2()` running-score curve for a single
#' pathway in a GSEA enrichment bundle. Requires the `enrichplot`
#' Bioconductor package.
#'
#' @param bundle An [`analysis_bundle`][is_analysis_bundle()] produced by
#'   [run_enrichment()] with `type = "gsea"`.
#' @param pathway_id Pathway ID (matches `pathway_id` in the standardized
#'   table) or pathway name.
#' @param database Database key. Required when the bundle holds multiple
#'   databases; ignored when there's only one.
#'
#' @return A `ggplot` object.
#' @export
#' @family enrich
plot_gsea <- function(bundle, pathway_id, database = NULL) {
  if (!is_analysis_bundle(bundle) ||
      !identical(bundle$analysis_name, "run_enrichment")) {
    stop("`bundle` must be an analysis_bundle from run_enrichment().")
  }
  if (!identical(bundle$params$type, "gsea")) {
    stop("plot_gsea() requires a bundle with type = 'gsea'.")
  }
  if (!is_installed("enrichplot")) {
    stop(
      "Package 'enrichplot' is required for plot_gsea(). ",
      "Install with: BiocManager::install('enrichplot').",
      call. = FALSE
    )
  }

  objects <- bundle$results$enrich_object
  if (length(objects) == 0L) {
    stop("Bundle does not contain any GSEA objects.")
  }

  obj <- if (length(objects) == 1L) {
    objects[[1L]]
  } else {
    if (is.null(database)) {
      stop("`database` must be supplied when the bundle covers multiple databases: ",
           paste(names(objects), collapse = ", "))
    }
    objects[[normalize_enrich_database(database)]]
  }
  if (is.null(obj)) {
    stop("No GSEA object available for the requested database.")
  }

  enrichplot::gseaplot2(obj, geneSetID = pathway_id)
}

# ---- internal helpers --------------------------------------------------

enrich_result_from_bundle <- function(bundle) {
  if (!is_analysis_bundle(bundle) ||
      !identical(bundle$analysis_name, "run_enrichment")) {
    stop("`bundle` must be an analysis_bundle from run_enrichment().")
  }
  df <- bundle$results$enrich_result_df
  if (is.null(df)) stop("Bundle is missing `results$enrich_result_df`.")
  check_enrich_result_schema(df)
  df
}

pick_top_per_database <- function(df, top_n, p_col) {
  df <- df[!is.na(df[[p_col]]), , drop = FALSE]
  if (nrow(df) == 0L) return(df)
  df <- df[order(df$database, df[[p_col]]), , drop = FALSE]
  split_df <- split(df, df$database, drop = TRUE)
  picked <- lapply(split_df, function(x) {
    # Up and down lists take turns -- best of each, then second best of
    # each -- so a strong up list does not leave the down list unseen,
    # and a short list hands its unused places to the other.
    if (all(c("up", "down") %in% x$direction) && all(x$direction %in% c("up", "down")) &&
        ora_list_facet(x)) {
      turn <- stats::ave(seq_len(nrow(x)), x$direction, FUN = seq_along)
      x <- x[order(turn, x[[p_col]]), , drop = FALSE]
    }
    utils::head(x, top_n)
  })
  out <- do.call(rbind, picked)
  rownames(out) <- NULL
  out
}

# ORA rows that say which gene list they came from: the panels are then
# split by list. ORA has no effect size; a result that has one (GSEA's
# NES) carries its direction as that sign, which the x axis already
# shows.
ora_list_facet <- function(df) {
  !any(is.finite(df$effect)) && any(df$direction %in% c("up", "down"))
}

ORA_LIST_LABELS <- c(up = "Up-regulated genes", down = "Down-regulated genes")

with_list_label <- function(df) {
  lab <- ifelse(df$direction %in% names(ORA_LIST_LABELS),
                ORA_LIST_LABELS[df$direction], "Up and down pooled")
  df$.list <- factor(lab, levels = c(unname(ORA_LIST_LABELS), "Up and down pooled"))
  df
}

enrich_facets <- function(df, scales) {
  if (ora_list_facet(df)) {
    ggplot2::facet_wrap(ggplot2::vars(.data$database, .data$.list),
                        scales = scales, ncol = 1)
  } else {
    ggplot2::facet_wrap(~ .data$database, scales = scales, ncol = 1)
  }
}


plot_enrich_dot <- function(df, p_col) {
  if (ora_list_facet(df)) df <- with_list_label(df)
  df$.label <- wrap_pathway_name(df$pathway_name)
  df$.label <- factor(df$.label, levels = unique(df$.label[order(-df[[p_col]])]))
  df$.signif <- neg_log10_p(df[[p_col]])

  usable <- function(col) col %in% names(df) && any(is.finite(df[[col]]))

  # Each aesthetic carries a different number. ORA once put the overlap
  # on x *and* on point size: the same count twice, with a size legend
  # that said nothing new and took a third of a phone's width.
  #
  # Effect (NES for GSEA) says which way a pathway moved, so it has x
  # when the result carries it. ORA has no effect; its x is the share of
  # the pathway's genes that turned up in the list -- a 10-gene overlap
  # means more in a 15-gene pathway than in a 500-gene one, and the
  # count alone cannot say which. The count stays as point size.
  has_effect <- usable("effect")
  ratio_ok <- !has_effect && usable("overlap_size") && usable("gene_set_size")
  if (ratio_ok) {
    df$.ratio <- df$overlap_size / df$gene_set_size
    # A pathway without its size would drop out of the plot without a
    # word; better the whole panel falls back than loses a row.
    ratio_ok <- all(is.finite(df$.ratio)) && all(df$gene_set_size > 0)
  }
  x_aes <- if (has_effect) "effect" else if (ratio_ok) ".ratio" else ".signif"
  colour_signif <- !identical(x_aes, ".signif")

  # The size aesthetic needs a fallback too, for GSEA. GSEA has no
  # overlap: standardize_enrich_results only fills overlap_size from
  # `Count` or `GeneRatio`, which fgsea emits neither of, so the column
  # is entirely NA. Mapping size to it made every point NA-sized, and
  # geom_point drops those -- the panel drew its axes and facet strips
  # and not one dot.
  #
  # That reads as "enrichment found nothing", which is the same thing an
  # empty result looks like, so the failure hides as a result.
  size_aes <- if (usable("overlap_size")) "overlap_size" else "gene_set_size"
  has_size <- usable(size_aes)

  aes_args <- list(x = quote(.data[[x_aes]]), y = quote(.data$.label))
  if (colour_signif) aes_args$color <- quote(.data$.signif)
  if (has_size) aes_args$size <- quote(.data[[size_aes]])

  # When x already is the significance (no effect, no usable ratio),
  # colouring by it as well would be the same number twice again: the
  # points take one plain colour instead.
  point <- if (colour_signif) {
    ggplot2::geom_point(na.rm = TRUE)
  } else {
    ggplot2::geom_point(na.rm = TRUE, color = omics_colors$scale_high)
  }

  x_label <- if (has_effect) {
    nes <- "effect_type" %in% names(df) && all(df$effect_type %in% "nes")
    if (nes) "normalized enrichment score (NES)" else "effect"
  } else if (ratio_ok) {
    "% of pathway genes in the list"
  } else {
    p_axis_label(p_col)
  }

  p <- ggplot2::ggplot(df, do.call(ggplot2::aes, aes_args)) +
    point +
    # x is shared between panels so the up and down lists, or two
    # databases, can be compared along it; only the names differ.
    enrich_facets(df, scales = "free_y") +
    ggplot2::labs(title = "Enrichment", x = x_label, y = NULL) +
    theme_omics_labelled() +
    enrich_narrow_theme()

  if (colour_signif) p <- p + signif_colour_scale(df$.signif, p_col)

  if (ratio_ok) {
    # From zero, so a share reads as a share; the extra room on the
    # right keeps the largest dot from being cut by the panel edge. Few
    # ticks: beside long pathway names the panel can be 120 px wide, and
    # "0% 25% 50% 75% 100%" ran together there.
    p <- p + ggplot2::expand_limits(x = 0) +
      ggplot2::scale_x_continuous(labels = scales::label_percent(),
                                  breaks = scales::breaks_extended(n = 3),
                                  expand = ggplot2::expansion(mult = c(0.03, 0.1)))
  } else {
    p <- p + ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.08, 0.1)))
  }

  if (has_size) {
    # Named for the column actually mapped, so a GSEA panel does not
    # label a gene-set size as genes found in a list. Three keys at
    # most: under a phone-width panel the legend runs across, and a
    # fourth large dot ran off the edge.
    p <- p + ggplot2::scale_size_continuous(
      range = c(2, 6), breaks = size_breaks,
      guide = ggplot2::guide_legend(order = 2L),
      name = if (identical(size_aes, "overlap_size")) "genes in list" else "set size")
  }

  if (has_effect) {
    p <- p + ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                                 color = omics_colors$ns)
  }
  p
}

# Pathway names take more than half of a phone-width plot, which leaves
# the panel -- and anything centred on it -- a narrow strip at the right.
# The x title, centred there, ran off the edge; so did the colour bar's
# top label and the size key once the app moved the legends below the
# panel. The title ends where the panel ends, and legends under the plot
# start at its left edge and use its whole width. On a desktop card the
# legends stay at the side and nothing changes but the title's place.
enrich_narrow_theme <- function() {
  th <- ggplot2::theme(axis.title.x = ggplot2::element_text(hjust = 1))
  tree <- names(ggplot2::get_element_tree())
  # Both appeared in ggplot2 3.5.
  if (all(c("legend.location", "legend.justification.bottom") %in% tree)) {
    th <- th + ggplot2::theme(legend.location = "plot",
                              legend.justification.bottom = "left")
  }
  th
}

# Two or three round breaks for a size key. Asked for three, the
# pretty breaks of 10..91 are 0, 50, 100, and only one of those falls
# inside the range: a key of a single dot says nothing about size. So
# ask for more until two land inside, and keep at most three of them.
# When every point has the same size the range has no width, no round
# break lands in it, and the key shows that one value.
size_breaks <- function(limits) {
  b <- numeric(0)
  for (n in 3:6) {
    b <- scales::breaks_extended(n = n)(limits)
    b <- unique(b[b >= limits[1L] & b <= limits[2L]])
    if (length(b) >= 2L) break
  }
  # Every other one (or every third), so the keys stay evenly spaced.
  if (length(b) > 3L) b <- b[seq(1L, length(b), by = ceiling((length(b) - 1L) / 2L))]
  if (length(b) == 0L) unique(limits) else b
}

neg_log10_p <- function(p) -log10(pmax(p, .Machine$double.xmin))

# Colour for significance, -log10(p), low to high in the same colours as
# the volcano: a reader who has learnt one scale reads the other without
# relearning it.
#
# The top of the scale is capped. Left to the data, one pathway at
# -log10 p = 200 set the top colour, and every other pathway -- the
# ones between 2 and 15 a reader actually has to tell apart -- came out
# the same grey. Values above the cap are drawn in the top colour
# (squished, not dropped) and the top label is the cap with a
# greater-or-equal sign, so the legend says the scale stops there.
signif_colour_scale <- function(values, p_col) {
  lim <- signif_limits(values)
  breaks <- ggplot2::waiver()
  labels <- ggplot2::waiver()
  if (lim$capped) {
    lo <- lim$limits[1L]
    cap <- lim$limits[2L]
    # Round breaks below the cap, kept clear of it so the two labels do
    # not run together, then the cap itself.
    b <- scales::breaks_extended(n = 4)(lim$limits)
    b <- b[b >= lo & b <= cap - 0.25 * (cap - lo)]
    breaks <- c(b, cap)
    labels <- c(format_signif_break(b),
                paste0("\u2265 ", format_signif_break(cap)))
  }
  ggplot2::scale_color_gradient(low = omics_colors$scale_low,
                                high = omics_colors$scale_high,
                                name = p_axis_label(p_col),
                                limits = lim$limits, oob = scales::squish,
                                breaks = breaks, labels = labels,
                                guide = ggplot2::guide_colourbar(order = 1L))
}

format_signif_break <- function(x) as.character(round(x, 1))

# The cap is the largest value that is not an outlier: Tukey's upper
# fence (Q3 + 1.5 IQR) on the plotted values, or, when the values are
# too few for quartiles to find it, a top value more than four times the
# next. A quantile on its own does not work here: the panel shows a
# dozen pathways, and the 90% quantile of twelve values interpolates
# most of the way to the extreme one.
#
# Nothing is capped when nothing stands out, so the ">=" only appears
# when it is true of some point. One pathway, or all at the same p,
# gives a range of zero; the scale then starts at 0 so the points still
# take a colour.
signif_limits <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(list(limits = NULL, capped = FALSE))
  lo <- min(x)
  hi <- max(x)
  q <- stats::quantile(x, c(0.25, 0.75), names = FALSE)
  cap <- max(x[x <= q[2L] + 1.5 * (q[2L] - q[1L])])
  if (cap == hi && length(x) >= 3L && any(x < hi)) {
    second <- max(x[x < hi])
    if (hi > 4 * second) cap <- second
  }
  capped <- cap < hi && cap > lo
  if (capped) {
    # A round number for the label: ">= 15", not ">= 15.35".
    nice <- if (cap >= 10) floor(cap) else floor(cap * 10) / 10
    if (nice > lo) cap <- nice
  } else {
    cap <- hi
  }
  if (cap <= lo) {
    lo <- 0
    if (cap <= 0) cap <- 1
  }
  list(limits = c(lo, cap), capped = capped)
}

plot_enrich_bar <- function(df, p_col) {
  if (ora_list_facet(df)) df <- with_list_label(df)
  df$.label <- wrap_pathway_name(df$pathway_name)
  df$.label <- factor(df$.label, levels = unique(df$.label[order(df[[p_col]], decreasing = TRUE)]))
  ggplot2::ggplot(
    df,
    ggplot2::aes(x = -log10(pmax(.data[[p_col]], .Machine$double.xmin)),
                 y = .data$.label, fill = .data$database)
  ) +
    ggplot2::geom_col() +
    enrich_facets(df, scales = "free_y") +
    ggplot2::guides(fill = "none") +
    ggplot2::labs(
      title = "Enrichment",
      x = p_axis_label(p_col),
      y = NULL
    ) +
    theme_omics_labelled()
}

plot_enrich_gsea_dot <- function(df, p_col) {
  if (!"direction" %in% colnames(df) || all(is.na(df$direction))) {
    stop("`gsea_dot` view requires a `direction` column with non-NA values.")
  }
  df$.label <- wrap_pathway_name(df$pathway_name)
  df$.label <- factor(df$.label, levels = unique(df$.label[order(-df[[p_col]])]))
  df$.signif <- neg_log10_p(df[[p_col]])
  ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data$effect, y = .data$.label,
                 size = .data$gene_set_size, color = .data$.signif)
  ) +
    ggplot2::geom_point(na.rm = TRUE) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = omics_colors$ns) +
    ggplot2::facet_grid(rows = ggplot2::vars(.data$database),
                        cols = ggplot2::vars(.data$direction),
                        scales = "free_y") +
    # The same capped significance scale as the dot view. It coloured
    # the raw p before, red for small and grey for large -- the reverse
    # of every other plot in the app -- under the column name.
    signif_colour_scale(df$.signif, p_col) +
    ggplot2::scale_size_continuous(range = c(2, 6), breaks = size_breaks,
                                   guide = ggplot2::guide_legend(order = 2L),
                                   name = "set size") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.08, 0.1))) +
    ggplot2::labs(
      title = "GSEA",
      x = "NES",
      y = NULL
    ) +
    theme_omics_labelled() +
    enrich_narrow_theme()
}

# Wrapped, not cut short. A truncated GO BP term is often
# indistinguishable from three others that share its opening words,
# which is the one thing a pathway label has to avoid; wrapping keeps
# the whole name at the cost of vertical space, and vertical space is
# what this panel has.
#
# `max_lines` is the backstop. Without it a single 200-character
# Reactome term would claim the height of four others, and the plot
# would be about the label rather than the result.
# Kept alongside wrap_pathway_name() for plot-integration.R, which
# applies it to gene symbols. A symbol is short and atomic: wrapping one
# across two lines would be worse than shortening it.
truncate_pathway_name <- function(x, max_chars = 45L) {
  x <- as.character(x)
  too_long <- !is.na(x) & nchar(x) > max_chars
  x[too_long] <- paste0(substr(x[too_long], 1L, max_chars - 1L), "\u2026")
  x
}

wrap_pathway_name <- function(x, width = 34L, max_lines = 3L) {
  wrap_label(x, width = width, max_lines = max_lines)
}

# Used by the enrichment and integration panels -- named for why, not
# for which module, since it now serves both.
#
# These panels are read, not glanced at: the y
# axis carries pathway names or gene symbols, and those *are* the
# result. theme_omicsCore()'s 11pt base leaves them smaller than the
# title, which is backwards, and smaller still once the panel is one of
# two sharing a row.
#
# The base size moves rather than only axis.text, so the legend, the
# axis titles and the facet strips come with it -- raising one and
# leaving the rest is how a plot ends up looking mismatched.
theme_omics_labelled <- function(base_size = 13) {
  theme_omicsCore(base_size = base_size) +
    ggplot2::theme(
      axis.text.y  = ggplot2::element_text(size = ggplot2::rel(1.05),
                                           lineheight = 0.95),
      legend.text  = ggplot2::element_text(size = ggplot2::rel(0.9)),
      legend.title = ggplot2::element_text(size = ggplot2::rel(0.9))
    )
}
