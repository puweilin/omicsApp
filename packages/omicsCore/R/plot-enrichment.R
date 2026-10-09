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
#' `"bar"` view draws a horizontal bar chart of the same selection, bar
#' length -log10 p, coloured by gene list (ORA: up red, down blue) or by
#' the sign of the NES (GSEA); a bar more than three times the longest
#' ordinary one is cut short, marked with a break and labelled with its
#' true value. For
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
#' The running enrichment score of one pathway in a GSEA result: the
#' genes ranked from most up to most down along x, a tick for each gene
#' of the pathway in a band under the curve, and the curve climbing at
#' each tick and falling between them. Its peak (marked) is the
#' enrichment score: a peak near the left means the pathway's genes
#' gather among the most up-regulated, a trough near the right among the
#' most down-regulated.
#'
#' The ranking and the gene sets are read from the clusterProfiler object
#' the run keeps (`results$enrich_object`). A bundle without it is drawn
#' from the differential result it was run on, passed as `diff_bundle`:
#' the ranking is rebuilt as [run_enrichment()] builds it and the gene
#' set read again from the database.
#'
#' @param bundle An [`analysis_bundle`][is_analysis_bundle()] produced by
#'   [run_enrichment()] with `type = "gsea"`.
#' @param pathway_id Pathway ID (matches `pathway_id` in the standardized
#'   table) or pathway name.
#' @param database Database key. Required when the bundle holds multiple
#'   databases; ignored when there's only one.
#' @param diff_bundle Optional: the differential bundle the enrichment
#'   was run on. Used only when `bundle` does not keep its GSEA objects.
#'
#' @return A `ggplot` object.
#' @export
#' @family enrich
plot_gsea <- function(bundle, pathway_id, database = NULL, diff_bundle = NULL) {
  if (!is_analysis_bundle(bundle) ||
      !identical(bundle$analysis_name, "run_enrichment")) {
    stop("`bundle` must be an analysis_bundle from run_enrichment().")
  }
  if (!identical(bundle$params$type, "gsea")) {
    stop("plot_gsea() requires a bundle with type = 'gsea'.")
  }
  assert_string(pathway_id, "pathway_id")
  dbs <- bundle$params$database %||% names(bundle$results$enrich_object)
  if (is.null(database)) {
    if (length(dbs) > 1L) {
      stop("`database` must be supplied when the bundle covers multiple databases: ",
           paste(dbs, collapse = ", "), call. = FALSE)
    }
    database <- dbs[[1L]]
  }
  database <- normalize_enrich_database(database)

  # The pathway's row of the table, for its name, NES and p. Found by ID
  # or by name, as the caller has it.
  df <- bundle$results$enrich_result_df
  row <- NULL
  if (is.data.frame(df) && nrow(df)) {
    in_db <- is.na(df$database) | df$database == database
    hit <- which(in_db & (df$pathway_id %in% pathway_id | df$pathway_name %in% pathway_id))
    if (length(hit)) row <- df[hit[[1L]], , drop = FALSE]
  }
  set_id <- if (!is.null(row)) row$pathway_id else pathway_id

  run <- gsea_ranking_and_set(bundle, database, set_id, diff_bundle)
  curve <- gsea_running_score(run$ranking, run$gene_set, exponent = run$exponent)
  if (!any(curve$hit)) {
    stop("None of the genes of ", pathway_id, " are in the ranked list.", call. = FALSE)
  }
  gsea_curve_plot(curve, row = row, title = pathway_id,
                  metric = bundle$params$rank_metric %||% attr(run$ranking, "metric"))
}

# The ranked statistics and the pathway's genes: from the stored
# clusterProfiler object, or without one rebuilt from the differential
# result the way run_enrichment() built them (symbol case included).
gsea_ranking_and_set <- function(bundle, database, set_id, diff_bundle = NULL) {
  objects <- bundle$results$enrich_object
  obj <- if (length(objects)) {
    objects[[database]] %||% (if (length(objects) == 1L) objects[[1L]])
  }
  if (!is.null(obj) && methods::is(obj, "gseaResult")) {
    sets <- methods::slot(obj, "geneSets")
    if (!set_id %in% names(sets)) {
      stop("No gene set called ", set_id, " in this result.", call. = FALSE)
    }
    return(list(ranking = methods::slot(obj, "geneList"),
                gene_set = sets[[set_id]],
                exponent = methods::slot(obj, "params")$exponent %||% 1))
  }
  if (is.null(diff_bundle)) {
    stop("This result does not keep its ranked gene list; pass the differential ",
         "result it was run on (`diff_bundle`) to draw the curve.", call. = FALSE)
  }
  assert_diff_bundle(diff_bundle)
  organism <- bundle$params$organism %||% "Hs"
  if (identical(bundle$params$symbol_case, "ignored")) {
    fix <- match_symbol_case(diff_bundle, database, organism)
    if (!is.null(fix)) diff_bundle <- fix$bundle
  }
  res <- diff_result_from_bundle(diff_bundle)
  col <- if ("feature_symbol" %in% colnames(res)) "feature_symbol" else "feature_id"
  t2g <- build_term_tables(database = database, organism = organism)$term2gene
  genes <- unique(t2g$gene[t2g$term == set_id])
  if (!length(genes)) {
    stop("No gene set called ", set_id, " in ", database, ".", call. = FALSE)
  }
  list(ranking = gsea_rank_vector(res, col), gene_set = genes, exponent = 1)
}

# The running score as GSEA computes it (Subramanian et al. 2005; the
# same sums as DOSE::gseaScores(), which clusterProfiler uses). Walking
# down the ranking, a gene of the pathway adds its share of the
# pathway's summed |statistic|^exponent and any other gene takes away
# 1 / (genes outside the pathway), so the walk ends at 0 and the
# enrichment score is its largest excursion either way.
gsea_running_score <- function(ranking, gene_set, exponent = 1) {
  ranking <- sort(ranking[is.finite(ranking)], decreasing = TRUE)
  hit <- names(ranking) %in% gene_set
  n <- length(ranking)
  n_hit <- sum(hit)
  w <- ifelse(hit, abs(ranking)^exponent, 0)
  up <- if (sum(w) > 0) cumsum(w / sum(w)) else cumsum(hit / max(n_hit, 1L))
  down <- cumsum((!hit) / max(n - n_hit, 1L))
  data.frame(position = seq_len(n), gene = names(ranking), stat = unname(ranking),
             score = up - down, hit = hit, stringsAsFactors = FALSE)
}

# One panel, not gseaplot2's three. Its bright green curve, rainbow bar
# and "Ranked List Metric" panel used none of the app's colours, and at
# phone width its two axis titles ran into each other. Here the curve
# takes the direction's colour, the pathway's genes are ticks in a band
# under it, and a strip under the ticks says which end of the ranking is
# up and which down -- the one thing a reader needs to read the rest.
gsea_curve_plot <- function(curve, row = NULL, title = NULL, metric = NULL) {
  n <- nrow(curve)
  peak <- which.max(abs(curve$score))
  es <- curve$score[peak]
  nes <- if (!is.null(row)) row$effect else NA_real_
  colour <- if (if (is.finite(nes)) nes >= 0 else es >= 0) omics_colors$up else omics_colors$down

  lo <- min(0, curve$score)
  hi <- max(0, curve$score)
  span <- max(hi - lo, 1e-6)
  tick_top <- lo - 0.06 * span
  tick_bot <- tick_top - 0.16 * span
  strip_top <- tick_bot - 0.05 * span
  strip_bot <- strip_top - 0.08 * span
  n_up <- sum(curve$stat > 0)
  strip <- data.frame(xmin = c(0.5, n_up + 0.5), xmax = c(n_up + 0.5, n + 0.5),
                      fill = c(omics_colors$up, omics_colors$down),
                      label = c("up", "down"), x = c(1, n), hjust = c(-0.2, 1.2),
                      stringsAsFactors = FALSE)
  strip <- strip[strip$xmax > strip$xmin, , drop = FALSE]
  ticks <- curve[curve$hit, , drop = FALSE]

  name <- if (!is.null(row)) row$pathway_name else title %||% "Pathway"
  n_genes <- sum(curve$hit)
  sub <- if (!is.null(row) && is.finite(nes)) {
    sprintf("NES %s \u00B7 adjusted p %s \u00B7 %d genes",
            formatC(nes, digits = 2, format = "f"), format_enrich_p(row$adj_p_value), n_genes)
  } else {
    sprintf("%d genes in the set", n_genes)
  }
  metric_txt <- switch(metric %||% "",
                       "signed test statistic" = "test statistic",
                       "signed sqrt(F)" = "signed \u221AF",
                       "sign(effect) * -log10(p)" = "signed -log10 p",
                       "change")
  score_breaks <- function(l) {
    b <- scales::breaks_extended(n = 4)(c(lo, hi))
    b[b >= lo - 0.05 * span & b <= hi + 0.05 * span]
  }

  ggplot2::ggplot(curve, ggplot2::aes(x = .data$position)) +
    ggplot2::geom_hline(yintercept = 0, colour = omics_colors$ns, linewidth = 0.4) +
    ggplot2::annotate("segment", x = peak, xend = peak, y = 0, yend = es,
                      colour = colour, linetype = "dashed", linewidth = 0.4) +
    ggplot2::geom_line(ggplot2::aes(y = .data$score), colour = colour, linewidth = 0.9) +
    ggplot2::annotate("point", x = peak, y = es, colour = colour, size = 2.2) +
    # Beside the dashed line, half-way up it: at the peak itself the label
    # sat on the curve wherever the score plateaus there.
    ggplot2::annotate("text", x = peak, y = es / 2,
                      label = paste("ES", formatC(es, digits = 2, format = "f")),
                      hjust = if (peak < n / 2) -0.25 else 1.25,
                      vjust = 0.5,
                      size = 3.4, colour = omics_colors$fg_dark) +
    ggplot2::geom_segment(data = ticks,
                          ggplot2::aes(x = .data$position, xend = .data$position),
                          y = tick_bot, yend = tick_top,
                          colour = omics_colors$fg_dark, linewidth = 0.3, alpha = 0.7) +
    ggplot2::geom_rect(data = strip,
                       ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax),
                       ymin = strip_bot, ymax = strip_top,
                       fill = strip$fill, alpha = 0.3, inherit.aes = FALSE) +
    ggplot2::geom_text(data = strip,
                       ggplot2::aes(x = .data$x, label = .data$label, hjust = .data$hjust),
                       y = (strip_top + strip_bot) / 2, size = 3,
                       colour = omics_colors$fg_dark, inherit.aes = FALSE) +
    # Numbers for the score only; the band and the strip carry none.
    ggplot2::scale_y_continuous(breaks = score_breaks,
                                expand = ggplot2::expansion(mult = c(0.02, 0.08))) +
    ggplot2::scale_x_continuous(breaks = NULL, expand = ggplot2::expansion(mult = 0.01)) +
    ggplot2::labs(title = wrap_label(name, width = 40L, max_lines = 2L),
                  subtitle = sub,
                  x = paste("genes ranked by", metric_txt),
                  y = "running enrichment score") +
    theme_omicsCore(base_size = 12) +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_text(hjust = 0),
                   axis.title.y = ggplot2::element_text(hjust = 1))
}

# A p-value for a subtitle: three digits, or two in scientific notation
# when small.
format_enrich_p <- function(p) {
  if (length(p) != 1L || !is.finite(p)) return("NA")
  if (p < 0.001) formatC(p, digits = 1, format = "e") else formatC(p, digits = 3, format = "fg")
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
    ggplot2::labs(title = "Pathways enriched", x = x_label, y = NULL) +
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
  bl <- signif_capped_breaks(lim, n = 4L)
  ggplot2::scale_color_gradient(low = omics_colors$scale_low,
                                high = omics_colors$scale_high,
                                name = p_axis_label(p_col),
                                limits = lim$limits, oob = scales::squish,
                                breaks = bl$breaks, labels = bl$labels,
                                guide = ggplot2::guide_colourbar(order = 1L))
}

# Breaks and labels for a scale capped by signif_limits(): round breaks
# below the cap, kept clear of it so the two labels do not run together,
# then the cap itself marked "at least". Left to ggplot (waiver) when
# nothing is capped.
signif_capped_breaks <- function(lim, n = 4L) {
  if (!isTRUE(lim$capped)) {
    return(list(breaks = ggplot2::waiver(), labels = ggplot2::waiver()))
  }
  lo <- lim$limits[1L]
  cap <- lim$limits[2L]
  b <- scales::breaks_extended(n = n)(lim$limits)
  b <- b[b >= lo & b <= cap - 0.25 * (cap - lo)]
  list(breaks = c(b, cap),
       labels = c(format_signif_break(b), paste0("\u2265 ", format_signif_break(cap))))
}

# The same capped significance as a point size, for a plot whose colour
# already says something else (the comparison plot's gene list, or the
# NES). Left to the data, one pathway at -log10 p = 200 took the largest
# dot and every pathway between 2 and 20 came out alike and small.
# Values past the cap take the largest size, and the key's top label
# says "at least". Three keys or so, as for the dot plot's size key:
# under a phone-width panel the key runs across, and more large dots ran
# off the edge.
signif_size_scale <- function(values, p_col, range = c(2, 6), guide = "legend") {
  lim <- signif_limits(values)
  bl <- signif_capped_breaks(lim, n = 3L)
  if (!isTRUE(lim$capped)) bl$breaks <- size_breaks
  # scale_size_continuous() takes no `oob`; this is what it builds, with
  # the squish that keeps an outlier as the largest dot instead of
  # dropping it.
  ggplot2::continuous_scale("size", palette = scales::area_pal(range),
                            name = p_axis_label(p_col),
                            limits = lim$limits, oob = scales::squish,
                            breaks = bl$breaks, labels = bl$labels, guide = guide)
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
#
# `top` is the largest ordinary value before that rounding: where the
# bar view cuts its axis. Cut at the rounded cap, a bar at 15.35 would
# be marked as broken at 15.
signif_limits <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(list(limits = NULL, capped = FALSE, top = NA_real_))
  lo <- min(x)
  hi <- max(x)
  q <- stats::quantile(x, c(0.25, 0.75), names = FALSE)
  cap <- max(x[x <= q[2L] + 1.5 * (q[2L] - q[1L])])
  if (cap == hi && length(x) >= 3L && any(x < hi)) {
    second <- max(x[x < hi])
    if (hi > 4 * second) cap <- second
  }
  capped <- cap < hi && cap > lo
  top <- if (capped) cap else hi
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
  list(limits = c(lo, cap), capped = capped, top = top)
}

plot_enrich_bar <- function(df, p_col) {
  if (ora_list_facet(df)) df <- with_list_label(df)
  df$.label <- wrap_pathway_name(df$pathway_name)
  df$.label <- factor(df$.label, levels = unique(df$.label[order(df[[p_col]], decreasing = TRUE)]))
  df$.signif <- neg_log10_p(df[[p_col]])

  # Colour means direction, as everywhere else in the app. The fill was
  # the database, in ggplot's default hue (salmon for the usual single
  # database) with its legend hidden: a colour that meant nothing.
  #  * ORA split by gene list: the list's direction. The facet strips
  #    already name the list, so no legend.
  #  * GSEA: the sign of the NES, with a legend above the panel, since
  #    nothing else on the plot says which way a bar's pathway moved.
  #  * Anything else, pooled lists included: one neutral ink colour.
  fill_by_sign <- FALSE
  if (ora_list_facet(df)) {
    df$.fill <- unname(c(up = omics_colors$up, down = omics_colors$down)[df$direction])
  } else if (any(is.finite(df$effect))) {
    df$.fill <- ifelse(df$effect >= 0, omics_colors$up, omics_colors$down)
    fill_by_sign <- TRUE
  } else {
    df$.fill <- NA_character_
  }
  df$.fill[is.na(df$.fill)] <- omics_colors$fg_dark

  # Bar length is -log10(p), and one pathway at 248 left every other bar
  # a sliver. The axis stops at the same cap as the dot plot's colour
  # scale (the largest value that is not an outlier); a bar beyond it is
  # drawn to the cap, cut by a break mark near its end, with its true
  # value printed past the end -- an axis label alone would not tell the
  # reader that this one bar is longer than drawn.
  #
  # A bar is cut only when it is more than three times the longest
  # ordinary one. A colour scale loses its spread to any outlier, but a
  # bar twice as long as the next still leaves the others readable, and
  # a cut bar is harder to read than a long one.
  lim <- signif_limits(df$.signif)
  cap <- if (lim$capped && max(df$.signif) > 3 * lim$top) lim$top else Inf
  df$.bar <- pmin(df$.signif, cap)
  over <- df[df$.signif > cap, , drop = FALSE]

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$.bar, y = .data$.label)) +
    ggplot2::geom_col(ggplot2::aes(fill = .data$.fill), width = 0.75) +
    enrich_facets(df, scales = "free_y") +
    ggplot2::labs(title = "Pathways enriched", x = p_axis_label(p_col), y = NULL) +
    theme_omics_labelled() +
    enrich_narrow_theme()

  if (fill_by_sign) {
    eff <- if (all(df$effect_type %in% "nes")) "NES" else "effect"
    labs_sign <- stats::setNames(paste(eff, c("> 0 (up)", "< 0 (down)")),
                                 c(omics_colors$up, omics_colors$down))
    shown <- intersect(names(labs_sign), df$.fill)
    p <- p + ggplot2::scale_fill_identity(guide = "legend", name = NULL,
                                          breaks = shown, labels = unname(labs_sign[shown])) +
      ggplot2::theme(legend.position = "top", legend.justification = "left")
  } else {
    p <- p + ggplot2::scale_fill_identity()
  }

  if (nrow(over) > 0L) {
    # A slanted line needs a numeric y. Each panel's y axis holds only
    # its own pathways, in level order, so a bar's position is its rank
    # among the labels of its panel, not its place among all labels.
    panel <- if (".list" %in% names(df)) paste(df$database, df$.list) else df$database
    over_panel <- panel[df$.signif > cap]
    over$.y <- vapply(seq_len(nrow(over)), function(i) {
      present <- sort(unique(as.integer(df$.label[panel == over_panel[i]])))
      match(as.integer(over$.label[i]), present)
    }, numeric(1))
    # Two short white slashes across the bar, a little before its end:
    # the usual mark for an axis that has been cut. Spaced in proportion
    # to the axis, wide enough apart that the two stay two at phone
    # width, where the panel is about 110 px.
    at <- cap * 0.86
    d <- cap * 0.025
    marks <- rbind(transform(over, .x0 = at - d, .x1 = at),
                   transform(over, .x0 = at + 1.2 * d, .x1 = at + 2.2 * d))
    p <- p +
      ggplot2::geom_segment(
        data = marks,
        ggplot2::aes(x = .data$.x0, xend = .data$.x1,
                     y = .data$.y - 0.42, yend = .data$.y + 0.42),
        color = "white", linewidth = 0.9, inherit.aes = FALSE) +
      ggplot2::geom_text(
        data = over,
        ggplot2::aes(x = cap, label = format_bar_value(.data$.signif)),
        hjust = -0.15, size = 3.2, color = omics_colors$fg_dark)
  }

  # Bars start at the axis; on the right, room for the printed value.
  # Clipping off, so a label a few pixels wider than that room on a
  # narrow panel is not cut in half.
  p + ggplot2::scale_x_continuous(
    expand = ggplot2::expansion(mult = c(0, if (nrow(over) > 0L) 0.22 else 0.05))) +
    ggplot2::coord_cartesian(clip = "off")
}

format_bar_value <- function(x) {
  as.character(ifelse(x >= 10, round(x), round(x, 1)))
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
      title = "Pathways enriched in the ranked gene list",
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
