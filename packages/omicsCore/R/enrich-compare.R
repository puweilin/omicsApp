# Enrichment of several comparisons, side by side.
#
# With one control and several treatments the question is rarely "which
# pathways does TreatA hit" alone but "which pathways do the treatments
# share, and which are specific to one" -- clusterProfiler's
# compareCluster() question. Each comparison is enriched on its own, with
# the same settings, and the tables are stacked under their `comparison`
# labels.

#' Enrich every comparison of a differential run
#'
#' @param diff_bundle An [`analysis_bundle`][is_analysis_bundle()] from
#'   [run_diff()], usually holding several comparisons.
#' @param comparisons Which comparisons; default all of
#'   [diff_comparisons()].
#' @param ... Passed to [run_enrichment()] for every comparison (`type`,
#'   `database`, `direction`, thresholds). As there, ORA tests the up- and
#'   down-regulated hits as separate lists unless `direction` says
#'   otherwise, and each row's `direction` says which list it came from.
#' @return An [`analysis_bundle`][is_analysis_bundle()] named
#'   `"compare_enrichment"` whose `results$enrich_result_df` stacks the
#'   per-comparison tables (the `comparison` column says which), and whose
#'   `results$per_comparison` keeps each comparison's own bundle.
#'   A comparison whose enrichment fails (no hits at the thresholds, say)
#'   contributes no rows and a line in `warnings`.
#' @export
#' @family enrich
compare_enrichment <- function(diff_bundle, comparisons = NULL, ...) {
  assert_diff_bundle(diff_bundle)
  assert_character(comparisons, "comparisons", allow_null = TRUE)
  available <- diff_comparisons(diff_bundle)
  comparisons <- comparisons %||% available
  unknown <- setdiff(comparisons, available)
  if (length(unknown)) {
    stop("Not in the differential result: ", paste(unknown, collapse = ", "),
         ". Available: ", paste(available, collapse = ", "), ".", call. = FALSE)
  }
  warns <- character(0)
  per <- list()
  for (cmp in comparisons) {
    one <- if (length(available) > 1L) select_comparison(diff_bundle, cmp) else diff_bundle
    res <- tryCatch(run_enrichment(one, ...), error = function(e) e)
    if (inherits(res, "error")) {
      warns <- c(warns, sprintf("%s: %s", cmp, conditionMessage(res)))
      next
    }
    # A note from the run itself (gene names matched ignoring case) is
    # about this comparison's genes as much as its result.
    if (length(res$warnings)) {
      warns <- c(warns, sprintf("%s: %s", cmp, res$warnings))
    }
    df <- res$results$enrich_result_df
    if (nrow(df)) df$comparison <- cmp
    res$results$enrich_result_df <- df
    if (!nrow(df)) {
      warns <- c(warns, sprintf("%s: no pathway passed the thresholds.", cmp))
    }
    per[[cmp]] <- res
  }
  if (!length(per)) {
    stop("No comparison could be enriched: ", paste(warns, collapse = "; "),
         call. = FALSE)
  }
  df <- dplyr::bind_rows(lapply(per, function(b) b$results$enrich_result_df))
  check_enrich_result_schema(df)
  first <- per[[1L]]$params
  new_analysis_bundle(
    analysis_name = "compare_enrichment",
    input_info = diff_bundle$input_info,
    params = c(first[setdiff(names(first), "comparison")],
               list(comparison = names(per))),
    results = list(enrich_result_df = df, per_comparison = per),
    warnings = warns
  )
}

#' Pathways across comparisons (a compareCluster-style dot plot)
#'
#' One column per comparison, one row per pathway: the pathways shown are
#' the `top_n` strongest of each comparison, so a pathway that matters in
#' only one of them still appears -- and its empty cells in the other
#' columns are the finding.
#'
#' Point size is -log10 of the chosen p-value, capped as [plot_enrichment()]
#' caps its colour: a pathway beyond the cap takes the largest dot and the
#' key's top label reads "at least" the cap. Colour is the gene list a
#' pathway was found among (ORA run on up and down separately) or the NES
#' (GSEA). Pathways that miss `p_cutoff` are drawn hollow, and only then
#' does the plot carry a key for filled and hollow. When every comparison
#' has the same control the columns are named by the treatment and the
#' control is named once, in the subtitle.
#'
#' @param bundle A bundle from [compare_enrichment()].
#' @param top_n Pathways taken from each comparison.
#' @param p_cutoff Pathways above this (adjusted, unless `p_preference`
#'   says raw) p-value are drawn hollow.
#' @param p_preference `"adjusted"` or `"raw"`.
#' @return A `ggplot` object.
#' @export
#' @family enrich
plot_enrichment_comparison <- function(bundle, top_n = 8L, p_cutoff = 0.05,
                                       p_preference = c("adjusted", "raw")) {
  if (!is_analysis_bundle(bundle) ||
      !identical(bundle$analysis_name, "compare_enrichment")) {
    stop("`bundle` must come from compare_enrichment().", call. = FALSE)
  }
  assert_count(top_n, "top_n")
  assert_number(p_cutoff, "p_cutoff", lower = 0, upper = 1)
  p_preference <- match.arg(p_preference)
  p_col <- if (p_preference == "adjusted") "adj_p_value" else "p_value"
  df <- bundle$results$enrich_result_df
  df <- df[!is.na(df[[p_col]]), , drop = FALSE]
  if (!nrow(df)) {
    return(ggplot2::ggplot() + ggplot2::theme_void() +
             ggplot2::annotate("text", x = 0.5, y = 0.5,
                               label = "No enriched pathways to compare.") +
             ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1))
  }
  cmps <- unique(bundle$params$comparison %||% df$comparison)
  pick <- unique(unlist(lapply(cmps, function(cmp) {
    sub <- df[df$comparison == cmp, , drop = FALSE]
    utils::head(sub$pathway_id[order(sub[[p_col]])], top_n)
  })))
  df <- df[df$pathway_id %in% pick, , drop = FALSE]
  is_gsea <- any(df$effect_type %in% "nes")
  # ORA with up and down tested separately: which list found the pathway
  # is the colour, since a pathway enriched among the up-regulated genes
  # in one comparison and the down-regulated in another is the opposite
  # of shared. Read from the rows that pass `p_cutoff` when there are
  # any, so a list that found the pathway only weakly does not turn a
  # clear "up" into "up and down".
  directional_ora <- !is_gsea && any(df$direction %in% c("up", "down"))
  key <- paste(df$comparison, df$pathway_id, sep = "\r")
  if (directional_ora) {
    sig_row <- df[[p_col]] < p_cutoff
    lists <- vapply(split(seq_len(nrow(df)), key), function(i) {
      j <- if (any(sig_row[i])) i[sig_row[i]] else i
      d <- unique(stats::na.omit(df$direction[j]))
      if (all(c("up", "down") %in% d)) "Up and down"
      else if ("up" %in% d) "Up-regulated genes"
      else if ("down" %in% d) "Down-regulated genes"
      else "Pooled"
    }, character(1))
    df$.list <- unname(lists[key])
  }
  # One row per pathway and comparison: ORA run in both directions can
  # list a pathway twice for one comparison; the stronger one is kept.
  ord <- order(df[[p_col]])
  df <- df[ord, , drop = FALSE]
  key <- key[ord]
  df <- df[!duplicated(key), , drop = FALSE]
  # Rows ordered by how many comparisons a pathway is significant in,
  # then by its best p: the shared pathways at the top.
  sig <- df[[p_col]] < p_cutoff
  n_sig <- tapply(sig, df$pathway_id, sum)
  best <- tapply(df[[p_col]], df$pathway_id, min)
  ord <- names(sort(-n_sig * 1e6 + rank(best)[names(n_sig)]))
  labels <- stats::setNames(truncate_pathway_name(df$pathway_name), df$pathway_id)
  df$.row <- factor(df$pathway_id, levels = rev(ord))
  cols <- comparison_columns(cmps)
  df$.col <- factor(unname(cols$labels[df$comparison]), levels = unique(unname(cols$labels)))
  df$.neglog <- neg_log10_p(df[[p_col]])
  df$.sig <- ifelse(sig, "yes", "no")

  # Filled or hollow says whether a dot passes `p_cutoff`. With every dot
  # filled, a key of one entry ("adjusted p < 0.05: yes") told the reader
  # nothing and took a legend's room; the key appears only when there is
  # a hollow dot to explain, and then in words.
  p_name <- if (p_preference == "adjusted") "adjusted p" else "p"
  shape_scale <- ggplot2::scale_shape_manual(
    values = c(yes = 16, no = 1), breaks = c("yes", "no"),
    labels = c(yes = sprintf("%s < %s", p_name, format(p_cutoff)),
               no = "not significant"),
    name = NULL,
    guide = if (any(!sig)) {
      ggplot2::guide_legend(order = 3L, override.aes = list(size = 3, colour = omics_colors$fg_dark))
    } else "none")

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$.col, y = .data$.row)) +
    ggplot2::scale_y_discrete(labels = labels[levels(df$.row)]) +
    shape_scale +
    ggplot2::labs(title = "Pathways across comparisons",
                  subtitle = comparison_subtitle(bundle$params, cols$control),
                  x = NULL, y = NULL) +
    theme_omics_labelled() +
    enrich_narrow_theme() +
    comparison_legend_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1))
  # Size is the significance in every case but the pooled ORA, where the
  # colour carries it: the same -log10 of the p the plot is told to use,
  # capped as the dot plot caps its colour, so one extreme pathway does
  # not leave every other dot the same small size.
  size_scale <- signif_size_scale(df$.neglog, p_col, guide = ggplot2::guide_legend(order = 2L))
  if (is_gsea) {
    p + ggplot2::geom_point(ggplot2::aes(size = .data$.neglog, color = .data$effect,
                                         shape = .data$.sig)) +
      ggplot2::scale_color_gradient2(low = omics_colors$down, mid = "#F2F2F2",
                                     high = omics_colors$up, midpoint = 0,
                                     name = "NES",
                                     guide = ggplot2::guide_colourbar(order = 1L)) +
      size_scale
  } else if (directional_ora) {
    lists <- c(`Up-regulated genes` = omics_colors$up,
               `Down-regulated genes` = omics_colors$down,
               `Up and down` = omics_colors$conc_down_up,
               Pooled = omics_colors$ns)
    shown_lists <- intersect(names(lists), df$.list)
    p + ggplot2::geom_point(ggplot2::aes(size = .data$.neglog, color = .data$.list,
                                         shape = .data$.sig)) +
      # Up first, as everywhere else. One entry per line wherever the key
      # is: side by side under a 290 px phone panel, "Down-regulated
      # genes" ran off the right edge.
      ggplot2::scale_color_manual(
        values = lists, breaks = shown_lists,
        name = "found among",
        guide = ggplot2::guide_legend(order = 1L, ncol = 1L,
                                      override.aes = list(size = 3.5))) +
      size_scale
  } else {
    p + ggplot2::geom_point(ggplot2::aes(size = .data$overlap_size,
                                         color = .data$.neglog,
                                         shape = .data$.sig)) +
      signif_colour_scale(df$.neglog, p_col) +
      ggplot2::scale_size_continuous(name = "genes in list", range = c(2, 6),
                                     breaks = size_breaks,
                                     guide = ggplot2::guide_legend(order = 2L))
  }
}

# What was run, under the title. The full line ("ORA hallmark - up- and
# down-regulated genes tested separately") is wider than a phone-width
# plot and, at report size, ran into the legend beside the panel; split
# at its dot it is two short lines that fit either.
comparison_subtitle <- function(params, control = NULL) {
  method <- trimws(paste(toupper(params$type %||% ""),
                         paste(params$database, collapse = ", ")))
  if (!is.null(control)) method <- paste0(method, " \u00B7 each vs ", control)
  # Shorter than ora_list_caption(), whose wording filled the whole width
  # of a 310 px phone plot; the legend names the two lists anyway.
  separate <- identical(params$type, "ora") && identical(params$direction, "separate")
  lists <- if (separate) "up and down genes tested separately"
           else sub("^ \u00B7 ", "", ora_list_caption(params))
  lines <- c(method, if (nzchar(lists)) lists)
  paste(wrap_label(lines, width = 48L, max_lines = 2L), collapse = "\n")
}

# Column labels. Several treatments against one control is the usual
# design, and "vs Vehicle control" repeated under every column took most
# of each label: wrapped over two lines and slanted, five long ones ran
# into each other. When every comparison shares its control, the
# columns are named by the treatment alone (shortened past 24
# characters, on one line) and the control is said once, in the
# subtitle. Otherwise each column keeps its full "A vs B", wrapped.
comparison_columns <- function(cmps) {
  parts <- strsplit(cmps, "_vs_", fixed = TRUE)
  two <- all(lengths(parts) == 2L)
  control <- if (two) unique(vapply(parts, `[`, character(1), 2L))
  if (length(cmps) > 1L && length(control) == 1L) {
    cases <- vapply(parts, `[`, character(1), 1L)
    short <- truncate_pathway_name(cases, max_chars = 24L)
    # Two treatments that only differ past the cut keep their full names.
    if (anyDuplicated(short)) short <- cases
    return(list(labels = stats::setNames(short, cmps), control = control))
  }
  list(labels = stats::setNames(wrap_comparison(cmps, width = 18L), cmps), control = NULL)
}

# The legends sit at the top of the space beside the panel, under the
# subtitle, not centred on the panel's height: centred, a tall legend
# stack reached up into the subtitle. (Below the panel -- where the app
# puts them on a phone -- enrich_narrow_theme() starts the stack at the
# left, and each legend in it starts there too rather than being centred
# under the widest.)
comparison_legend_theme <- function() {
  th <- ggplot2::theme(legend.box.just = "left")
  if (!"legend.justification.right" %in% names(ggplot2::get_element_tree())) {
    return(th + ggplot2::theme(legend.justification = "top"))
  }
  th + ggplot2::theme(legend.justification.right = "top")
}
