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
#'   `database`, `direction`, thresholds).
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
  # One row per pathway and comparison: ORA run in both directions can
  # list a pathway twice for one comparison; the stronger one is kept.
  df <- df[order(df[[p_col]]), , drop = FALSE]
  df <- df[!duplicated(paste(df$comparison, df$pathway_id)), , drop = FALSE]
  # Rows ordered by how many comparisons a pathway is significant in,
  # then by its best p: the shared pathways at the top.
  sig <- df[[p_col]] < p_cutoff
  n_sig <- tapply(sig, df$pathway_id, sum)
  best <- tapply(df[[p_col]], df$pathway_id, min)
  ord <- names(sort(-n_sig * 1e6 + rank(best)[names(n_sig)]))
  labels <- stats::setNames(truncate_pathway_name(df$pathway_name), df$pathway_id)
  df$.row <- factor(df$pathway_id, levels = rev(ord))
  df$.col <- factor(gsub("_vs_", " vs ", df$comparison, fixed = TRUE),
                    levels = gsub("_vs_", " vs ", cmps, fixed = TRUE))
  df$.neglog <- -log10(pmax(df[[p_col]], .Machine$double.xmin))
  df$.sig <- ifelse(sig, "yes", "no")
  is_gsea <- any(df$effect_type %in% "nes")

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$.col, y = .data$.row)) +
    ggplot2::scale_y_discrete(labels = labels[levels(df$.row)]) +
    ggplot2::scale_shape_manual(values = c(yes = 16, no = 1),
                                name = sprintf("%s < %s",
                                               if (p_preference == "adjusted") "adj. p" else "p",
                                               format(p_cutoff))) +
    ggplot2::labs(title = "Pathways across comparisons",
                  subtitle = paste(toupper(bundle$params$type %||% ""),
                                   paste(bundle$params$database, collapse = ", ")),
                  x = NULL, y = NULL) +
    theme_omics_labelled() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
  if (is_gsea) {
    p + ggplot2::geom_point(ggplot2::aes(size = .data$.neglog, color = .data$effect,
                                         shape = .data$.sig)) +
      ggplot2::scale_color_gradient2(low = omics_colors$down, mid = "#F2F2F2",
                                     high = omics_colors$up, midpoint = 0,
                                     name = "NES") +
      ggplot2::scale_size_continuous(name = "-log10 p", range = c(2, 7))
  } else {
    p + ggplot2::geom_point(ggplot2::aes(size = .data$overlap_size,
                                         color = .data$.neglog,
                                         shape = .data$.sig)) +
      ggplot2::scale_color_gradient(low = "#9DB2D9", high = omics_colors$up,
                                    name = "-log10 p") +
      ggplot2::scale_size_continuous(name = "genes", range = c(2, 7))
  }
}
