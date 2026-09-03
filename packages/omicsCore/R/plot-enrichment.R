#' Enrichment summary plot
#'
#' Visualises the standardized enrichment table produced by
#' [run_enrichment()]. The `"dot"` view draws a dotplot of the top `top_n`
#' pathways per database (size = overlap, color = adjusted p-value); the
#' `"bar"` view draws a horizontal bar chart of the same selection. For
#' GSEA bundles a `"gsea_dot"` view is available that splits pathways by
#' direction (up / down).
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
    return(empty_enrich_plot("No pathways to plot."))
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
  if (!requireNamespace("enrichplot", quietly = TRUE)) {
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
  picked <- lapply(split_df, function(x) utils::head(x, top_n))
  out <- do.call(rbind, picked)
  rownames(out) <- NULL
  out
}

empty_enrich_plot <- function(label) {
  ggplot2::ggplot() +
    ggplot2::theme_void() +
    ggplot2::annotate("text", x = 0.5, y = 0.5, label = label,
                      color = "#4D4D4D", size = 4) +
    ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1)
}

plot_enrich_dot <- function(df, p_col) {
  df$.label <- truncate_pathway_name(df$pathway_name)
  df$.label <- factor(df$.label, levels = unique(df$.label[order(-df[[p_col]])]))

  # Effect (NES for GSEA, log2 odds for ORA) says which direction a
  # pathway moved; overlap size only says how many genes were in it.
  # Prefer the former when the result carries it, and fall back
  # otherwise so an ORA result without an effect column still plots.
  has_effect <- "effect" %in% names(df) && any(is.finite(df$effect))
  x_aes <- if (has_effect) "effect" else "overlap_size"
  df$.signif <- -log10(pmax(df[[p_col]], .Machine$double.xmin))

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data[[x_aes]], y = .data$.label,
                 size = .data$overlap_size, color = .data$.signif)
  ) +
    ggplot2::geom_point(na.rm = TRUE) +
    ggplot2::facet_wrap(~ .data$database, scales = "free", ncol = 1) +
    # Low to high significance, matching the volcano: a reader who has
    # learnt one colour scale reads the other without relearning it.
    ggplot2::scale_color_gradient(low = omics_colors$scale_low,
                                  high = omics_colors$scale_high,
                                  name = paste0("-log10(", p_col, ")")) +
    ggplot2::scale_size_continuous(range = c(2, 6), name = "overlap") +
    ggplot2::labs(
      title = "Enrichment",
      x = if (has_effect) "effect" else "overlap size",
      y = NULL
    ) +
    theme_omicsCore()

  if (has_effect) {
    p <- p + ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                                 color = omics_colors$ns)
  }
  p
}

plot_enrich_bar <- function(df, p_col) {
  df$.label <- truncate_pathway_name(df$pathway_name)
  df$.label <- factor(df$.label, levels = unique(df$.label[order(df[[p_col]], decreasing = TRUE)]))
  ggplot2::ggplot(
    df,
    ggplot2::aes(x = -log10(pmax(.data[[p_col]], .Machine$double.xmin)),
                 y = .data$.label, fill = .data$database)
  ) +
    ggplot2::geom_col() +
    ggplot2::facet_wrap(~ .data$database, scales = "free_y", ncol = 1) +
    ggplot2::guides(fill = "none") +
    ggplot2::labs(
      title = "Enrichment",
      x = paste0("-log10(", p_col, ")"),
      y = NULL
    ) +
    theme_omicsCore()
}

plot_enrich_gsea_dot <- function(df, p_col) {
  if (!"direction" %in% colnames(df) || all(is.na(df$direction))) {
    stop("`gsea_dot` view requires a `direction` column with non-NA values.")
  }
  df$.label <- truncate_pathway_name(df$pathway_name)
  df$.label <- factor(df$.label, levels = unique(df$.label[order(-df[[p_col]])]))
  ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data$effect, y = .data$.label,
                 size = .data$gene_set_size, color = .data[[p_col]])
  ) +
    ggplot2::geom_point(na.rm = TRUE) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = omics_colors$ns) +
    ggplot2::facet_grid(rows = ggplot2::vars(.data$database),
                        cols = ggplot2::vars(.data$direction),
                        scales = "free_y") +
    ggplot2::scale_color_gradient(low = omics_colors$up, high = omics_colors$ns, name = p_col) +
    ggplot2::scale_size_continuous(range = c(2, 6), name = "set size") +
    ggplot2::labs(
      title = "GSEA",
      x = "NES",
      y = NULL
    ) +
    theme_omicsCore()
}

# 45, not 60. Hallmark names are short enough that 60 never bit, but a
# GO BP or Reactome term routinely runs past it, and at 60 characters
# the y-axis labels take more of the panel than the points do. 45 keeps
# a term identifiable -- the discriminating words are at the front --
# without the plot becoming a caption with dots attached.
truncate_pathway_name <- function(x, max_chars = 45L) {
  x <- as.character(x)
  too_long <- !is.na(x) & nchar(x) > max_chars
  x[too_long] <- paste0(substr(x[too_long], 1L, max_chars - 1L), "\u2026")
  x
}
