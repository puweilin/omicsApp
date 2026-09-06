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
  df$.label <- wrap_pathway_name(df$pathway_name)
  df$.label <- factor(df$.label, levels = unique(df$.label[order(-df[[p_col]])]))

  # Effect (NES for GSEA, log2 odds for ORA) says which direction a
  # pathway moved; overlap size only says how many genes were in it.
  # Prefer the former when the result carries it, and fall back
  # otherwise so an ORA result without an effect column still plots.
  usable <- function(col) col %in% names(df) && any(is.finite(df[[col]]))

  has_effect <- usable("effect")
  x_aes <- if (has_effect) "effect" else "overlap_size"
  df$.signif <- -log10(pmax(df[[p_col]], .Machine$double.xmin))

  # The size aesthetic needs the same fallback as x, and for the same
  # kind of result. GSEA has no overlap: standardize_enrich_results only
  # fills overlap_size from `Count` or `GeneRatio`, which fgsea emits
  # neither of, so the column is entirely NA. Mapping size to it made
  # every point NA-sized, and geom_point drops those -- the panel drew
  # its axes and facet strips and not one dot.
  #
  # That reads as "enrichment found nothing", which is the same thing an
  # empty result looks like, so the failure hides as a result.
  size_aes <- if (usable("overlap_size")) "overlap_size" else "gene_set_size"
  has_size <- usable(size_aes)

  aes_args <- list(x = quote(.data[[x_aes]]), y = quote(.data$.label),
                   color = quote(.data$.signif))
  if (has_size) aes_args$size <- quote(.data[[size_aes]])

  p <- ggplot2::ggplot(df, do.call(ggplot2::aes, aes_args)) +
    ggplot2::geom_point(na.rm = TRUE) +
    ggplot2::facet_wrap(~ .data$database, scales = "free", ncol = 1) +
    # Low to high significance, matching the volcano: a reader who has
    # learnt one colour scale reads the other without relearning it.
    ggplot2::scale_color_gradient(low = omics_colors$scale_low,
                                  high = omics_colors$scale_high,
                                  name = paste0("-log10(", p_col, ")")) +
    ggplot2::labs(
      title = "Enrichment",
      x = if (has_effect) "effect" else "overlap size",
      y = NULL
    ) +
    theme_omics_labelled()

  if (has_size) {
    # Named for the column actually mapped, so a GSEA panel does not
    # label a gene-set size "overlap".
    p <- p + ggplot2::scale_size_continuous(
      range = c(2, 6),
      name = if (identical(size_aes, "overlap_size")) "overlap" else "set size")
  }

  if (has_effect) {
    p <- p + ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                                 color = omics_colors$ns)
  }
  p
}

plot_enrich_bar <- function(df, p_col) {
  df$.label <- wrap_pathway_name(df$pathway_name)
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
    theme_omics_labelled()
}

plot_enrich_gsea_dot <- function(df, p_col) {
  if (!"direction" %in% colnames(df) || all(is.na(df$direction))) {
    stop("`gsea_dot` view requires a `direction` column with non-NA values.")
  }
  df$.label <- wrap_pathway_name(df$pathway_name)
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
    theme_omics_labelled()
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
  x <- as.character(x)
  vapply(x, function(s) {
    if (is.na(s) || !nzchar(s)) return(s)
    lines <- strwrap(s, width = width)
    if (length(lines) > max_lines) {
      lines <- lines[seq_len(max_lines)]
      lines[max_lines] <- paste0(lines[max_lines], "\u2026")
    }
    paste(lines, collapse = "\n")
  }, character(1), USE.NAMES = FALSE)
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
