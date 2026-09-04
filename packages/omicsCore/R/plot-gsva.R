#' GSVA heatmap
#'
#' Draws a pathways-by-samples heatmap from a [run_gsva()] bundle. When
#' `ComplexHeatmap` is installed the function returns a `Heatmap` object
#' (recommended); otherwise it falls back to a ggplot2 tile plot so the
#' function is always usable from a clean install.
#'
#' @param bundle An [`analysis_bundle`][is_analysis_bundle()] produced by
#'   [run_gsva()].
#' @param top_n Number of pathways to plot when `pathways` is `NULL`.
#'   Selected by variance across samples.
#' @param pathways Optional character vector of pathway names to plot.
#' @param meta_df Optional sample metadata `data.frame` (or named list of
#'   columns) used for column annotations. Rows must be indexable by sample
#'   id (column names of the score matrix).
#' @param annotation_cols Optional character vector of columns from
#'   `meta_df` to surface as a `HeatmapAnnotation`.
#' @param scale One of `"row"` (default), `"none"`, or `"column"`.
#' @param cluster_rows,cluster_cols Whether to cluster rows / columns.
#' @param title Heatmap title.
#'
#' @return A `ComplexHeatmap::Heatmap` object when `ComplexHeatmap` is
#'   available, otherwise a `ggplot` tile plot.
#' @export
#' @family enrich
plot_gsva_heatmap <- function(
  bundle,
  top_n = 25L,
  pathways = NULL,
  meta_df = NULL,
  annotation_cols = NULL,
  scale = c("row", "none", "column"),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  title = "GSVA"
) {
  if (!is_analysis_bundle(bundle) ||
      !identical(bundle$analysis_name, "run_gsva")) {
    stop("`bundle` must be an analysis_bundle from run_gsva().")
  }
  scale <- match.arg(scale)
  mat <- bundle$results$gsva_matrix
  if (is.null(mat) || !is.matrix(mat) || nrow(mat) == 0L) {
    stop("Bundle is missing a non-empty `results$gsva_matrix`.")
  }

  if (!is.null(pathways)) {
    selected <- intersect(pathways, rownames(mat))
    if (length(selected) == 0L) {
      stop("None of the requested pathways were found in the GSVA matrix.")
    }
  } else {
    row_var <- apply(mat, 1L, stats::var, na.rm = TRUE)
    selected <- rownames(mat)[order(row_var, decreasing = TRUE)]
    selected <- utils::head(selected, top_n)
  }
  mat <- mat[selected, , drop = FALSE]

  mat <- switch(scale,
    row    = t(scale(t(mat))),
    column = scale(mat),
    none   = mat
  )

  display_rownames <- tidy_pathway_names(rownames(mat))
  rownames(mat) <- display_rownames

  if (is_installed("ComplexHeatmap")) {
    plot_gsva_complexheatmap(
      mat = mat,
      meta_df = meta_df,
      annotation_cols = annotation_cols,
      cluster_rows = cluster_rows,
      cluster_cols = cluster_cols,
      title = title
    )
  } else {
    plot_gsva_ggplot(mat = mat, title = title)
  }
}

# ---- internal helpers --------------------------------------------------

tidy_pathway_names <- function(x) {
  # Strip MSigDB collection prefix (e.g. "HALLMARK_") and prettify spacing.
  x <- gsub("^(HALLMARK|KEGG|REACTOME|WP|GOBP|GOMF|GOCC|GO|BIOCARTA|PID)_",
            "", x)
  x <- gsub("_", " ", x)
  stringr::str_to_title(x)
}

plot_gsva_complexheatmap <- function(mat, meta_df, annotation_cols,
                                     cluster_rows, cluster_cols, title) {
  top_anno <- NULL
  if (!is.null(meta_df) && !is.null(annotation_cols)) {
    meta_df <- as.data.frame(meta_df)
    cols <- intersect(annotation_cols, colnames(meta_df))
    if (length(cols) > 0L) {
      anno_df <- meta_df[colnames(mat), cols, drop = FALSE]
      top_anno <- ComplexHeatmap::HeatmapAnnotation(df = anno_df)
    }
  }

  col_fn <- if (is_installed("circlize")) {
    circlize::colorRamp2(c(-2, 0, 2), c("#3C5488B2", "white", "#E64B35B2"))
  } else {
    NULL
  }

  ComplexHeatmap::Heatmap(
    mat,
    name = "Z-score",
    top_annotation = top_anno,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_cols,
    show_row_names = nrow(mat) <= 60L,
    show_column_names = TRUE,
    column_title = title,
    col = col_fn
  )
}

plot_gsva_ggplot <- function(mat, title) {
  long <- data.frame(
    pathway = rep(rownames(mat), times = ncol(mat)),
    sample  = rep(colnames(mat), each = nrow(mat)),
    score   = as.numeric(mat),
    stringsAsFactors = FALSE
  )
  long$pathway <- factor(long$pathway, levels = rev(rownames(mat)))
  long$sample <- factor(long$sample, levels = colnames(mat))
  ggplot2::ggplot(
    long,
    ggplot2::aes(x = .data$sample, y = .data$pathway, fill = .data$score)
  ) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(
      low = "#3C5488B2", mid = "white", high = "#E64B35B2",
      midpoint = 0, name = "score"
    ) +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    theme_omicsCore() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
}
