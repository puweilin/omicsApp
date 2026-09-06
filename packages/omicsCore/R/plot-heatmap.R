# Top-features heatmap. The function is overloaded:
#   * input is an `omics_input` -> heatmap of `top_n` most-variable features
#     across all samples;
#   * input is a `run_diff` analysis_bundle -> heatmap of the top features by
#     (adjusted) p-value, drawn from the original `omics_input` carried in
#     `bundle$input_info` ... except diff bundles don't carry the expression
#     matrix. So callers must also pass `input` for the bundle case.
#
# Renders via ComplexHeatmap when installed, otherwise a ggplot tile plot.

#' Top-features expression heatmap
#'
#' Two dispatch modes:
#'
#' * **`omics_input` mode** -- pass an `omics_input` as the first argument.
#'   The function selects the top `n_top` features by row variance (after
#'   `coerce_to_continuous()`) and draws a samples-by-features heatmap.
#' * **`analysis_bundle` mode** -- pass a [`run_diff()`] bundle and the
#'   corresponding `omics_input` via the `input` argument. Features are
#'   selected by ascending adjusted p-value.
#'
#' When `ComplexHeatmap` is installed the function returns a
#' `ComplexHeatmap::Heatmap`; otherwise it falls back to a ggplot2 tile
#' plot so the function is always usable from a clean install.
#'
#' @param x Either an [`omics_input`][omics_input()] or an
#'   [`analysis_bundle`][is_analysis_bundle()] from [run_diff()].
#' @param input Required when `x` is an analysis_bundle: the
#'   `omics_input` whose `expr_mat` should be used for the tiles.
#' @param n_top Number of features to display.
#' @param features Optional character vector of `feature_id`s to force.
#'   Overrides `n_top`.
#' @param scale One of `"row"` (default), `"none"`, or `"column"`.
#' @param annotation_cols Optional character vector of columns from
#'   `input$meta_df` to surface as a column annotation.
#' @param cluster_rows,cluster_cols Whether to cluster rows / columns.
#' @param show_rownames If `NULL`, auto-decide based on row count.
#' @param title Plot title.
#'
#' @return A `ComplexHeatmap::Heatmap` object when `ComplexHeatmap` is
#'   available, otherwise a `ggplot` tile plot.
#' @export
#' @family diff
#' @examples
#' \dontrun{
#'   # From an omics_input
#'   plot_heatmap(input, n_top = 30)
#'
#'   # From a diff bundle
#'   b <- run_diff(input, ...)
#'   plot_heatmap(b, input = input, n_top = 50)
#' }
plot_heatmap <- function(
  x,
  input = NULL,
  n_top = 50L,
  features = NULL,
  scale = c("row", "none", "column"),
  annotation_cols = NULL,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = NULL,
  title = NULL
) {
  if (!inherits(x, "omics_input") && !is_analysis_bundle(x)) {
    arg_stop("x", "an `omics_input` or an analysis_bundle from run_diff()", x)
  }
  assert_count(n_top, "n_top", lower = 1L)
  assert_character(features, "features", allow_null = TRUE)
  assert_character(annotation_cols, "annotation_cols", allow_null = TRUE)
  assert_flag(cluster_rows, "cluster_rows")
  assert_flag(cluster_cols, "cluster_cols")
  assert_flag(show_rownames, "show_rownames", allow_null = TRUE)
  assert_string(title, "title", allow_null = TRUE, allow_empty = TRUE)
  scale <- match.arg(scale)
  sel <- resolve_heatmap_selection(x, input = input, n_top = n_top,
                                    features = features)
  mat <- sel$mat
  meta <- sel$meta
  title <- title %||% sel$default_title

  mat <- switch(scale,
    row    = scale_rows(mat),
    column = scale(mat),
    none   = mat
  )

  if (is.null(show_rownames)) {
    show_rownames <- nrow(mat) <= 80L
  }

  if (is_installed("ComplexHeatmap")) {
    plot_heatmap_complex(
      mat = mat, meta = meta, annotation_cols = annotation_cols,
      cluster_rows = cluster_rows, cluster_cols = cluster_cols,
      show_rownames = show_rownames, title = title
    )
  } else {
    plot_heatmap_ggplot(mat = mat, title = title)
  }
}

# ---- internal helpers --------------------------------------------------

resolve_heatmap_selection <- function(x, input, n_top, features) {
  if (inherits(x, "omics_input")) {
    selected <- pick_features_by_variance(x$expr_mat, x$assay_type,
                                          n_top = n_top, features = features)
    mat <- coerce_to_continuous(x$expr_mat[selected, , drop = FALSE], x$assay_type)
    rownames(mat) <- map_feature_symbols(x$feature_df, selected)
    list(
      mat = mat,
      meta = x$meta_df,
      default_title = paste0("Top ", nrow(mat), " variable features")
    )
  } else if (is_analysis_bundle(x) &&
             identical(x$analysis_name, "run_diff")) {
    if (is.null(input) || !inherits(input, "omics_input")) {
      stop("`input` (an omics_input) is required when `x` is a diff bundle.")
    }
    res <- x$results$diff_result_df
    if (is.null(res) || nrow(res) == 0L) {
      stop("Diff bundle does not contain a non-empty `results$diff_result_df`.")
    }
    selected <- pick_features_by_diff(res, n_top = n_top, features = features)
    keep <- intersect(selected, rownames(input$expr_mat))
    if (length(keep) == 0L) {
      stop("None of the top features from the diff bundle are present in `input$expr_mat`.")
    }
    mat <- coerce_to_continuous(input$expr_mat[keep, , drop = FALSE], input$assay_type)
    rownames(mat) <- map_feature_symbols(input$feature_df, keep)
    list(
      mat = mat,
      meta = input$meta_df,
      default_title = paste0("Top ", length(keep), " features by adj. p")
    )
  } else {
    stop("`x` must be an `omics_input` or a `run_diff` analysis_bundle.")
  }
}

pick_features_by_variance <- function(mat, assay_type, n_top, features) {
  if (!is.null(features)) {
    selected <- intersect(features, rownames(mat))
    if (length(selected) == 0L) {
      stop("None of the requested features were found in `expr_mat`.")
    }
    return(selected)
  }
  mat_c <- coerce_to_continuous(mat, assay_type)
  row_var <- apply(mat_c, 1L, stats::var, na.rm = TRUE)
  row_var[is.na(row_var)] <- -Inf
  order_idx <- order(row_var, decreasing = TRUE)
  utils::head(rownames(mat)[order_idx], n_top)
}

pick_features_by_diff <- function(res, n_top, features) {
  if (!is.null(features)) {
    sel <- intersect(features, res$feature_id)
    if (length(sel) == 0L) {
      stop("None of the requested features were found in the diff result table.")
    }
    return(sel)
  }
  ord <- order(res$adj_p_value, res$p_value, na.last = NA)
  utils::head(res$feature_id[ord], n_top)
}

map_feature_symbols <- function(feature_df, feature_ids) {
  if (!is.data.frame(feature_df) ||
      !"feature_id" %in% colnames(feature_df) ||
      !"feature_symbol" %in% colnames(feature_df)) {
    return(feature_ids)
  }
  idx <- match(feature_ids, feature_df$feature_id)
  sym <- feature_df$feature_symbol[idx]
  ifelse(is.na(sym) | !nzchar(sym), feature_ids, sym)
}

scale_rows <- function(mat) {
  out <- t(scale(t(mat)))
  out[!is.finite(out)] <- 0
  out
}

plot_heatmap_complex <- function(mat, meta, annotation_cols, cluster_rows,
                                  cluster_cols, show_rownames, title) {
  top_anno <- NULL
  if (!is.null(meta) && !is.null(annotation_cols)) {
    cols <- intersect(annotation_cols, colnames(meta))
    if (length(cols) > 0L) {
      anno_df <- as.data.frame(meta)[colnames(mat), cols, drop = FALSE]
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
    show_row_names = show_rownames,
    show_column_names = TRUE,
    column_title = title,
    col = col_fn
  )
}

plot_heatmap_ggplot <- function(mat, title) {
  long <- data.frame(
    feature = rep(rownames(mat), times = ncol(mat)),
    sample  = rep(colnames(mat), each = nrow(mat)),
    value   = as.numeric(mat),
    stringsAsFactors = FALSE
  )
  long$feature <- factor(long$feature, levels = rev(rownames(mat)))
  long$sample <- factor(long$sample, levels = colnames(mat))
  ggplot2::ggplot(
    long,
    ggplot2::aes(x = .data$sample, y = .data$feature, fill = .data$value)
  ) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(
      low = "#3C5488B2", mid = "white", high = "#E64B35B2",
      midpoint = 0, name = "Z-score"
    ) +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    theme_omicsCore() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
}
