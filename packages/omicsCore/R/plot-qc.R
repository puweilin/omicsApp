#' QC visualizations
#'
#' Builds standard ggplot panels from a [run_qc()] bundle.
#'
#' @param bundle An `analysis_bundle` produced by [run_qc()].
#' @param view One of `"missing"`, `"pca"`, `"connectivity"`, `"imputation"`.
#' @param color_by Optional name of a column in the cleaned input's `meta_df`
#'   used to color samples in the `"pca"` view.
#' @param ... Reserved for future arguments.
#'
#' @return A `ggplot` object.
#' @export
#' @family qc
plot_qc <- function(bundle,
                    view = c("missing", "pca", "connectivity", "imputation"),
                    color_by = NULL,
                    ...) {
  if (!is_analysis_bundle(bundle) || !identical(bundle$analysis_name, "run_qc")) {
    stop("`bundle` must be an analysis_bundle from run_qc().")
  }
  view <- match.arg(view)

  switch(view,
    missing      = plot_qc_missing(bundle),
    pca          = plot_qc_pca(bundle, color_by = color_by),
    connectivity = plot_qc_connectivity(bundle),
    imputation   = plot_qc_imputation(bundle)
  )
}

# ---- views -------------------------------------------------------------

plot_qc_missing <- function(bundle) {
  miss <- bundle$results$qc_summary$missingness
  sample_df <- miss$sample_metrics
  feature_df <- miss$feature_metrics

  combined <- rbind(
    data.frame(
      axis = "sample",
      id = sample_df$sample_id,
      missing_rate = sample_df$missing_rate,
      stringsAsFactors = FALSE
    ),
    data.frame(
      axis = "feature",
      id = feature_df$feature_id,
      missing_rate = feature_df$missing_rate,
      stringsAsFactors = FALSE
    )
  )

  ggplot2::ggplot(combined,
                  ggplot2::aes(x = .data$missing_rate)) +
    ggplot2::geom_histogram(bins = 30, fill = "#2C3E99", color = "white") +
    ggplot2::facet_wrap(~ axis, scales = "free", ncol = 1) +
    ggplot2::scale_x_continuous(labels = scales::label_percent()) +
    ggplot2::labs(
      title = "Missingness distribution",
      x = "Missing rate",
      y = "Count"
    ) +
    theme_omicsCore()
}

plot_qc_pca <- function(bundle, color_by = NULL) {
  cleaned <- bundle$results$cleaned_input
  mat <- mean_impute_rows(cleaned$expr_mat)
  if (ncol(mat) < 2L) {
    stop("Need at least 2 samples to draw a PCA scatter.")
  }
  pca <- stats::prcomp(t(mat), scale. = TRUE)
  scores <- as.data.frame(pca$x[, 1:2, drop = FALSE])
  scores$sample_id <- rownames(scores)
  rownames(scores) <- NULL

  if (!is.null(color_by)) {
    if (!color_by %in% colnames(cleaned$meta_df)) {
      stop("`color_by` not found in cleaned_input$meta_df: ", color_by)
    }
    scores[[color_by]] <- cleaned$meta_df[scores$sample_id, color_by]
  }

  var_pct <- (pca$sdev^2) / sum(pca$sdev^2) * 100

  mapping <- if (is.null(color_by)) {
    ggplot2::aes(x = .data$PC1, y = .data$PC2)
  } else {
    ggplot2::aes(x = .data$PC1, y = .data$PC2,
                 color = .data[[color_by]])
  }

  ggplot2::ggplot(scores, mapping) +
    ggplot2::geom_point(size = 2.5, alpha = 0.9) +
    ggplot2::labs(
      title = "PCA of cleaned input",
      x = sprintf("PC1 (%.1f%%)", var_pct[1L]),
      y = sprintf("PC2 (%.1f%%)", var_pct[2L])
    ) +
    theme_omicsCore()
}

plot_qc_connectivity <- function(bundle) {
  out <- bundle$results$qc_summary$outliers
  stats_df <- if (!is.null(out$by_method) && "connectivity" %in% names(out$by_method)) {
    out$by_method$connectivity$stats
  } else if (identical(out$method, "connectivity")) {
    out$stats
  } else {
    # outlier_method was something else; recompute connectivity on the
    # cleaned input so users always see this view.
    cleaned <- bundle$results$cleaned_input
    qc_outliers_connectivity(cleaned$expr_mat,
                             sd_threshold = bundle$params$outlier_sd_threshold)$stats
  }
  stats_df <- stats_df[order(stats_df$mean_correlation), , drop = FALSE]
  stats_df$sample_id <- factor(stats_df$sample_id, levels = stats_df$sample_id)

  ggplot2::ggplot(stats_df,
                  ggplot2::aes(x = .data$sample_id,
                               y = .data$mean_correlation,
                               fill = .data$is_outlier)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = c(`TRUE` = "#C0392B", `FALSE` = "#2C3E99")) +
    ggplot2::labs(
      title = "Sample connectivity",
      x = NULL,
      y = "Mean pairwise correlation",
      fill = "Outlier"
    ) +
    theme_omicsCore() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1))
}

plot_qc_imputation <- function(bundle) {
  cleaned <- bundle$results$cleaned_input
  if (is.null(cleaned$raw_mat)) {
    stop("This bundle has no imputation step (run_qc was called with impute_method='none').")
  }
  before <- as.numeric(cleaned$raw_mat)
  after <- as.numeric(cleaned$expr_mat)
  df <- data.frame(
    value = c(before, after),
    type  = c(rep("raw", length(before)), rep("imputed", length(after))),
    stringsAsFactors = FALSE
  )
  df <- df[is.finite(df$value), , drop = FALSE]

  ggplot2::ggplot(df,
                  ggplot2::aes(x = .data$value, fill = .data$type, color = .data$type)) +
    ggplot2::geom_density(alpha = 0.35) +
    ggplot2::scale_fill_manual(values = c(raw = "#9AA3AE", imputed = "#1FBF9E")) +
    ggplot2::scale_color_manual(values = c(raw = "#9AA3AE", imputed = "#1FBF9E")) +
    ggplot2::labs(
      title = "Imputation effect on intensity distribution",
      x = "Value",
      y = "Density",
      fill = NULL, color = NULL
    ) +
    theme_omicsCore()
}

#' Project-standard ggplot2 theme
#'
#' Matches the visual identity defined in the omicsApp design tokens.
#'
#' @param base_size Base font size.
#' @param base_family Base font family.
#'
#' @return A ggplot2 theme.
#' @export
theme_omicsCore <- function(base_size = 11, base_family = "") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.05)),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "#E5E7EB"),
      axis.line = ggplot2::element_line(color = "#9AA3AE"),
      strip.text = ggplot2::element_text(face = "bold")
    )
}
