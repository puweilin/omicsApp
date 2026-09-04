#' Volcano plot for a diff bundle
#'
#' x-axis is `effect` (log2FC, beta, correlation, depending on the backend);
#' y-axis is `-log10(p_value)` by default. Significant features
#' (`is_significant`) are highlighted; the top `top_n` rows by (adjusted)
#' p-value are labelled. Optionally supply `label_features` to force-label a
#' specific set of feature symbols.
#'
#' @param bundle An `analysis_bundle` produced by [run_diff()].
#' @param top_n Number of top features to label, ranked by `p_basis`.
#' @param label_features Optional character vector of `feature_symbol` values
#'   to always label.
#' @param p_basis Which p-value column to use for the y-axis,
#'   `"adjusted"` or `"raw"`.
#' @param effect_threshold Vertical reference line for absolute effect.
#' @param p_threshold Horizontal reference p-value line (significance cutoff).
#'
#' @return A `ggplot` object.
#' @export
#' @family diff
plot_volcano <- function(
  bundle,
  top_n = 20,
  label_features = NULL,
  p_basis = c("adjusted", "raw"),
  effect_threshold = NULL,
  p_threshold = 0.05
) {
  result_df <- diff_result_from_bundle(bundle)
  p_basis <- match.arg(p_basis)
  p_col <- resolve_p_col(result_df, p_preference = p_basis)

  df <- result_df
  df$.neglog10p <- -log10(pmax(df[[p_col]], .Machine$double.xmin))

  df$.sig <- factor(
    ifelse(diff_significance(df, p_col, p_threshold, effect_threshold),
           "significant", "ns"),
    levels = c("ns", "significant"))

  # Choose label set: union of forced labels and top_n by p-basis.
  ranked <- df[order(df[[p_col]], na.last = NA), , drop = FALSE]
  top_ids <- utils::head(ranked$feature_id, top_n)
  forced_ids <- if (is.null(label_features)) character(0) else {
    df$feature_id[df$feature_symbol %in% label_features]
  }
  df$.label <- ifelse(df$feature_id %in% unique(c(top_ids, forced_ids)),
                      df$feature_symbol, NA_character_)

  p <- ggplot2::ggplot(df,
                       ggplot2::aes(x = .data$effect,
                                    y = .data$.neglog10p,
                                    color = .data$.sig)) +
    ggplot2::geom_point(alpha = 0.75, size = 1.6, na.rm = TRUE) +
    ggplot2::scale_color_manual(
      values = c(ns = omics_colors$ns, significant = omics_colors$up),
      name = NULL
    ) +
    ggplot2::labs(
      title = "Volcano",
      subtitle = volcano_subtitle(bundle),
      # What "significant" meant here travels with the figure. Without
      # it a reader has a two-coloured cloud and no way to know which
      # cut produced it -- and a screenshot outlives the session that
      # set the controls.
      caption = threshold_caption(p_col, p_threshold, effect_threshold),
      x = volcano_xlab(bundle),
      y = paste0("-log10(", p_col, ")")
    ) +
    theme_omicsCore()

  # Both rules are conditional. `p_threshold = NULL` means "do not draw
  # one", and drawing it unconditionally made an explicit NULL crash in
  # log10() rather than doing the obvious thing.
  if (!is.null(p_threshold)) {
    p <- p + ggplot2::geom_hline(
      yintercept = -log10(p_threshold),
      linetype = "dashed", color = omics_colors$ns
    )
  }
  if (!is.null(effect_threshold)) {
    p <- p + ggplot2::geom_vline(
      xintercept = c(-effect_threshold, effect_threshold),
      linetype = "dashed", color = omics_colors$ns
    )
  }

  p + add_repel_layer(df, "effect", ".neglog10p", ".label")
}

#' MA plot for a diff bundle
#'
#' Plots `effect` against `base_mean` so the user can spot effect-size
#' biases concentrated at low- or high-expressed features.
#'
#' Which features count as significant is decided here, from the
#' thresholds this call is given — the same contract as [plot_volcano()],
#' and for the same reason: `run_diff()` applies no cutoff, so its
#' `is_significant` column is `NA` and has nothing to colour by.
#'
#' @param bundle An `analysis_bundle` produced by [run_diff()].
#' @param top_n Number of top features to label.
#' @param label_features Optional character vector of `feature_symbol` values
#'   to always label.
#' @param p_basis Whether to threshold on the adjusted or raw p-value.
#' @param effect_threshold Optional absolute-effect cutoff.
#' @param p_threshold P-value cutoff, or `NULL` to make no distinction.
#'
#' @return A `ggplot` object.
#' @export
#' @family diff
plot_ma <- function(bundle, top_n = 20, label_features = NULL,
                    p_basis = c("adjusted", "raw"),
                    effect_threshold = NULL,
                    p_threshold = 0.05) {
  result_df <- diff_result_from_bundle(bundle)
  if (all(is.na(result_df$base_mean))) {
    stop("`base_mean` is all NA in this bundle; MA plot is not available.")
  }
  p_basis <- match.arg(p_basis)
  p_col <- resolve_p_col(result_df, p_preference = p_basis)

  df <- result_df
  df$.sig <- factor(
    ifelse(diff_significance(df, p_col, p_threshold, effect_threshold),
           "significant", "ns"),
    levels = c("ns", "significant")
  )

  ranked <- df[order(df$p_value, na.last = NA), , drop = FALSE]
  top_ids <- utils::head(ranked$feature_id, top_n)
  forced_ids <- if (is.null(label_features)) character(0) else {
    df$feature_id[df$feature_symbol %in% label_features]
  }
  df$.label <- ifelse(df$feature_id %in% unique(c(top_ids, forced_ids)),
                      df$feature_symbol, NA_character_)

  p <- ggplot2::ggplot(df,
                       ggplot2::aes(x = .data$base_mean,
                                    y = .data$effect,
                                    color = .data$.sig)) +
    ggplot2::geom_point(alpha = 0.75, size = 1.6, na.rm = TRUE) +
    ggplot2::scale_color_manual(
      values = c(ns = omics_colors$ns, significant = omics_colors$up),
      name = NULL
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = omics_colors$ns) +
    ggplot2::labs(
      title = "MA plot",
      subtitle = volcano_subtitle(bundle),
      caption = threshold_caption(p_col, p_threshold, effect_threshold),
      x = ma_xlab(bundle),
      y = volcano_xlab(bundle)
    ) +
    theme_omicsCore()

  p + add_repel_layer(df, "base_mean", "effect", ".label")
}

#' PCA scores plot for an omics_input
#'
#' Mean-imputes NA cells, optionally log2-transforms raw counts, then runs
#' `prcomp(scale. = TRUE)` on samples-by-features and returns a scatter on
#' the first two principal components. Sample metadata supplies optional
#' `color_by` and `shape_by` aesthetics.
#'
#' @param input A validated `omics_input`.
#' @param color_by Optional column in `meta_df` used to color samples.
#' @param shape_by Optional column in `meta_df` used as point shape.
#' @param log2 If `TRUE`, apply `log2(x + 1)` before PCA. Defaults to `TRUE`
#'   when `assay_type == "raw_count"`.
#'
#' @return A `ggplot` object.
#' @export
#' @family diff
plot_pca <- function(input, color_by = NULL, shape_by = NULL, log2 = NULL) {
  validate_omics_input(input)
  if (is.null(log2)) {
    log2 <- identical(input$assay_type, "raw_count")
  }
  mat <- input$expr_mat
  if (isTRUE(log2)) mat <- log2(mat + 1)
  mat <- mean_impute_rows(mat)

  if (ncol(mat) < 2L) {
    stop("Need at least 2 samples to draw a PCA scatter.")
  }
  pca <- pca_over_samples(mat)
  scores <- as.data.frame(pca$x[, 1:2, drop = FALSE])
  scores$sample_id <- rownames(scores)
  rownames(scores) <- NULL
  meta <- input$meta_df

  if (!is.null(color_by)) {
    if (!color_by %in% colnames(meta)) {
      stop("`color_by` not found in `meta_df`: ", color_by)
    }
    scores[[color_by]] <- meta[scores$sample_id, color_by]
  }
  if (!is.null(shape_by)) {
    if (!shape_by %in% colnames(meta)) {
      stop("`shape_by` not found in `meta_df`: ", shape_by)
    }
    scores[[shape_by]] <- factor(meta[scores$sample_id, shape_by])
  }

  var_pct <- (pca$sdev^2) / sum(pca$sdev^2) * 100
  mapping <- ggplot2::aes(x = .data$PC1, y = .data$PC2)
  if (!is.null(color_by)) mapping$colour <- ggplot2::aes(color = .data[[color_by]])$colour
  if (!is.null(shape_by)) mapping$shape  <- ggplot2::aes(shape = .data[[shape_by]])$shape

  ggplot2::ggplot(scores, mapping) +
    ggplot2::geom_point(size = 2.5, alpha = 0.9) +
    ggplot2::labs(
      title = "PCA scores",
      x = sprintf("PC1 (%.1f%%)", var_pct[1L]),
      y = sprintf("PC2 (%.1f%%)", var_pct[2L])
    ) +
    theme_omicsCore()
}

#' Per-feature expression box / violin
#'
#' Plots one or more features as boxplots (with jittered points) split by a
#' grouping column. RNA-seq raw counts are log2-transformed automatically.
#'
#' @param input A validated `omics_input`.
#' @param features Character vector of feature IDs or symbols to plot.
#' @param group_by Column in `meta_df` used for the x-axis grouping.
#' @param color_by Optional column for additional color aesthetic.
#'
#' @return A `ggplot` object.
#' @export
#' @family diff
plot_feature_expression <- function(input, features, group_by, color_by = NULL) {
  validate_omics_input(input)
  if (length(features) == 0L) {
    stop("`features` must contain at least one feature.")
  }
  if (!group_by %in% colnames(input$meta_df)) {
    stop("`group_by` not found in `meta_df`: ", group_by)
  }

  feat_df <- input$feature_df
  if (!"feature_symbol" %in% colnames(feat_df)) {
    feat_df$feature_symbol <- feat_df$feature_id
  }

  # Resolve user-supplied features against both feature_id and feature_symbol.
  feature_ids <- character(0)
  for (f in features) {
    hit <- feat_df$feature_id[feat_df$feature_id == f | feat_df$feature_symbol == f]
    feature_ids <- c(feature_ids, hit)
  }
  feature_ids <- unique(feature_ids)
  if (length(feature_ids) == 0L) {
    stop("None of the requested features were found in `feature_df`.")
  }

  mat <- input$expr_mat
  if (identical(input$assay_type, "raw_count")) {
    mat <- log2(mat + 1)
  }
  mat <- mat[feature_ids, , drop = FALSE]

  long <- data.frame(
    feature_id = rep(rownames(mat), times = ncol(mat)),
    sample_id = rep(colnames(mat), each = nrow(mat)),
    value = as.numeric(mat),
    stringsAsFactors = FALSE
  )
  long$feature_symbol <- feat_df$feature_symbol[match(long$feature_id, feat_df$feature_id)]
  long[[group_by]] <- input$meta_df[long$sample_id, group_by]
  if (!is.null(color_by)) {
    if (!color_by %in% colnames(input$meta_df)) {
      stop("`color_by` not found in `meta_df`: ", color_by)
    }
    long[[color_by]] <- input$meta_df[long$sample_id, color_by]
  }

  base_aes <- ggplot2::aes(x = .data[[group_by]], y = .data$value)
  if (!is.null(color_by)) {
    base_aes$colour <- ggplot2::aes(color = .data[[color_by]])$colour
  }

  ggplot2::ggplot(long, base_aes) +
    ggplot2::geom_boxplot(outlier.shape = NA, fill = "#EAEEF4", color = "#3F4A5A") +
    ggplot2::geom_jitter(width = 0.18, height = 0, alpha = 0.7, size = 1.6) +
    ggplot2::facet_wrap(~ .data$feature_symbol, scales = "free_y") +
    ggplot2::labs(
      title = "Feature expression",
      x = group_by,
      y = if (identical(input$assay_type, "raw_count")) "log2(count + 1)" else "value"
    ) +
    theme_omicsCore()
}

# ---- internal helpers --------------------------------------------------

diff_result_from_bundle <- function(bundle) {
  if (!is_analysis_bundle(bundle) || !identical(bundle$analysis_name, "run_diff")) {
    stop("`bundle` must be an analysis_bundle from run_diff().")
  }
  result_df <- bundle$results$diff_result_df
  if (is.null(result_df)) {
    stop("Bundle is missing `results$diff_result_df`.")
  }
  check_diff_result_schema(result_df)
  result_df
}

volcano_subtitle <- function(bundle) {
  comparison <- bundle$params$comparison
  method <- bundle$params$method
  if (is.null(comparison) || is.null(method)) return(NULL)
  paste0("method = ", method, "  |  comparison = ", comparison)
}

# Which features count as significant, decided from the thresholds the
# caller gave rather than read off `is_significant`.
#
# `run_diff()` applies no cutoff -- it is not told one -- so it writes
# that column NA: "no threshold was applied", which is the truth. It
# used to write FALSE, which asserts "this feature is not significant"
# about a feature with adj.P = 1e-30, and both the volcano and the MA
# plot believed it and drew a single colour over data where a fifth of
# the features cleared 0.05.
#
# The threshold belongs to whoever is asking. `run_integration()` takes
# a `p_cutoff` and so its results carry a real answer; nothing hands one
# to `run_diff()`, so the question is open until a figure or a filter
# closes it -- and then that figure says which cut it used.
diff_significance <- function(df, p_col, p_threshold, effect_threshold) {
  if (is.null(p_threshold) && is.null(effect_threshold)) {
    # Asked to make no distinction: the stored column is all there is,
    # and NA there means the question was never answered.
    return(!is.na(df$is_significant) & df$is_significant)
  }
  sig <- rep(TRUE, nrow(df))
  if (!is.null(p_threshold)) {
    sig <- sig & df[[p_col]] < p_threshold
  }
  if (!is.null(effect_threshold)) {
    sig <- sig & abs(df$effect) > effect_threshold
  }
  sig[is.na(sig)] <- FALSE
  sig
}

threshold_caption <- function(p_col, p_threshold, effect_threshold) {
  bits <- character(0)
  if (!is.null(p_threshold)) {
    bits <- c(bits, sprintf("%s < %g", p_col, p_threshold))
  }
  if (!is.null(effect_threshold)) {
    bits <- c(bits, sprintf("|effect| > %g", effect_threshold))
  }
  if (length(bits) == 0L) {
    return("significance as recorded on the result")
  }
  paste("significant:", paste(bits, collapse = ", "))
}

# An axis label is read off the first row, which is not there when every
# feature was filtered out upstream. That is an empty plot to draw, not
# an error to raise -- the caller asked for a figure of nothing, and a
# subscript error tells them nothing about why.
first_or_na <- function(x) {
  if (is.null(x) || length(x) == 0L) return(NA_character_)
  x[[1L]]
}

volcano_xlab <- function(bundle) {
  effect_type <- first_or_na(bundle$results$diff_result_df$effect_type)
  if (is.na(effect_type)) return("effect")
  switch(effect_type,
    log2FC = "log2 fold change",
    log2FC_per_unit = "log2 fold change / unit",
    beta = "beta",
    mean_diff = "mean difference",
    correlation = "Spearman rho",
    F_statistic = "F statistic",
    effect_type
  )
}

ma_xlab <- function(bundle) {
  omics_type <- first_or_na(bundle$results$diff_result_df$omics_type)
  if (identical(omics_type, "rnaseq")) "log CPM (base_mean)" else "mean expression"
}

add_repel_layer <- function(df, x, y, label_col) {
  if (!any(!is.na(df[[label_col]]))) return(NULL)
  if (is_installed("ggrepel")) {
    ggrepel::geom_text_repel(
      data = df,
      mapping = ggplot2::aes(x = .data[[x]], y = .data[[y]],
                             label = .data[[label_col]]),
      size = 3, color = omics_colors$fg_dark, max.overlaps = Inf,
      na.rm = TRUE, inherit.aes = FALSE
    )
  } else {
    ggplot2::geom_text(
      data = df,
      mapping = ggplot2::aes(x = .data[[x]], y = .data[[y]],
                             label = .data[[label_col]]),
      size = 3, color = omics_colors$fg_dark,
      vjust = -0.6, na.rm = TRUE, inherit.aes = FALSE
    )
  }
}
