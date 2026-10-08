#' Volcano plot for a diff bundle
#'
#' x-axis is `effect` (log2FC, beta, correlation, depending on the backend);
#' y-axis is `-log10(p)` of the adjusted p-value by default. Features that
#' pass `p_threshold` (and `effect_threshold`, when given) are coloured by
#' direction -- up (red) or down (blue); `"positive"` / `"negative"` for a
#' continuous variable -- and drawn over the features that do not, which
#' are smaller and grey. The legend counts each class. A result whose
#' effect has no sign (the F statistic of a global test) has a single
#' "significant" class instead. For a signed effect the x axis is
#' symmetric about zero, so up and down magnitudes compare at a glance.
#' The top `top_n` rows by (adjusted) p-value are labelled. Optionally
#' supply `label_features` to force-label a specific set of feature
#' symbols.
#'
#' @param bundle An `analysis_bundle` produced by [run_diff()].
#' @param top_n Number of top features to label, ranked by `p_basis`.
#' @param label_features Optional character vector of `feature_symbol` values
#'   to always label.
#' @param p_basis Which p-value column to use for the y-axis and the
#'   significance cut, `"adjusted"` or `"raw"`.
#' @param effect_threshold Optional absolute-effect cutoff: a feature must
#'   also reach it to count as significant, and dashed vertical lines mark
#'   it on both sides of zero.
#' @param p_threshold P-value cutoff, drawn as a dashed horizontal line, or
#'   `NULL` for none.
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
  assert_count(top_n, "top_n")
  assert_character(label_features, "label_features", allow_null = TRUE)
  assert_number(effect_threshold, "effect_threshold", lower = 0, allow_null = TRUE)
  assert_number(p_threshold, "p_threshold", lower = 0, upper = 1, allow_null = TRUE)
  result_df <- diff_result_from_bundle(bundle)
  p_basis <- match.arg(p_basis)
  p_col <- resolve_p_col(result_df, p_preference = p_basis)

  df <- result_df
  df$.neglog10p <- -log10(pmax(df[[p_col]], .Machine$double.xmin))

  sig <- diff_significance(df, p_col, p_threshold, effect_threshold)
  df$.sig <- factor(ifelse(sig, "significant", "ns"),
                    levels = c("ns", "significant"))
  kind <- diff_direction_kind(df)
  classes <- diff_direction_classes(df$effect, sig, kind)
  df$.class <- classes$class

  # Choose label set: union of forced labels and top_n by p-basis.
  ranked <- df[order(df[[p_col]], na.last = NA), , drop = FALSE]
  top_ids <- utils::head(ranked$feature_id, top_n)
  forced_ids <- if (is.null(label_features)) character(0) else {
    df$feature_id[df$feature_symbol %in% label_features]
  }
  df$.label <- ifelse(df$feature_id %in% unique(c(top_ids, forced_ids)),
                      df$feature_symbol, NA_character_)

  # What an interactive viewer (plotly) shows on hover: the gene first.
  # The static figure ignores it.
  # In the app's words (log2FC, adjusted p), not the column names.
  df$.hover <- sprintf("%s<br>%s: %.3f<br>%s: %.3g",
                       ifelse(is.na(df$feature_symbol) | !nzchar(df$feature_symbol),
                              df$feature_id, df$feature_symbol),
                       effect_label(bundle), df$effect,
                       switch(p_col, adj_p_value = "adjusted p", p_value = "p", p_col),
                       df[[p_col]])
  p <- ggplot2::ggplot(df,
                       ggplot2::aes(x = .data$effect,
                                    y = .data$.neglog10p,
                                    color = .data$.class,
                                    text = .data$.hover)) +
    diff_point_layers(df, classes) +
    ggplot2::labs(
      title = "Volcano",
      subtitle = volcano_subtitle(bundle),
      # What "significant" meant here travels with the figure. Without
      # it a reader has a two-coloured cloud and no way to know which
      # cut produced it -- and a screenshot outlives the session that
      # set the controls.
      caption = threshold_caption(p_col, p_threshold, effect_threshold,
                                effect_label(bundle)),
      x = volcano_xlab(bundle),
      y = p_axis_label(p_col)
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
    # An F statistic is never negative: a line at -cut would mark
    # nothing and push the axis below zero.
    p <- p + ggplot2::geom_vline(
      xintercept = if (identical(kind, "unsigned")) effect_threshold
                   else c(-effect_threshold, effect_threshold),
      linetype = "dashed", color = omics_colors$ns
    )
  }
  # Symmetric about zero, so a reader compares how far up and how far
  # down things moved without reading two different tick ranges. The
  # threshold is included so its lines are never cut off when nothing
  # reaches it. A coordinate limit drops no point (a scale limit would),
  # and unlike expand_limits() it adds no invisible layer for ggplotly()
  # to turn into a trace that toWebGL() then cannot convert.
  # An F statistic has no sign and nothing to mirror.
  if (!identical(kind, "unsigned")) {
    reach <- c(abs(df$effect), effect_threshold)
    reach <- reach[is.finite(reach)]
    if (length(reach)) {
      p <- p + ggplot2::coord_cartesian(xlim = c(-max(reach), max(reach)))
    }
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
#' Significant features are coloured by direction, as in [plot_volcano()].
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
  assert_count(top_n, "top_n")
  assert_character(label_features, "label_features", allow_null = TRUE)
  assert_number(effect_threshold, "effect_threshold", lower = 0, allow_null = TRUE)
  assert_number(p_threshold, "p_threshold", lower = 0, upper = 1, allow_null = TRUE)
  result_df <- diff_result_from_bundle(bundle)
  if (all(is.na(result_df$base_mean))) {
    stop("`base_mean` is all NA in this bundle; MA plot is not available.")
  }
  p_basis <- match.arg(p_basis)
  p_col <- resolve_p_col(result_df, p_preference = p_basis)

  df <- result_df
  sig <- diff_significance(df, p_col, p_threshold, effect_threshold)
  df$.sig <- factor(ifelse(sig, "significant", "ns"),
                    levels = c("ns", "significant"))
  # Coloured by direction like the volcano beside it, so the same
  # feature is the same colour in both.
  classes <- diff_direction_classes(df$effect, sig, diff_direction_kind(df))
  df$.class <- classes$class

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
                                    color = .data$.class)) +
    diff_point_layers(df, classes) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = omics_colors$ns) +
    ggplot2::labs(
      title = "MA plot",
      subtitle = volcano_subtitle(bundle),
      caption = threshold_caption(p_col, p_threshold, effect_threshold,
                                effect_label(bundle)),
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
  assert_string(color_by, "color_by", allow_null = TRUE)
  assert_string(shape_by, "shape_by", allow_null = TRUE)
  assert_flag(log2, "log2", allow_null = TRUE)
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
  group_vals <- if (is.null(color_by)) NULL else scores[[color_by]]
  redundant <- is.null(shape_by) && use_group_shape(group_vals)
  if (redundant) mapping$shape <- ggplot2::aes(shape = .data[[color_by]])$shape
  legend <- group_legend_scales(group_vals, redundant_shape = redundant)
  if (!is.null(shape_by)) {
    legend <- c(legend, list(ggplot2::scale_shape_discrete(
      labels = function(x) wrap_label(x, width = 18L, max_lines = 3L))))
  }

  ggplot2::ggplot(scores, mapping) +
    ggplot2::geom_point(size = 2.5, alpha = 0.9) +
    ggplot2::labs(
      title = "PCA scores",
      x = sprintf("PC1 (%.1f%%)", var_pct[1L]),
      y = sprintf("PC2 (%.1f%%)", var_pct[2L])
    ) +
    legend +
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
  assert_character(features, "features")
  assert_string(group_by, "group_by")
  assert_string(color_by, "color_by", allow_null = TRUE)
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
    # A pinned jitter: the figure is the same on every render, and the
    # draw leaves the caller's random stream where it was.
    ggplot2::geom_point(
      position = ggplot2::position_jitter(width = 0.18, height = 0, seed = 1L),
      alpha = 0.7, size = 1.6) +
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
  # One figure, one comparison. A bundle holding several drew every
  # feature once per comparison -- 400 points for 200 genes, labels
  # repeated -- under a subtitle naming all of them.
  n_cmp <- length(unique(stats::na.omit(result_df$comparison)))
  if (n_cmp > 1L) {
    stop("The bundle holds ", n_cmp, " comparisons; plot one at a time with ",
         "select_comparison(bundle, \"...\").", call. = FALSE)
  }
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
    # >=, as filter_diff_results() has it: the figure and the hit table
    # used to disagree about a feature exactly at the cutoff.
    sig <- sig & abs(df$effect) >= effect_threshold
  }
  sig[is.na(sig)] <- FALSE
  sig
}

# What the sign of `effect` means for this result. A group comparison
# goes up or down; a continuous variable correlates positively or
# negatively (the words the result's own `direction` column uses); a
# global test's F statistic is never negative and has no direction, so
# calling its hits "up" would assert something the test did not say.
diff_direction_kind <- function(df) {
  effect_type <- first_or_na(df$effect_type)
  analysis_type <- first_or_na(df$analysis_type)
  if (identical(analysis_type, "anova") ||
      (!is.na(effect_type) && grepl("_statistic$", effect_type))) {
    return("unsigned")
  }
  if (!is.na(analysis_type) && startsWith(analysis_type, "continuous")) {
    return("continuous")
  }
  "group"
}

# Colour classes for the points of a volcano or MA plot, as a factor
# whose levels are the legend labels ("up (412)", "down (388)", "not
# significant") and the colour for each. The counts are in the label
# itself, not added by the scale, so ggplotly() -- which names its
# traces from the data, not from scale labels -- shows the same words.
#
# A significant feature with no direction (effect exactly 0, or NA) is
# grouped with "not significant": it cannot be placed on either side,
# and the app's up/down counts leave it out the same way.
diff_direction_classes <- function(effect, sig, kind) {
  words <- switch(kind,
    unsigned   = c(hit = "significant"),
    continuous = c(up = "positive", down = "negative"),
    c(up = "up", down = "down"))
  if (identical(kind, "unsigned")) {
    key <- ifelse(sig, "hit", "ns")
    colours <- c(hit = omics_colors$up, ns = omics_colors$ns)
  } else {
    key <- ifelse(sig & !is.na(effect) & effect > 0, "up",
                  ifelse(sig & !is.na(effect) & effect < 0, "down", "ns"))
    colours <- c(up = omics_colors$up, down = omics_colors$down,
                 ns = omics_colors$ns)
  }
  hits <- names(words)
  labels <- c(sprintf("%s (%s)", words,
                      format(vapply(hits, function(h) sum(key == h), integer(1)),
                             big.mark = ",", trim = TRUE)),
              "not significant")
  names(labels) <- c(hits, "ns")
  # Legend order: the hits first, as the reader looks for them.
  class <- factor(unname(labels[key]), levels = unname(labels))
  list(class = class,
       colours = stats::setNames(unname(colours[names(labels)]), labels),
       ns_label = labels[["ns"]])
}

# The points of a volcano or MA plot: the features that pass nothing
# first, smaller and fainter, then the hits on top -- one layer drew
# the grey cloud over red points wherever they overlapped, and at full
# strength the cloud competed with the hits for attention.
diff_point_layers <- function(df, classes) {
  is_ns <- df$.class == classes$ns_label
  # A layer with no rows is left out: ggplotly() turns it into an
  # invisible trace that toWebGL() cannot convert, and plotly warns.
  # The legend does not need it -- the scale's limits keep every key.
  ns_layer <- if (any(is_ns) || !any(!is_ns)) {
    ggplot2::geom_point(data = df[is_ns, , drop = FALSE],
                        alpha = 0.45, size = 1.1, na.rm = TRUE,
                        show.legend = TRUE)
  }
  hit_layer <- if (any(!is_ns)) {
    ggplot2::geom_point(data = df[!is_ns, , drop = FALSE],
                        alpha = 0.85, size = 1.7, na.rm = TRUE,
                        show.legend = TRUE)
  }
  list(
    ns_layer,
    hit_layer,
    # Every class keeps its legend entry and swatch, even with nothing
    # in it (show.legend = TRUE above draws the key of an empty class):
    # "down (0)" says nothing went down, which a missing entry does not.
    ggplot2::scale_color_manual(
      values = classes$colours, breaks = names(classes$colours),
      limits = names(classes$colours), drop = FALSE, name = NULL
    ),
    ggplot2::guides(color = ggplot2::guide_legend(
      override.aes = list(size = 2.4, alpha = 1)))
  )
}

# Axis titles say "adjusted p", the words the controls use, not the
# column name (adj_p_value) the value was read from.
p_axis_label <- function(p_col) {
  paste0("-log10(", switch(p_col, adj_p_value = "adjusted p", p_value = "p", p_col), ")")
}

threshold_caption <- function(p_col, p_threshold, effect_threshold,
                              effect_name = "effect") {
  bits <- character(0)
  if (!is.null(p_threshold)) {
    p_name <- switch(p_col, adj_p_value = "adjusted p", p_value = "p", p_col)
    bits <- c(bits, sprintf("%s < %g", p_name, p_threshold))
  }
  if (!is.null(effect_threshold)) {
    # ">=", as the cut is applied (diff_significance()) and as the app
    # states it; ">" described a rule one feature short of the real one.
    bits <- c(bits, sprintf("|%s| >= %g", effect_name, effect_threshold))
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
  # What base_mean is depends on the engine: DESeq2 reports the mean of
  # normalized counts (linear), edgeR log-CPM, limma the average log value.
  switch(first_or_na(bundle$results$diff_result_df$method) %||% "",
         deseq2 = "mean of normalized counts (base_mean)",
         edger = "average log CPM (base_mean)",
         limma = "average log expression (base_mean)",
         "mean expression (base_mean)")
}

add_repel_layer <- function(df, x, y, label_col, ...) {
  if (!any(!is.na(df[[label_col]]))) return(NULL)
  if (is_installed("ggrepel")) {
    ggrepel::geom_text_repel(
      data = df,
      mapping = ggplot2::aes(x = .data[[x]], y = .data[[y]],
                             label = .data[[label_col]]),
      size = 3, color = omics_colors$fg_dark, max.overlaps = Inf,
      na.rm = TRUE, inherit.aes = FALSE, ...
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

#' Short name of a differential result's effect
#'
#' What the `effect` column of a result is, in the words the app, the
#' plots and the report all use: `"log2FC"` for a group comparison,
#' `"slope"` for a linear trend, `"rho"` for a rank correlation. One
#' name, so a threshold set as "|log2FC| 0.26" is not described as
#' "effect > 0.26" on the next card and "|effect| >= 0.263" on the
#' next page.
#'
#' @param x An `analysis_bundle` from [run_diff()], or an `effect_type`
#'   string from its result table.
#' @return A length-one character string.
#' @export
#' @family diff
#' @examples
#' effect_label("log2FC")
effect_label <- function(x) {
  type <- if (is_analysis_bundle(x)) {
    first_or_na(x$results$diff_result_df$effect_type)
  } else {
    first_or_na(as.character(x))
  }
  if (is.na(type)) return("log2FC")
  switch(type,
    log2FC = "log2FC",
    log2FC_per_unit = "slope",
    beta = "slope",
    mean_diff = "mean difference",
    correlation = "rho",
    F_statistic = "F",
    "effect")
}
