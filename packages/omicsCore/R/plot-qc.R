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
                    view = c("missing", "depth", "pca", "connectivity",
                             "imputation"),
                    color_by = NULL,
                    ...) {
  if (!is_analysis_bundle(bundle) || !identical(bundle$analysis_name, "run_qc")) {
    stop("`bundle` must be an analysis_bundle from run_qc().")
  }
  view <- match.arg(view)

  switch(view,
    missing      = plot_qc_missing(bundle),
    depth        = plot_qc_depth(bundle),
    pca          = plot_qc_pca(bundle, color_by = color_by),
    connectivity = plot_qc_connectivity(bundle),
    imputation   = plot_qc_imputation(bundle)
  )
}

# ---- views -------------------------------------------------------------

MISSING_FILL <- "#2C3E99"

# How many samples the per-sample panel will draw before it starts
# showing only the worst. Past this the bars are too thin to read in
# the height the app gives the panel, and the question the panel
# answers -- which samples are worst -- only needs the top of the list.
MISSING_MAX_SAMPLE_BARS <- 30L

# Samples and features are not the same kind of question, and were
# previously drawn as though they were: one 30-bin histogram, faceted.
# With a dozen samples most of those bins are empty and the few that are
# not read as unexplained spikes, and a histogram discards the one thing
# the sample panel is for -- *which* sample is bad.
#
# So: a bar per sample, worst first; a density for features, where the
# count is large and the shape is the point.
plot_qc_missing <- function(bundle) {
  miss <- bundle$results$qc_summary$missingness
  # One x range for both panels. They measure the same quantity, and
  # the previous facet_wrap(scales = "free") gave them separate ones,
  # so a sample panel and a feature panel that looked alike could be an
  # order of magnitude apart. Scaled to the data rather than fixed at
  # [0, 1]: real missingness is usually a few percent, and a panel that
  # is 95% empty space hides the shape it exists to show.
  upper <- missing_axis_upper(c(miss$sample_metrics$missing_rate,
                                miss$feature_metrics$missing_rate))
  patchwork::wrap_plots(
    plot_missing_by_sample(miss$sample_metrics, upper),
    plot_missing_by_feature(miss$feature_metrics, upper),
    ncol = 1
  )
}

# Always anchored at 0, so the reader can see where the floor is, and
# never past 1. The floor of 5% stops an all-but-complete dataset from
# being magnified into what looks like a problem.
missing_axis_upper <- function(rates) {
  rates <- rates[!is.na(rates)]
  if (!length(rates)) return(1)
  max(0.05, min(1, max(rates) * 1.15 + 0.01))
}

# How many of the worst samples to name when there are too many to name
# them all. Five fits in the tail of a sorted curve without the labels
# colliding.
MISSING_LABEL_WORST <- 5L

# One bar per sample stops working somewhere past a few dozen: the bars
# get thinner than their labels are tall. Truncating to the worst N
# answers "which sample is bad" but throws away "how bad is this cohort
# overall", which is the other half of what the panel is for.
#
# So past the cap it becomes a sorted curve: every sample is still a
# point, rank on y, and the shape of the curve is the distribution -- a
# flat line with a short tail reads very differently from a steady
# slope. The worst few keep their names, which is all anyone reads off
# the top of a bar chart anyway.
#
# Missing rate stays on x either way, so this panel and the feature
# panel below it remain directly comparable.
plot_missing_by_sample <- function(sample_df, upper = 1) {
  df <- sample_df[order(sample_df$missing_rate, decreasing = TRUE), ,
                  drop = FALSE]
  if (nrow(df) > MISSING_MAX_SAMPLE_BARS) {
    return(plot_missing_sample_curve(df, upper))
  }
  # Reversed, because a discrete y axis is drawn bottom-up and the
  # worst sample belongs at the top.
  df$sample_id <- factor(df$sample_id, levels = rev(df$sample_id))

  ggplot2::ggplot(df, ggplot2::aes(x = .data$missing_rate,
                                   y = .data$sample_id)) +
    ggplot2::geom_col(fill = MISSING_FILL, width = 0.7) +
    ggplot2::scale_x_continuous(
      labels = scales::label_percent(),
      limits = c(0, upper),
      expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::labs(title = "Missing rate per sample",
                  subtitle = sprintf("%d samples", nrow(df)),
                  x = NULL, y = NULL) +
    theme_omicsCore()
}

# `df` is already sorted worst-first.
plot_missing_sample_curve <- function(df, upper = 1) {
  df$rank <- seq_len(nrow(df))
  worst <- df[seq_len(min(MISSING_LABEL_WORST, nrow(df))), , drop = FALSE]

  ggplot2::ggplot(df, ggplot2::aes(x = .data$missing_rate, y = .data$rank)) +
    ggplot2::geom_point(colour = MISSING_FILL, size = 1.1, alpha = 0.7) +
    # Rank 1 is the worst, and belongs at the top.
    ggplot2::scale_y_reverse(breaks = scales::breaks_pretty(4)) +
    ggplot2::scale_x_continuous(labels = scales::label_percent(),
                                limits = c(0, upper)) +
    ggplot2::labs(
      title = "Missing rate per sample",
      # The names go in the subtitle rather than beside their points.
      # The worst few sit at almost the same rank, so on the plot their
      # labels land on top of one another; and putting them here costs
      # no axis room, which keeps this panel's x range identical to the
      # feature panel's.
      subtitle = sprintf("%d samples, ranked. Worst: %s",
                         nrow(df), paste(worst$sample_id, collapse = ", ")),
      x = NULL, y = "Rank"
    ) +
    theme_omicsCore()
}

plot_missing_by_feature <- function(feature_df, upper = 1) {
  rate <- feature_df$missing_rate
  rate <- rate[!is.na(rate)]
  base <- ggplot2::labs(
    title = "Missing rate across features",
    subtitle = sprintf("%d features", length(rate)),
    # Named, because an unlabelled axis reading 2.5, 5.0, 7.5 next to a
    # panel of counts invites being read as one.
    x = "Missing rate", y = "Density"
  )

  # density() needs spread to estimate a bandwidth from. All-complete
  # data -- the good case -- has none, so say so rather than error.
  if (length(rate) < 3L || length(unique(rate)) < 2L) {
    return(
      ggplot2::ggplot(data.frame(x = rate), ggplot2::aes(x = .data$x)) +
        ggplot2::geom_rug(colour = MISSING_FILL) +
        ggplot2::scale_x_continuous(labels = scales::label_percent(),
                                    limits = c(0, upper)) +
        base +
        ggplot2::labs(subtitle = sprintf(
          "%d features, all at %s", length(rate),
          scales::label_percent()(if (length(rate)) rate[1] else 0))) +
        theme_omicsCore()
    )
  }

  # Evaluated on [0, 1] rather than left to spill past either end: a
  # missing rate below 0 or above 1 is not a thing, and a smooth that
  # draws one invites the reader to believe it.
  dens <- stats::density(rate, from = 0, to = 1)
  curve <- data.frame(x = dens$x, y = dens$y)

  ggplot2::ggplot(curve, ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_area(fill = MISSING_FILL, alpha = 0.25) +
    ggplot2::geom_line(colour = MISSING_FILL, linewidth = 0.6) +
    # The rug keeps the individual features visible under the smooth,
    # which matters when there are few of them and the curve is mostly
    # bandwidth.
    ggplot2::geom_rug(data = data.frame(x = rate, y = 0),
                      ggplot2::aes(x = .data$x), sides = "b",
                      alpha = 0.4, colour = MISSING_FILL) +
    # coord_cartesian, not scale limits: the density is estimated over
    # the whole [0, 1] and then viewed, rather than re-estimated from a
    # truncated sample, which would change the curve's shape.
    ggplot2::scale_x_continuous(labels = scales::label_percent()) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::coord_cartesian(xlim = c(0, upper)) +
    base +
    theme_omicsCore()
}

plot_qc_pca <- function(bundle, color_by = NULL) {
  cleaned <- bundle$results$cleaned_input
  mat <- mean_impute_rows(cleaned$expr_mat)
  if (ncol(mat) < 2L) {
    stop("Need at least 2 samples to draw a PCA scatter.")
  }
  pca <- pca_over_samples(mat)
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

  # Said on the plot rather than left implicit. A counts matrix carries
  # every annotated gene, so a fifth of them being zero in every sample
  # is ordinary -- but "PCA of 63,241 features" and "PCA of the 49,000
  # that vary" are different claims, and only one of them is true.
  dropped <- attr(pca, "n_dropped") %||% 0L
  subtitle <- if (dropped > 0L) {
    sprintf("%s features; %s constant across all samples, excluded",
            format(nrow(mat), big.mark = ","),
            format(dropped, big.mark = ","))
  } else NULL

  ggplot2::ggplot(scores, mapping) +
    ggplot2::geom_point(size = 2.5, alpha = 0.9) +
    ggplot2::labs(
      title = "PCA of cleaned input",
      subtitle = subtitle,
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

# ---- depth (RNA-seq) ---------------------------------------------------
# What replaces the missingness panels for a counts matrix, where
# missingness is 0% for every feature and the panel says nothing.
#
# Two questions, drawn together because the pair is what distinguishes
# the two problems: a shallow library drops both total counts and genes
# detected, while a degraded sample drops detection with the total
# holding up.

DEPTH_LOW_RATIO <- 0.3

plot_qc_depth <- function(bundle) {
  depth <- bundle$results$qc_summary$depth
  if (is.null(depth) || nrow(depth) == 0L) {
    stop("This QC bundle carries no depth summary.", call. = FALSE)
  }
  patchwork::wrap_plots(
    plot_depth_library(depth),
    plot_depth_detection(depth),
    ncol = 1
  )
}

# Ordered worst-first for the same reason the missingness panel is: the
# question is which sample is bad, and sorting is what answers it
# without reading every label.
depth_ordered <- function(depth, col) {
  depth <- depth[order(depth[[col]]), , drop = FALSE]
  depth$sample_id <- factor(depth$sample_id, levels = depth$sample_id)
  depth
}

plot_depth_library <- function(depth) {
  d <- depth_ordered(depth, "library_size")
  # Coloured against the median rather than an absolute count: what
  # counts as shallow depends on the experiment, and a fixed cutoff
  # would be wrong for every study but one.
  d$.low <- !is.na(d$library_size_ratio) &
    d$library_size_ratio < DEPTH_LOW_RATIO

  ggplot2::ggplot(d, ggplot2::aes(x = .data$library_size,
                                  y = .data$sample_id,
                                  fill = .data$.low)) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::scale_fill_manual(
      values = c(`FALSE` = MISSING_FILL, `TRUE` = omics_colors$up),
      guide = "none") +
    ggplot2::scale_x_continuous(labels = depth_axis_labels) +
    ggplot2::labs(
      title = "Library size per sample",
      subtitle = depth_library_subtitle(d),
      x = "Total counts", y = NULL
    ) +
    theme_omicsCore() +
    depth_sample_axis(nrow(d))
}

plot_depth_detection <- function(depth) {
  d <- depth_ordered(depth, "n_detected")
  ggplot2::ggplot(d, ggplot2::aes(x = .data$n_detected,
                                  y = .data$sample_id)) +
    ggplot2::geom_col(width = 0.75, fill = MISSING_FILL) +
    ggplot2::scale_x_continuous(labels = depth_axis_labels) +
    ggplot2::labs(
      title = "Features detected per sample",
      subtitle = sprintf("of %s in the matrix",
                         format(round(max(d$n_detected) /
                                        max(d$detection_rate, na.rm = TRUE)),
                                big.mark = ",")),
      x = "Features with any signal", y = NULL
    ) +
    theme_omicsCore() +
    depth_sample_axis(nrow(d))
}

depth_library_subtitle <- function(d) {
  n_low <- sum(d$.low)
  if (n_low == 0L) {
    "no sample below 30% of the median"
  } else {
    sprintf("%d sample%s below 30%% of the median: %s",
            n_low, if (n_low == 1L) "" else "s",
            paste(as.character(d$sample_id[d$.low]), collapse = ", "))
  }
}

# Past a few dozen samples the labels stop being legible in the height
# the app gives the panel, and the ordering is what carries the meaning
# anyway.
depth_sample_axis <- function(n) {
  if (n <= MISSING_MAX_SAMPLE_BARS) return(NULL)
  ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                 axis.ticks.y = ggplot2::element_blank())
}

depth_axis_labels <- function(x) {
  ifelse(is.na(x), "",
         ifelse(abs(x) >= 1e6, paste0(round(x / 1e6, 1), "M"),
                ifelse(abs(x) >= 1e3, paste0(round(x / 1e3), "k"),
                       format(x, big.mark = ","))))
}
