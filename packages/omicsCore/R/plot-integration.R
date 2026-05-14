#' Integration summary plot
#'
#' Visualises an `analysis_bundle` produced by [run_integration()]. The
#' available views depend on the method that produced the bundle:
#'
#' * `"scatter"` -- always available. For `correlation` plots the
#'   correlation coefficient against `-log10(adj_p_value)`. For
#'   `concordance` plots `effect_a` vs `effect_b` reconstructed from the
#'   `effect` column (which carries `effect_a - effect_b`). For
#'   `active_pathways` plots `-log10(adj_p_value)` against pathway rank.
#' * `"dual_volcano"` -- concordance-only. Plots `effect` (the difference
#'   of effects) on the x-axis against `-log10(p)` on the y-axis and
#'   colors by quadrant.
#' * `"quadrant"` -- concordance-only. Bar count of the four
#'   `(direction_a, direction_b)` sign quadrants.
#' * `"dotplot"` -- active_pathways-only. Dotplot of top pathways, with
#'   color = adjusted p-value and shape = shared / unique evidence.
#'
#' @param bundle An [`analysis_bundle`][is_analysis_bundle()] produced by
#'   [run_integration()].
#' @param view One of `"scatter"`, `"dual_volcano"`, `"quadrant"`,
#'   `"dotplot"`.
#' @param top_n Number of features / pathways to label or display.
#' @param label_features Optional character vector of `feature_symbol`
#'   values to force-label (scatter / dual_volcano views).
#' @param p_cutoff Significance cutoff used for highlight color in scatter
#'   and dual_volcano views.
#'
#' @return A `ggplot` object.
#' @export
#' @family integration
plot_integration <- function(
  bundle,
  view = c("scatter", "dual_volcano", "quadrant", "dotplot"),
  top_n = 20L,
  label_features = NULL,
  p_cutoff = 0.05
) {
  view <- match.arg(view)
  df <- integration_result_from_bundle(bundle)
  method <- bundle$params$method
  if (is.null(method)) {
    stop("Bundle is missing `params$method`.")
  }

  if (view == "dual_volcano" && method != "concordance") {
    stop("`dual_volcano` view requires method = 'concordance'.")
  }
  if (view == "quadrant" && method != "concordance") {
    stop("`quadrant` view requires method = 'concordance'.")
  }
  if (view == "dotplot" && method != "active_pathways") {
    stop("`dotplot` view requires method = 'active_pathways'.")
  }

  if (nrow(df) == 0L) {
    return(empty_integration_plot("No integration rows to plot."))
  }

  switch(view,
    scatter      = plot_integration_scatter(df, bundle, top_n, label_features, p_cutoff),
    dual_volcano = plot_integration_dual_volcano(df, bundle, top_n, label_features, p_cutoff),
    quadrant     = plot_integration_quadrant(df, bundle),
    dotplot      = plot_integration_dotplot(df, bundle, top_n)
  )
}

# ---- internal helpers --------------------------------------------------

integration_result_from_bundle <- function(bundle) {
  if (!is_analysis_bundle(bundle) ||
      !identical(bundle$analysis_name, "run_integration")) {
    stop("`bundle` must be an analysis_bundle from run_integration().")
  }
  df <- bundle$results$integration_df
  if (is.null(df)) stop("Bundle is missing `results$integration_df`.")
  check_integration_result_schema(df)
  df
}

empty_integration_plot <- function(label) {
  ggplot2::ggplot() +
    ggplot2::theme_void() +
    ggplot2::annotate("text", x = 0.5, y = 0.5, label = label,
                      color = "#4D4D4D", size = 4) +
    ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1)
}

integration_axis_label <- function(bundle, side) {
  exps <- bundle$params$experiments
  if (is.null(exps) || length(exps) < 2L) {
    return(if (side == "a") "experiment A" else "experiment B")
  }
  exps[[if (side == "a") 1L else 2L]]
}

pick_label_ids <- function(df, top_n, label_features, p_col = "adj_p_value") {
  ranked <- df[!is.na(df[[p_col]]), , drop = FALSE]
  ranked <- ranked[order(ranked[[p_col]]), , drop = FALSE]
  top_ids <- utils::head(ranked$feature_id, top_n)
  forced_ids <- if (is.null(label_features)) character(0) else {
    df$feature_id[df$feature_symbol %in% label_features]
  }
  unique(c(top_ids, forced_ids))
}

plot_integration_scatter <- function(df, bundle, top_n, label_features, p_cutoff) {
  method <- bundle$params$method
  df$.neglog10p <- -log10(pmax(df$adj_p_value, .Machine$double.xmin))
  df$.sig <- factor(
    ifelse(!is.na(df$is_significant) & df$is_significant, "significant", "ns"),
    levels = c("ns", "significant")
  )

  label_ids <- pick_label_ids(df, top_n, label_features)
  df$.label <- ifelse(df$feature_id %in% label_ids, df$feature_symbol, NA_character_)

  if (method == "correlation") {
    x_aes <- "effect"
    xlab <- paste0("correlation (", df$effect_type[[1L]], ")")
    title <- "Integration: correlation"
  } else if (method == "concordance") {
    x_aes <- "effect"
    xlab <- paste0("effect difference (", integration_axis_label(bundle, "a"),
                   " - ", integration_axis_label(bundle, "b"), ")")
    title <- "Integration: concordance"
  } else {
    df$.rank <- seq_len(nrow(df))
    x_aes <- ".rank"
    xlab <- "pathway rank"
    title <- "Integration: ActivePathways"
  }

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data[[x_aes]], y = .data$.neglog10p, color = .data$.sig)
  ) +
    ggplot2::geom_point(alpha = 0.75, size = 1.6, na.rm = TRUE) +
    ggplot2::scale_color_manual(
      values = c(ns = "#9AA3AE", significant = "#C0392B"),
      name = NULL
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(p_cutoff), linetype = "dashed", color = "#9AA3AE"
    ) +
    ggplot2::labs(
      title = title,
      subtitle = paste(bundle$params$experiments, collapse = " vs "),
      x = xlab,
      y = "-log10(adj_p_value)"
    ) +
    theme_omicsCore()

  p + add_repel_layer(df, x_aes, ".neglog10p", ".label")
}

plot_integration_dual_volcano <- function(df, bundle, top_n, label_features, p_cutoff) {
  df$.neglog10p <- -log10(pmax(df$p_value, .Machine$double.xmin))
  df$.quad <- ifelse(is.na(df$quadrant), "n/a", df$quadrant)
  label_ids <- pick_label_ids(df, top_n, label_features, p_col = "p_value")
  df$.label <- ifelse(df$feature_id %in% label_ids, df$feature_symbol, NA_character_)

  quadrant_colors <- c(
    up_up     = "#C0392B",
    down_down = "#2C7BB6",
    up_down   = "#8E44AD",
    down_up   = "#16A085",
    `n/a`     = "#9AA3AE"
  )

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data$effect, y = .data$.neglog10p, color = .data$.quad)
  ) +
    ggplot2::geom_point(alpha = 0.8, size = 1.6, na.rm = TRUE) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "#9AA3AE") +
    ggplot2::geom_hline(
      yintercept = -log10(p_cutoff), linetype = "dashed", color = "#9AA3AE"
    ) +
    ggplot2::scale_color_manual(values = quadrant_colors, name = "quadrant",
                                na.value = "#9AA3AE") +
    ggplot2::labs(
      title = "Integration: dual volcano",
      subtitle = paste(bundle$params$experiments, collapse = " vs "),
      x = paste0("effect (", integration_axis_label(bundle, "a"),
                 " - ", integration_axis_label(bundle, "b"), ")"),
      y = "-log10(combined p)"
    ) +
    theme_omicsCore()

  p + add_repel_layer(df, "effect", ".neglog10p", ".label")
}

plot_integration_quadrant <- function(df, bundle) {
  quads <- df$quadrant
  quads[is.na(quads)] <- "n/a"
  levels_order <- c("up_up", "down_down", "up_down", "down_up", "n/a")
  counts <- as.data.frame(table(factor(quads, levels = levels_order)),
                          stringsAsFactors = FALSE)
  names(counts) <- c("quadrant", "n")

  quadrant_colors <- c(
    up_up     = "#C0392B",
    down_down = "#2C7BB6",
    up_down   = "#8E44AD",
    down_up   = "#16A085",
    `n/a`     = "#9AA3AE"
  )

  ggplot2::ggplot(
    counts,
    ggplot2::aes(x = .data$quadrant, y = .data$n, fill = .data$quadrant)
  ) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = quadrant_colors, guide = "none") +
    ggplot2::labs(
      title = "Integration: concordance quadrants",
      subtitle = paste(bundle$params$experiments, collapse = " vs "),
      x = NULL, y = "features"
    ) +
    theme_omicsCore()
}

plot_integration_dotplot <- function(df, bundle, top_n) {
  df <- df[!is.na(df$adj_p_value), , drop = FALSE]
  df <- df[order(df$adj_p_value), , drop = FALSE]
  df <- utils::head(df, top_n)
  if (nrow(df) == 0L) {
    return(empty_integration_plot("No pathways to plot."))
  }
  df$.label <- truncate_pathway_name(df$feature_symbol)
  df$.label <- factor(df$.label, levels = unique(df$.label[order(-df$adj_p_value)]))

  ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data$effect, y = .data$.label,
                 color = .data$adj_p_value, shape = .data$direction)
  ) +
    ggplot2::geom_point(size = 4, na.rm = TRUE) +
    ggplot2::scale_color_gradient(low = "#C0392B", high = "#9AA3AE",
                                  name = "adj p") +
    ggplot2::scale_shape_manual(values = c(shared = 16, unique = 1),
                                na.value = 4, name = "evidence") +
    ggplot2::labs(
      title = "Integration: ActivePathways",
      subtitle = paste(bundle$params$experiments, collapse = " vs "),
      x = "-log10(adj p)", y = NULL
    ) +
    theme_omicsCore()
}
