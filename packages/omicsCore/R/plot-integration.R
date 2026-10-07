#' Integration summary plot
#'
#' Visualises an `analysis_bundle` produced by [run_integration()]. The
#' available views depend on the method that produced the bundle:
#'
#' * `"scatter"` -- always available. For `correlation` plots the
#'   correlation coefficient against `-log10(adj_p_value)`. For
#'   `concordance` plots the effect difference (`effect_a - effect_b`)
#'   against `-log10(adj_p_value)`. For `active_pathways` plots
#'   `-log10(adj_p_value)` against pathway rank.
#' * `"dual_volcano"` -- concordance-only. Plots `effect` (the difference
#'   of effects) on the x-axis against `-log10(p)` on the y-axis and
#'   colors by quadrant.
#' * `"effect_pair"` -- concordance-only. `effect_a` against `effect_b`,
#'   features that are hits in both layers coloured by quadrant.
#' * `"quadrant"` -- concordance-only. Bar count of the four
#'   `(direction_a, direction_b)` sign quadrants.
#' * `"dotplot"` -- active_pathways-only. Dotplot of top pathways, with
#'   colour = direction (up in both layers, down in both, or layers
#'   disagree) and shape = evidence (shared by both layers, unique to one,
#'   or found only by the combined p-value). Results made before the
#'   directional method carry no direction and are coloured by adjusted p.
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
  view = c("scatter", "dual_volcano", "effect_pair", "quadrant", "dotplot"),
  top_n = 20L,
  label_features = NULL,
  p_cutoff = 0.05
) {
  assert_count(top_n, "top_n")
  assert_character(label_features, "label_features", allow_null = TRUE)
  assert_number(p_cutoff, "p_cutoff", lower = 0, upper = 1)
  view <- match.arg(view)
  df <- integration_result_from_bundle(bundle)
  method <- bundle$params$method
  if (is.null(method)) {
    stop("Bundle is missing `params$method`.")
  }

  if (view == "dual_volcano" && method != "concordance") {
    stop("`dual_volcano` view requires method = 'concordance'.")
  }
  if (view == "effect_pair" && method != "concordance") {
    stop("`effect_pair` view requires method = 'concordance'.")
  }
  if (view == "quadrant" && method != "concordance") {
    stop("`quadrant` view requires method = 'concordance'.")
  }
  if (view == "dotplot" && method != "active_pathways") {
    stop("`dotplot` view requires method = 'active_pathways'.")
  }

  if (nrow(df) == 0L) {
    return(empty_plot("No integration rows to plot."))
  }

  switch(view,
    scatter      = plot_integration_scatter(df, bundle, top_n, label_features, p_cutoff),
    dual_volcano = plot_integration_dual_volcano(df, bundle, top_n, label_features, p_cutoff),
    effect_pair  = plot_integration_effect_pair(df, bundle),
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


integration_axis_label <- function(bundle, side) {
  exps <- bundle$params$experiments
  if (is.null(exps) || length(exps) < 2L) {
    return(if (side == "a") "experiment A" else "experiment B")
  }
  exps[[if (side == "a") 1L else 2L]]
}

# The text a labelled point carries: the gene, or -- where a gene has
# several protein-gene pairs and so several points -- the pair's own id
# ("TP53 (P04637-2)"), so two points do not both read "TP53".
feature_point_label <- function(df) {
  sym <- df$feature_symbol
  shared <- !is.na(sym) & (duplicated(sym) | duplicated(sym, fromLast = TRUE))
  ifelse(shared, df$feature_id, sym)
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
  df$.label <- ifelse(df$feature_id %in% label_ids, feature_point_label(df), NA_character_)

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
      values = c(ns = omics_colors$ns, significant = omics_colors$up),
      name = NULL
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(p_cutoff), linetype = "dashed", color = omics_colors$ns
    ) +
    ggplot2::labs(
      title = title,
      subtitle = paste(bundle$params$experiments, collapse = " vs "),
      x = xlab,
      y = p_axis_label("adj_p_value")
    ) +
    theme_omics_labelled()

  p + add_repel_layer(df, x_aes, ".neglog10p", ".label")
}

plot_integration_dual_volcano <- function(df, bundle, top_n, label_features, p_cutoff) {
  df$.neglog10p <- -log10(pmax(df$p_value, .Machine$double.xmin))
  df$.quad <- ifelse(is.na(df$quadrant), "n/a", df$quadrant)
  # As in the effect-pair view: colour what both layers call a hit, when
  # the result says which those are.
  if (all(c("significant_a", "significant_b") %in% names(df))) {
    both <- df$significant_a %in% TRUE & df$significant_b %in% TRUE
    df$.quad[!both] <- "n/a"
    df <- df[order(both), , drop = FALSE]
  }
  label_ids <- pick_label_ids(df, top_n, label_features, p_col = "p_value")
  df$.label <- ifelse(df$feature_id %in% label_ids, feature_point_label(df), NA_character_)

  # Shared with the Shiny front end so the same comparison is tinted
  # identically on screen and in an exported report.
  quadrant_colors <- quadrant_palette()

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data$effect, y = .data$.neglog10p, color = .data$.quad)
  ) +
    ggplot2::geom_point(alpha = 0.8, size = 1.6, na.rm = TRUE) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = omics_colors$ns) +
    ggplot2::geom_hline(
      yintercept = -log10(p_cutoff), linetype = "dashed", color = omics_colors$ns
    ) +
    ggplot2::scale_color_manual(values = quadrant_colors, name = "quadrant",
                                na.value = omics_colors$ns) +
    ggplot2::labs(
      title = "Integration: dual volcano",
      subtitle = paste(bundle$params$experiments, collapse = " vs "),
      x = paste0("effect (", integration_axis_label(bundle, "a"),
                 " - ", integration_axis_label(bundle, "b"), ")"),
      y = "-log10(combined p)"
    ) +
    theme_omics_labelled()

  p + add_repel_layer(df, "effect", ".neglog10p", ".label")
}

# Effect against effect, one axis per layer. The dual volcano answers
# "how much do the two layers disagree?"; this answers "where does each
# feature sit in both?" -- concordant features fall on the diagonal, and
# the off-diagonal quadrants are the ones worth reading.
plot_integration_effect_pair <- function(df, bundle) {
  df <- integration_fill_effects(df)
  if (all(is.na(df$effect_a)) || all(is.na(df$effect_b))) {
    return(empty_plot(paste(
      "This integration result does not carry the per-layer effects.",
      "Re-run the integration to draw them.", sep = "\n")))
  }
  quad <- if ("quadrant" %in% names(df)) {
    ifelse(is.na(df$quadrant), "n/a", df$quadrant)
  } else {
    integration_derive_quadrant(df)
  }
  # Colour only what both layers call a hit; the rest of the cloud is
  # context. Every feature has *some* sign pair, so colouring them all
  # painted half of an unrelated background "concordant".
  both <- if (all(c("significant_a", "significant_b") %in% names(df))) {
    df$significant_a %in% TRUE & df$significant_b %in% TRUE
  } else {
    df$is_significant %in% TRUE
  }
  df$.quad <- ifelse(both, quad, "not significant in both")
  df <- df[order(both), , drop = FALSE]
  palette <- c(quadrant_palette(), `not significant in both` = omics_colors$ns)

  ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data$effect_a, y = .data$effect_b, color = .data$.quad)
  ) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         color = omics_colors$ns) +
    ggplot2::geom_hline(yintercept = 0, color = omics_colors$border) +
    ggplot2::geom_vline(xintercept = 0, color = omics_colors$border) +
    ggplot2::geom_point(alpha = 0.85, size = 2, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = palette, name = NULL,
                                na.value = omics_colors$ns) +
    ggplot2::labs(
      title = "Integration: effect pair",
      subtitle = paste(bundle$params$experiments, collapse = " vs "),
      x = paste0("effect (", integration_axis_label(bundle, "a"), ")"),
      y = paste0("effect (", integration_axis_label(bundle, "b"), ")")
    ) +
    theme_omics_labelled()
}

# The per-layer effects. Bundles written before the concordance table
# kept `effect_a` / `effect_b` store only their difference, from which
# the two cannot be recovered; those get NA (and the plot says so)
# rather than coordinates invented from the difference -- which put every
# point on the x axis at its difference and 0.
integration_fill_effects <- function(df) {
  if (all(c("effect_a", "effect_b") %in% names(df))) return(df)
  df$effect_a <- df$raw_a %||% rep(NA_real_, nrow(df))
  df$effect_b <- df$raw_b %||% rep(NA_real_, nrow(df))
  df
}

integration_derive_quadrant <- function(df) {
  dir_a <- if ("direction_a" %in% names(df)) df$direction_a
           else sign(df$effect_a %||% df$effect %||% 0)
  dir_b <- if ("direction_b" %in% names(df)) df$direction_b
           else sign(df$effect_b %||% 0)
  to_label <- function(s) ifelse(s > 0, "up", ifelse(s < 0, "down", "ns"))
  paste0(to_label(dir_a), "_", to_label(dir_b))
}

plot_integration_quadrant <- function(df, bundle) {
  # Among the features that are hits in both layers. Counted over every
  # feature, the bars reflected each backend's own direction labels at
  # its own thresholds, not the cutoffs this integration was run at.
  both <- if (all(c("significant_a", "significant_b") %in% names(df))) {
    df$significant_a %in% TRUE & df$significant_b %in% TRUE
  } else rep(TRUE, nrow(df))
  quads <- df$quadrant[both]
  quads <- quads[!is.na(quads)]
  levels_order <- c("up_up", "down_down", "up_down", "down_up")
  counts <- as.data.frame(table(factor(quads, levels = levels_order)),
                          stringsAsFactors = FALSE)
  names(counts) <- c("quadrant", "n")

  # Shared with the Shiny front end so the same comparison is tinted
  # identically on screen and in an exported report.
  quadrant_colors <- quadrant_palette()

  ggplot2::ggplot(
    counts,
    ggplot2::aes(x = .data$quadrant, y = .data$n, fill = .data$quadrant)
  ) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = quadrant_colors, guide = "none") +
    ggplot2::labs(
      title = "Integration: concordance quadrants",
      subtitle = sprintf("%s; features significant in both layers",
                         paste(bundle$params$experiments, collapse = " vs ")),
      x = NULL, y = "features"
    ) +
    theme_omics_labelled()
}

plot_integration_dotplot <- function(df, bundle, top_n) {
  df <- df[!is.na(df$adj_p_value), , drop = FALSE]
  df <- df[order(df$adj_p_value), , drop = FALSE]
  df <- utils::head(df, top_n)
  if (nrow(df) == 0L) {
    return(empty_plot("No pathways to plot."))
  }
  df$.label <- truncate_pathway_name(prettify_gene_set_name(df$feature_symbol))
  df$.label <- factor(df$.label, levels = unique(df$.label[order(-df$adj_p_value)]))
  # Bundles made before the directional rework kept the evidence class
  # in `direction` and had no pathway direction at all.
  directional <- "evidence" %in% names(df)
  df$.evidence <- if (directional) df$evidence else df$direction
  subtitle <- paste(bundle$params$experiments, collapse = " vs ")
  mm <- bundle$params$merge_method
  if (!is.null(mm)) {
    subtitle <- paste0(subtitle, " \u00B7 ",
                       if (mm %in% DIRECTIONAL_MERGE_METHODS)
                         "layers expected to agree in direction"
                       else "direction not used in the test")
  }

  shape_scale <- ggplot2::scale_shape_manual(
    values = c(shared = 16, unique = 1, combined = 17),
    labels = c(shared = "both layers", unique = "one layer",
               combined = "only combined"),
    na.value = 4, name = "found by")
  p <- if (directional) {
    df$.dir <- factor(ap_direction_label(df$direction),
                      levels = unname(AP_DIRECTION_LABELS))
    ggplot2::ggplot(df, ggplot2::aes(x = .data$effect, y = .data$.label,
                                     color = .data$.dir, shape = .data$.evidence)) +
      ggplot2::geom_point(size = 4, na.rm = TRUE) +
      ggplot2::scale_color_manual(
        values = stats::setNames(
          c(omics_colors$up, omics_colors$down, omics_colors$conc_up_down,
            omics_colors$ns),
          unname(AP_DIRECTION_LABELS)),
        drop = TRUE, name = "direction")
  } else {
    ggplot2::ggplot(df, ggplot2::aes(x = .data$effect, y = .data$.label,
                                     color = .data$adj_p_value, shape = .data$.evidence)) +
      ggplot2::geom_point(size = 4, na.rm = TRUE) +
      ggplot2::scale_color_gradient(low = omics_colors$up, high = omics_colors$ns,
                                    name = "adj p")
  }
  p + shape_scale +
    ggplot2::labs(
      title = "Integration: ActivePathways",
      subtitle = subtitle,
      x = "-log10(adj p)", y = NULL
    ) +
    theme_omics_labelled()
}

# A pathway's direction in words, for legends and tables.
AP_DIRECTION_LABELS <- c(up = "up in both layers", down = "down in both layers",
                         mixed = "mixed / layers disagree", none = "no direction")

ap_direction_label <- function(direction) {
  out <- unname(AP_DIRECTION_LABELS[direction])
  out[is.na(out)] <- AP_DIRECTION_LABELS[["none"]]
  out
}
