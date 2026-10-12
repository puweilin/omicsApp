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
#' * `"effect_pair"` -- concordance-only. Each layer's effect against the
#'   other's on equal axes, with the diagonal of perfect agreement. The
#'   features that are hits in both layers are coloured by which way each
#'   layer moved (the legend counts them), the rest stay a faint grey
#'   background, and the top hits (at most `min(top_n, 8)`, ranked by
#'   the combined adjusted p) are ringed and named.
#' * `"top_hits"` -- concordance-only. The top `top_n` hits in both
#'   layers, ranked by the combined adjusted p, one row per feature with
#'   a dot for each layer's effect joined by a line: the shorter the
#'   line, the closer the two layers agree.
#' * `"quadrant"` -- concordance-only. Bar count of the four
#'   `(direction_a, direction_b)` sign quadrants.
#' * `"dotplot"` -- active_pathways-only. Dotplot of top pathways, with
#'   colour = direction (up in both layers, down in both, or layers
#'   disagree) and shape = evidence (shared by both layers, unique to one,
#'   or found only by the combined p-value). Results made before the
#'   directional method carry no direction and are coloured by adjusted p.
#' * `"dual_volcano"` -- deprecated, kept so scripts exported by earlier
#'   versions still run; it warns. It plotted the difference of the two
#'   effects against the combined p, which put the features that agree
#'   best at the centre, where a volcano reader looks for "no change".
#'   The difference is each point's distance from the diagonal in
#'   `"effect_pair"`.
#'
#' @param bundle An [`analysis_bundle`][is_analysis_bundle()] produced by
#'   [run_integration()].
#' @param view One of `"scatter"`, `"effect_pair"`, `"top_hits"`,
#'   `"quadrant"`, `"dotplot"` (or the deprecated `"dual_volcano"`).
#' @param top_n Number of features / pathways to label or display. The
#'   scatter view names at most eight, and only significant features.
#' @param label_features Optional character vector of `feature_symbol`
#'   values to force-label (scatter / effect_pair views).
#' @param p_cutoff Significance cutoff used for highlight color in the
#'   scatter view.
#'
#' @return A `ggplot` object.
#' @export
#' @family integration
plot_integration <- function(
  bundle,
  view = c("scatter", "effect_pair", "top_hits", "quadrant", "dotplot",
           "dual_volcano"),
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

  if (view %in% c("dual_volcano", "effect_pair", "top_hits", "quadrant") &&
      method != "concordance") {
    stop("`", view, "` view requires method = 'concordance'.")
  }
  if (view == "dual_volcano") {
    warning(structure(class = c("deprecatedWarning", "warning", "condition"), list(
      message = paste(
        "`view = \"dual_volcano\"` is deprecated and will be removed.",
        "Use `view = \"effect_pair\"` (the difference of the two effects is",
        "each point's distance from the diagonal) or `view = \"top_hits\"`."),
      call = NULL)))
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
    effect_pair  = plot_integration_effect_pair(df, bundle, top_n, label_features),
    top_hits     = plot_integration_top_hits(df, bundle, top_n),
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

# The cap on names in the scatter views. Twenty, the old default, drew
# twenty leader lines into the densest corner of the plot and crossed
# most of them; past about eight the names stop being readable at the
# size the app draws this card.
SCATTER_MAX_LABELS <- 8L

plot_integration_scatter <- function(df, bundle, top_n, label_features, p_cutoff) {
  method <- bundle$params$method
  df$.neglog10p <- -log10(pmax(df$adj_p_value, .Machine$double.xmin))
  sig <- !is.na(df$is_significant) & df$is_significant
  df$.sig <- factor(ifelse(sig, "significant", "ns"), levels = c("ns", "significant"))

  # Only significant features are named, the strongest first: a name on
  # a grey point invites reading it as a finding. Names asked for by the
  # caller are added whatever their p.
  label_ids <- pick_label_ids(df[sig, , drop = FALSE], min(top_n, SCATTER_MAX_LABELS),
                              label_features = NULL)
  forced <- if (is.null(label_features)) character(0)
            else df$feature_id[df$feature_symbol %in% label_features]
  df$.label <- ifelse(df$feature_id %in% c(label_ids, forced),
                      feature_point_label(df), NA_character_)

  if (method == "correlation") {
    x_aes <- "effect"
    type <- first_or_na(df$effect_type)
    # Short enough for a half-width card: the longer "... across paired
    # samples" ran off both ends of the panel.
    xlab <- switch(type, spearman_r = "Spearman correlation (paired samples)",
                   pearson_r = "Pearson correlation (paired samples)",
                   paste0("correlation (", type, ")"))
    title <- "Correlation between layers, gene by gene"
  } else if (method == "concordance") {
    x_aes <- "effect"
    xlab <- paste0("effect difference (", integration_axis_label(bundle, "a"),
                   " - ", integration_axis_label(bundle, "b"), ")")
    title <- "Difference in fold change between layers"
  } else {
    df$.rank <- seq_len(nrow(df))
    x_aes <- ".rank"
    xlab <- "pathway rank"
    title <- "Pathways ranked across both layers"
  }
  # The background first and faint, the hits on top of it.
  df <- df[order(sig), , drop = FALSE]
  n_sig <- sum(sig)

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data[[x_aes]], y = .data$.neglog10p, color = .data$.sig)
  ) +
    ggplot2::geom_point(ggplot2::aes(size = .data$.sig, alpha = .data$.sig),
                        stroke = 0, na.rm = TRUE) +
    ggplot2::scale_size_manual(values = c(ns = 1, significant = 1.9), guide = "none") +
    ggplot2::scale_alpha_manual(values = c(ns = 0.35, significant = 0.9), guide = "none") +
    ggplot2::scale_color_manual(
      values = c(ns = omics_colors$ns, significant = omics_colors$up),
      labels = c(ns = "not significant",
                 significant = sprintf("significant (%s)", format(n_sig, big.mark = ","))),
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
    theme_omics_labelled() +
    # One short key line above the panel: beside it, it took the width the
    # x-axis title needed on a half-width card.
    ggplot2::theme(legend.position = "top", legend.justification = "left",
                   legend.location = "plot")

  p + add_repel_layer(df, x_aes, ".neglog10p", ".label",
                      min.segment.length = 0, box.padding = 0.4,
                      segment.color = omics_colors$ns, seed = 1L)
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
      title = "Difference between layers, by significance",
      subtitle = paste(bundle$params$experiments, collapse = " vs "),
      x = paste0("effect (", integration_axis_label(bundle, "a"),
                 " - ", integration_axis_label(bundle, "b"), ")"),
      y = "-log10(combined p)"
    ) +
    theme_omics_labelled()

  p + add_repel_layer(df, "effect", ".neglog10p", ".label")
}

# Effect against effect, one axis per layer, on equal axes: features the
# layers agree on fall on the diagonal, and a feature's distance from it
# is how far the layers disagree -- what the old dual volcano put on its
# x axis. Only the features both layers call a hit are coloured; the
# rest of the cloud is context, drawn faint and underneath.
plot_integration_effect_pair <- function(df, bundle, top_n = 10L,
                                         label_features = NULL) {
  df <- integration_fill_effects(df)
  if (all(is.na(df$effect_a)) || all(is.na(df$effect_b))) {
    return(empty_plot(NO_LAYER_EFFECTS_MSG))
  }
  quad <- if ("quadrant" %in% names(df)) {
    ifelse(is.na(df$quadrant), "n/a", df$quadrant)
  } else {
    integration_derive_quadrant(df)
  }
  # Every feature has *some* sign pair, so colouring them all painted
  # half of an unrelated background "concordant".
  both <- integration_both_hits(df)
  classes <- integration_hit_classes(bundle)
  df$.class <- ifelse(both & quad %in% names(classes), quad, "background")
  df$.class <- factor(df$.class, levels = c(names(classes), "background"))
  n <- table(df$.class)
  shown <- names(n)[n > 0 | names(n) == "background"]
  legend_labels <- sprintf("%s (%s)", c(classes, background = "not a hit in both")[shown],
                           vapply(as.integer(n[shown]), format, character(1),
                                  big.mark = ","))
  palette <- c(quadrant_palette()[names(classes)], background = omics_colors$ns)
  df$.bg <- df$.class == "background"

  # The top hits by the same ranking as the result table, so the names
  # here are the first rows there. At most eight: they crowd one corner,
  # and past that the labels need leader lines longer than the plot.
  top <- integration_top_hit_rows(df, min(top_n, 8L))
  ring <- df[unique(c(top, which(df$feature_symbol %in% label_features))), , drop = FALSE]
  ring$.label <- feature_point_label(df)[match(ring$feature_id, df$feature_id)]
  df <- df[order(!df$.bg), , drop = FALSE]

  lim <- max(abs(c(df$effect_a, df$effect_b)), na.rm = TRUE)
  lim <- if (is.finite(lim) && lim > 0) lim * 1.05 else 1
  lab <- integration_effect_label(bundle)

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data$effect_a, y = .data$effect_b, color = .data$.class)
  ) +
    ggplot2::geom_hline(yintercept = 0, color = omics_colors$border) +
    ggplot2::geom_vline(xintercept = 0, color = omics_colors$border) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         color = omics_colors$ns) +
    ggplot2::geom_point(ggplot2::aes(size = .data$.bg, alpha = .data$.bg),
                        stroke = 0, na.rm = TRUE) +
    ggplot2::scale_size_manual(values = c(`FALSE` = 1.9, `TRUE` = 1),
                               guide = "none") +
    ggplot2::scale_alpha_manual(values = c(`FALSE` = 0.9, `TRUE` = 0.35),
                                guide = "none") +
    ggplot2::scale_color_manual(values = palette, breaks = shown,
                                labels = legend_labels, name = NULL,
                                drop = FALSE) +
    ggplot2::guides(color = ggplot2::guide_legend(
      ncol = 1, override.aes = list(size = 2.6, alpha = 1))) +
    ggplot2::coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    ggplot2::labs(
      title = "Fold change in each layer",
      # The axes name the layers; a phone-width card has no room to
      # repeat them here.
      subtitle = if (nrow(ring)) "ringed: top hits in both layers"
                 else paste(bundle$params$experiments, collapse = " vs "),
      x = paste0(lab, " (", integration_axis_label(bundle, "a"), ")"),
      y = paste0(lab, " (", integration_axis_label(bundle, "b"), ")")
    ) +
    theme_omics_labelled()

  if (nrow(ring)) {
    # Labels move away from the diagonal, into the corners where the
    # layers disagree: for most results those are the emptiest part of
    # the plot, and a hit off the diagonal is pushed further out still.
    away <- ifelse(ring$effect_b >= ring$effect_a, 1, -1)
    p <- p +
      ggplot2::geom_point(data = ring, shape = 21, size = 3.2, stroke = 0.8,
                          color = omics_colors$fg_dark, fill = NA,
                          inherit.aes = TRUE, show.legend = FALSE) +
      add_repel_layer(ring, "effect_a", "effect_b", ".label",
                      nudge_x = -away * lim * 0.3, nudge_y = away * lim * 0.3,
                      min.segment.length = 0, box.padding = 0.35,
                      segment.color = omics_colors$ns, seed = 1L)
  }
  p
}

# The top hits as a ranked list: one row per feature, a dot for each
# layer's effect and a line between them. A scatter shows where the hits
# sit; this is where a reader gets their names and can compare a handful
# side by side -- the shorter the line, the closer the layers agree, and
# a line crossing zero is a feature the layers disagree on.
plot_integration_top_hits <- function(df, bundle, top_n = 15L) {
  df <- integration_fill_effects(df)
  if (all(is.na(df$effect_a)) || all(is.na(df$effect_b))) {
    return(empty_plot(NO_LAYER_EFFECTS_MSG))
  }
  rows <- integration_top_hit_rows(df, top_n)
  if (!length(rows)) {
    return(empty_plot("No feature is a hit in both layers at these cutoffs."))
  }
  top <- df[rows, , drop = FALSE]
  top$.label <- make.unique(feature_point_label(df)[rows], sep = " ")
  top$.label <- factor(top$.label, levels = rev(top$.label))

  exps <- c(integration_axis_label(bundle, "a"), integration_axis_label(bundle, "b"))
  if (identical(exps[[1L]], exps[[2L]])) exps <- paste(exps, c("(A)", "(B)"))
  long <- data.frame(
    .label = rep(top$.label, 2L),
    layer = factor(rep(exps, each = nrow(top)), levels = exps),
    effect = c(top$effect_a, top$effect_b)
  )
  layer_colors <- stats::setNames(c(omics_colors$layer_a, omics_colors$layer_b), exps)

  ggplot2::ggplot(long, ggplot2::aes(x = .data$effect, y = .data$.label)) +
    ggplot2::geom_vline(xintercept = 0, color = omics_colors$border) +
    ggplot2::geom_segment(
      data = top,
      ggplot2::aes(x = .data$effect_a, xend = .data$effect_b,
                   y = .data$.label, yend = .data$.label),
      color = "#C9CED6", linewidth = 1.1, na.rm = TRUE
    ) +
    ggplot2::geom_point(ggplot2::aes(color = .data$layer, shape = .data$layer),
                        size = 2.6, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = layer_colors, name = NULL) +
    ggplot2::scale_shape_manual(values = stats::setNames(c(16, 15), exps), name = NULL) +
    ggplot2::labs(
      title = "Top hits in both layers",
      subtitle = if (nrow(top) > 1L)
        "ranked by combined adjusted p\nshort line = layers agree"
      else "ranked by combined adjusted p",
      x = integration_effect_label(bundle), y = NULL
    ) +
    theme_omics_labelled() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
}

NO_LAYER_EFFECTS_MSG <- paste(
  "This integration result does not carry the per-layer effects.",
  "Re-run the integration to draw them.", sep = "\n")

# The features both layers call a hit. Results from before the per-layer
# columns only know the combined call.
integration_both_hits <- function(df) {
  if (all(c("significant_a", "significant_b") %in% names(df))) {
    df$significant_a %in% TRUE & df$significant_b %in% TRUE
  } else {
    df$is_significant %in% TRUE
  }
}

# Row numbers of the top `n` hits in both layers, in the order the result
# table lists them (combined adjusted p, then raw p).
integration_top_hit_rows <- function(df, n) {
  rows <- which(integration_both_hits(df) & !is.na(df$adj_p_value) &
                  !is.na(df$effect_a) & !is.na(df$effect_b))
  rows <- rows[order(df$adj_p_value[rows], df$p_value[rows])]
  utils::head(rows, n)
}

# Plain names for the four sign quadrants, with the layers' own names for
# the two where they disagree ("proteomics up, rnaseq down").
integration_hit_classes <- function(bundle) {
  a <- integration_axis_label(bundle, "a")
  b <- integration_axis_label(bundle, "b")
  c(up_up = "up in both", down_down = "down in both",
    up_down = sprintf("%s up, %s down", a, b),
    down_up = sprintf("%s down, %s up", a, b))
}

# The name of the per-layer effect: log2FC unless either layer's result
# came from a test whose effect is something else (a slope, an F).
# Results that do not record how the layers were tested are group
# comparisons -- the only kind the concordance view offered.
integration_effect_label <- function(bundle) {
  types <- vapply(bundle$params$diff_params %||% list(),
                  function(p) p$analysis_type %||% "group", character(1))
  if (any(types != "group")) "effect" else "log2FC"
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
      title = "Direction of change in both layers",
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
  # Wrapped, not cut: two pathways sharing their first words must stay
  # told apart, and the hover card gives the whole name anyway.
  df$.label <- wrap_pathway_name(prettify_gene_set_name(df$feature_symbol), width = 28L)
  df$.label <- factor(df$.label, levels = unique(df$.label[order(-df$adj_p_value)]))
  # Bundles made before the directional rework kept the evidence class
  # in `direction` and had no pathway direction at all.
  directional <- "evidence" %in% names(df)
  df$.evidence <- if (directional) df$evidence else df$direction
  subtitle <- paste(bundle$params$experiments, collapse = " vs ")
  mm <- bundle$params$merge_method
  if (!is.null(mm)) {
    subtitle <- paste0(subtitle, "\n",
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
      title = "Pathways from both layers combined",
      subtitle = subtitle,
      x = p_axis_label("adj_p_value"), y = NULL
    ) +
    ggplot2::scale_x_continuous(limits = c(0, NA),
                                expand = ggplot2::expansion(mult = c(0, 0.08))) +
    theme_omics_labelled() +
    # Keys under the panel: beside it, with pathway names on the left,
    # they left the points about 10 px of a 540 px card.
    ggplot2::guides(colour = ggplot2::guide_legend(ncol = 2, title.position = "top", order = 1),
                    shape = ggplot2::guide_legend(ncol = 2, title.position = "top", order = 2)) +
    ggplot2::theme(legend.position = "bottom", legend.box = "vertical",
                   legend.justification = "left", legend.location = "plot",
                   plot.title.position = "plot")
}

# A pathway's direction in words, for legends and tables.
AP_DIRECTION_LABELS <- c(up = "up in both layers", down = "down in both layers",
                         mixed = "mixed / layers disagree", none = "no direction")

ap_direction_label <- function(direction) {
  out <- unname(AP_DIRECTION_LABELS[direction])
  out[is.na(out)] <- AP_DIRECTION_LABELS[["none"]]
  out
}
