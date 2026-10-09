# Top-features heatmap. The function is overloaded:
#   * input is an `omics_input` -> heatmap of `top_n` most-variable features
#     across all samples;
#   * input is a `run_diff` analysis_bundle -> heatmap of the top features by
#     (adjusted) p-value, drawn from the original `omics_input` carried in
#     `bundle$input_info` ... except diff bundles don't carry the expression
#     matrix. So callers must also pass `input` for the bundle case.
#
# Renders via ComplexHeatmap when installed, otherwise a ggplot tile plot.
# `engine` picks one explicitly: the app asks for the ggplot, which takes
# its theme, fits a card's width and saves like every other figure.

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
#' Values are drawn on a comparable scale first: raw counts as
#' log2(CPM + 1) with library sizes from the whole matrix, linear
#' intensities, TPM and FPKM as log2(x + 1). With `scale = "row"` (the
#' default) each feature is then centred and scaled across the samples
#' shown, and the colours run blue (below the feature's mean) through
#' white to red (above it).
#'
#' With `group_by`, the samples are drawn group by group (in the order of
#' `group_levels`, or else the layer's reference group first and the rest
#' sorted), a gap between groups, under a bar coloured by group -- in the
#' colours the PCA gives the same groups ([group_palette()]). `highlight`
#' outlines a feature's row and names it in bold, also when there are too
#' many rows to name them all.
#'
#' When `ComplexHeatmap` is installed the function returns a
#' `ComplexHeatmap::Heatmap` unless `engine = "ggplot"`; otherwise it
#' falls back to a ggplot2 tile plot so the function is always usable
#' from a clean install.
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
#'   `input$meta_df` to surface as a column annotation. The ggplot
#'   heatmap draws one annotation bar: `group_by`, or else the first of
#'   these.
#' @param cluster_rows,cluster_cols Whether to cluster rows / columns.
#'   With `group_by` the columns are clustered within each group, so the
#'   groups stay together.
#' @param show_rownames,show_colnames If `NULL`, auto-decide: rows are
#'   named up to 80 of them (40 in the ggplot heatmap, where they would
#'   otherwise overlap at a screen's size), samples up to 40.
#' @param title Plot title.
#' @param group_by Optional column of `input$meta_df` to group and colour
#'   the samples by.
#' @param group_levels Optional groups of `group_by` to show, in order;
#'   samples in other groups are left out. Requires `group_by`.
#' @param highlight Optional `feature_id`s or symbols to mark.
#' @param engine `"auto"` (ComplexHeatmap when installed), `"ComplexHeatmap"`
#'   or `"ggplot"`.
#'
#' @return A `ComplexHeatmap::Heatmap` object when `ComplexHeatmap` is
#'   used, otherwise a `ggplot` tile plot.
#' @export
#' @family diff
#' @examples
#' \dontrun{
#'   # From an omics_input
#'   plot_heatmap(input, n_top = 30)
#'
#'   # From a diff bundle: its top 50 features, the two groups compared,
#'   # one gene marked
#'   b <- run_diff(input, ...)
#'   plot_heatmap(b, input = input, n_top = 50, group_by = "group",
#'                group_levels = c("Control", "TreatA"), highlight = "IL6",
#'                engine = "ggplot")
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
  title = NULL,
  group_by = NULL,
  group_levels = NULL,
  highlight = NULL,
  show_colnames = NULL,
  engine = c("auto", "ComplexHeatmap", "ggplot")
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
  assert_flag(show_colnames, "show_colnames", allow_null = TRUE)
  assert_string(title, "title", allow_null = TRUE, allow_empty = TRUE)
  assert_string(group_by, "group_by", allow_null = TRUE)
  assert_character(group_levels, "group_levels", allow_null = TRUE)
  assert_character(highlight, "highlight", allow_null = TRUE)
  scale <- match.arg(scale)
  engine <- match.arg(engine)
  if (identical(engine, "ComplexHeatmap") && !is_installed("ComplexHeatmap")) {
    stop("`engine = \"ComplexHeatmap\"` needs the ComplexHeatmap package, ",
         "which is not installed; use engine = \"ggplot\".", call. = FALSE)
  }
  if (!is.null(group_levels) && is.null(group_by)) {
    stop("`group_levels` needs `group_by`: the column whose groups they are.",
         call. = FALSE)
  }
  sel <- resolve_heatmap_selection(x, input = input, n_top = n_top,
                                    features = features)
  mat <- sel$mat
  meta <- sel$meta
  title <- title %||% sel$default_title

  # The samples to draw, group by group.
  groups <- NULL
  if (!is.null(group_by)) {
    if (is.null(meta) || !group_by %in% colnames(meta)) {
      stop("`group_by` not found in `meta_df`: ", group_by, call. = FALSE)
    }
    # Ordered and coloured over the column as it is (a factor's levels),
    # as the PCA and the boxplot order and colour it.
    g_raw <- as.data.frame(meta)[colnames(mat), group_by]
    g_all <- as.character(g_raw)
    reference <- design_reference(if (inherits(x, "omics_input")) x else input, group_by)
    lv <- if (is.null(group_levels)) group_order(g_raw, reference)
          else group_levels[group_levels %in% g_all]
    if (!length(lv)) {
      stop("None of `group_levels` is a group of `", group_by, "`: ",
           paste(group_levels, collapse = ", "), call. = FALSE)
    }
    keep <- !is.na(g_all) & g_all %in% lv
    ord <- order(match(g_all[keep], lv))
    mat <- mat[, which(keep)[ord], drop = FALSE]
    groups <- factor(g_all[keep][ord], levels = lv)
    palette <- group_colours(g_raw, reference)
  }

  # Scaled over the samples shown, so a gene's colours compare the
  # groups on the figure rather than groups left off it.
  mat <- switch(scale,
    row    = scale_rows(mat),
    column = scale(mat),
    none   = mat
  )

  hl <- if (is.null(highlight)) rep(FALSE, nrow(mat)) else
    rownames(mat) %in% highlight | sel$labels %in% highlight

  use_complex <- switch(engine,
    auto = is_installed("ComplexHeatmap"),
    ComplexHeatmap = TRUE,
    ggplot = FALSE)

  if (use_complex) {
    if (is.null(show_rownames)) show_rownames <- nrow(mat) <= 80L
    rownames(mat) <- sel$labels
    plot_heatmap_complex(
      mat = mat, meta = meta, annotation_cols = annotation_cols,
      cluster_rows = cluster_rows, cluster_cols = cluster_cols,
      show_rownames = show_rownames, title = title,
      groups = groups, group_by = group_by,
      palette = if (!is.null(groups)) palette, highlight = hl,
      show_colnames = show_colnames %||% (ncol(mat) <= 40L)
    )
  } else {
    # One annotation bar: the groups, or else the first annotation column.
    if (is.null(groups) && length(annotation_cols) && !is.null(meta) &&
        annotation_cols[[1L]] %in% colnames(meta)) {
      group_by <- annotation_cols[[1L]]
      g_raw <- as.data.frame(meta)[colnames(mat), group_by]
      reference <- design_reference(if (inherits(x, "omics_input")) x else input, group_by)
      groups <- factor(as.character(g_raw), levels = group_order(g_raw, reference))
      palette <- group_colours(g_raw, reference)
    }
    plot_heatmap_ggplot(
      mat = mat, labels = sel$labels, title = title,
      groups = groups, group_by = group_by,
      palette = if (!is.null(groups)) palette,
      cluster_rows = cluster_rows, cluster_cols = cluster_cols,
      show_rownames = show_rownames %||% (nrow(mat) <= 40L),
      show_colnames = show_colnames %||% (ncol(mat) <= 40L),
      highlight = hl, scaled = !identical(scale, "none")
    )
  }
}

# ---- internal helpers --------------------------------------------------

# The matrix comes back with feature ids as row names; `labels` are what
# to call each row (the symbol where there is one).
resolve_heatmap_selection <- function(x, input, n_top, features) {
  if (inherits(x, "omics_input")) {
    selected <- pick_features_by_variance(x$expr_mat, x$assay_type,
                                          n_top = n_top, features = features)
    mat <- coerce_to_continuous(x$expr_mat[selected, , drop = FALSE], x$assay_type,
                                lib_size = lib_sizes(x, colnames(x$expr_mat)))
    list(
      mat = mat,
      labels = map_feature_symbols(x$feature_df, selected),
      meta = x$meta_df,
      default_title = paste0("Top ", nrow(mat), " variable features")
    )
  } else if (is_analysis_bundle(x) &&
             identical(x$analysis_name, "run_diff")) {
    if (is.null(input) || !inherits(input, "omics_input")) {
      stop("`input` (an omics_input) is required when `x` is a diff bundle.")
    }
    res <- diff_result_from_bundle(x)
    if (is.null(res) || nrow(res) == 0L) {
      stop("Diff bundle does not contain a non-empty `results$diff_result_df`.")
    }
    selected <- pick_features_by_diff(res, n_top = n_top, features = features)
    keep <- intersect(selected, rownames(input$expr_mat))
    if (length(keep) == 0L) {
      stop("None of the top features from the diff bundle are present in `input$expr_mat`.")
    }
    # The library sizes of the whole matrix: counts-per-million over the
    # fifty genes drawn would rescale every sample by how much those
    # fifty happen to sum to.
    mat <- coerce_to_continuous(input$expr_mat[keep, , drop = FALSE], input$assay_type,
                                lib_size = lib_sizes(input, colnames(input$expr_mat)))
    list(
      mat = mat,
      labels = map_feature_symbols(input$feature_df, keep),
      meta = input$meta_df,
      default_title = paste0("Top ", length(keep), " features by adjusted p")
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

# A missing value stays missing (drawn grey by the ggplot heatmap) unless
# asked otherwise: filled with 0 it was drawn as exactly the feature's
# mean, a measurement nobody made.
scale_rows <- function(mat, keep_na = TRUE) {
  out <- t(scale(t(mat)))
  # A feature with one value, or the same value everywhere, has no
  # spread to scale by: its row is drawn at its mean.
  flat <- !is.finite(out) & !is.na(mat)
  out[flat] <- 0
  if (!keep_na) out[!is.finite(out)] <- 0
  out
}

plot_heatmap_complex <- function(mat, meta, annotation_cols, cluster_rows,
                                  cluster_cols, show_rownames, title,
                                  groups = NULL, group_by = NULL, palette = NULL,
                                  highlight = rep(FALSE, nrow(mat)),
                                  show_colnames = TRUE) {
  top_anno <- NULL
  anno_df <- NULL
  if (!is.null(meta) && !is.null(annotation_cols)) {
    cols <- intersect(annotation_cols, colnames(meta))
    if (length(cols) > 0L) {
      anno_df <- as.data.frame(meta)[colnames(mat), cols, drop = FALSE]
    }
  }
  anno_col <- list()
  # Other grouping columns in the group colours too: ComplexHeatmap's own
  # are random, and as likely as not a red or a blue.
  for (col in names(anno_df)) {
    v <- anno_df[[col]]
    if ((is.character(v) || is.factor(v)) && any(!is.na(v))) anno_col[[col]] <- group_colours(v)
  }
  if (!is.null(groups)) {
    if (is.null(anno_df)) anno_df <- data.frame(row.names = colnames(mat))
    anno_df[[group_by]] <- groups
    anno_col[[group_by]] <- palette[levels(groups)]
  }
  if (!is.null(anno_df)) {
    top_anno <- ComplexHeatmap::HeatmapAnnotation(df = anno_df, col = anno_col)
  }
  col_fn <- if (is_installed("circlize")) {
    circlize::colorRamp2(c(-2, 0, 2), c(omics_colors$down, "white", omics_colors$up))
  } else {
    NULL
  }
  ComplexHeatmap::Heatmap(
    # ComplexHeatmap clusters with missing values only when there are
    # few of them; a missing cell is drawn at the feature's mean here.
    ifelse(is.na(mat), 0, mat),
    name = "Z-score",
    top_annotation = top_anno,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_cols,
    column_split = groups,
    show_row_names = show_rownames,
    row_names_gp = grid::gpar(fontface = ifelse(highlight, "bold", "plain")),
    show_column_names = show_colnames,
    column_title = title,
    col = col_fn
  )
}

# The ggplot heatmap. Rows run from the top in cluster order, samples left
# to right group by group with a gap between the groups, the group bar
# above them. The row names are drawn as text to the left of the tiles,
# not as axis labels, so one of them can be bold: the highlighted row is
# outlined and named in bold even when the rest go unnamed.
plot_heatmap_ggplot <- function(mat, labels = rownames(mat), title = NULL,
                                groups = NULL, group_by = NULL, palette = NULL,
                                cluster_rows = FALSE, cluster_cols = FALSE,
                                show_rownames = TRUE, show_colnames = TRUE,
                                highlight = rep(FALSE, nrow(mat)),
                                scaled = TRUE) {
  n_row <- nrow(mat)
  n_col <- ncol(mat)
  filled <- mat
  filled[!is.finite(filled)] <- 0

  row_ord <- seq_len(n_row)
  if (cluster_rows && n_row >= 3L) {
    row_ord <- stats::hclust(stats::dist(filled))$order
  }
  col_ord <- seq_len(n_col)
  if (cluster_cols) {
    # Within each group, so the groups stay together.
    blocks <- if (is.null(groups)) list(seq_len(n_col)) else
      split(seq_len(n_col), groups, drop = TRUE)
    col_ord <- unlist(lapply(blocks, function(idx) {
      if (length(idx) < 3L) return(idx)
      idx[stats::hclust(stats::dist(t(filled[, idx, drop = FALSE])))$order]
    }), use.names = FALSE)
  }
  mat <- mat[row_ord, col_ord, drop = FALSE]
  labels <- labels[row_ord]
  highlight <- highlight[row_ord]
  if (!is.null(groups)) groups <- groups[col_ord]

  # x: one unit per sample, and a gap of a third of a sample between
  # groups. y: the first row at the top.
  gap <- 0.35
  xpos <- seq_len(n_col) +
    if (is.null(groups)) 0 else (as.integer(groups) - 1L) * gap
  ypos <- rev(seq_len(n_row))
  tiles <- data.frame(
    x = rep(xpos, each = n_row),
    y = rep(ypos, times = n_col),
    value = as.numeric(mat)
  )
  x_lo <- min(xpos) - 0.5
  x_hi <- max(xpos) + 0.5

  p <- ggplot2::ggplot(tiles, ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$value), width = 1, height = 1)
  p <- p + if (scaled) {
    ggplot2::scale_fill_gradient2(
      low = omics_colors$down, mid = "white", high = omics_colors$up,
      midpoint = 0, limits = c(-2, 2), breaks = -2:2, oob = scales::squish,
      na.value = "#D9DDE3", name = "z-score",
      guide = ggplot2::guide_colourbar(order = 1))
  } else {
    ggplot2::scale_fill_gradient(low = "#F3F5F8", high = omics_colors$fg_dark,
                                 na.value = "#D9DDE3", name = "value",
                                 guide = ggplot2::guide_colourbar(order = 1))
  }

  # The group bar, above the tiles: a constant fill per tile (the fill
  # scale is the z-score's), and its legend from an invisible point layer
  # mapped to colour.
  y_top <- n_row + 0.5
  label_y <- NULL
  if (!is.null(groups)) {
    bar_h <- max(0.8, n_row * 0.03)
    bar_y <- n_row + 0.5 + bar_h * 0.4 + bar_h / 2
    bar <- data.frame(x = xpos, y = bar_y, group = groups)
    p <- p +
      ggplot2::geom_tile(data = bar, width = 1, height = bar_h,
                         fill = unname(palette[as.character(groups)])) +
      ggplot2::geom_point(data = bar, ggplot2::aes(colour = .data$group), alpha = 0) +
      ggplot2::scale_colour_manual(
        values = palette, breaks = levels(groups), name = group_by,
        labels = function(x) wrap_label(x, width = 18L, max_lines = 3L),
        guide = ggplot2::guide_legend(ncol = 1, order = 2,
                                      override.aes = list(alpha = 1, shape = 15, size = 4)))
    y_top <- bar_y + bar_h / 2
    label_y <- data.frame(y = bar_y, text = group_by, face = "italic")
  }

  # Row names to the left of the tiles. Long ones are cut: a 40-character
  # identifier took half a phone's width.
  short <- ifelse(nchar(labels) > 20L, paste0(substr(labels, 1L, 19L), "…"), labels)
  named <- if (show_rownames) rep(TRUE, n_row) else highlight
  row_txt <- data.frame(y = ypos[named], text = short[named],
                        face = ifelse(highlight[named], "bold", "plain"))
  txt <- rbind(row_txt, label_y)
  if (nrow(txt)) {
    p <- p + ggplot2::geom_text(
      # A trailing space keeps the names off the tiles by the same few
      # px at any column width.
      data = txt, ggplot2::aes(x = x_lo, y = .data$y, label = paste0(.data$text, "  "),
                               fontface = .data$face),
      hjust = 1, size = if (n_row > 30L) 2.7 else 3, colour = omics_colors$fg_dark,
      inherit.aes = FALSE)
  }
  if (any(highlight)) {
    hy <- ypos[highlight]
    p <- p + ggplot2::annotate("rect", xmin = x_lo, xmax = x_hi,
                               ymin = hy - 0.5, ymax = hy + 0.5,
                               fill = NA, colour = omics_colors$fg_dark, linewidth = 0.7)
  }
  # Room on the left for the longest name drawn (about 5.4 pt a character
  # at this size).
  left_pt <- if (nrow(txt)) (max(nchar(txt$text)) + 2) * 5.4 + 6 else 5.5

  p +
    ggplot2::scale_x_continuous(
      breaks = if (show_colnames) xpos else NULL,
      labels = if (show_colnames) colnames(mat) else NULL,
      expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = NULL, expand = c(0, 0)) +
    ggplot2::coord_cartesian(xlim = c(x_lo, x_hi), ylim = c(0.5, y_top), clip = "off") +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    theme_omicsCore() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5,
                                          size = ggplot2::rel(0.7)),
      axis.ticks = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(5.5, 5.5, 5.5, left_pt)
    )
}
