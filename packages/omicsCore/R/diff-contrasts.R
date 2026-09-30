# Several contrasts in one run_diff() bundle.
#
# run_diff(case_group = c("B", "C", "D"), control_group = "A") fits every
# group into one model and stacks the contrasts in `diff_result_df` under
# their `comparison` labels. Everything downstream of a diff -- the
# volcano, enrichment, integration, the report -- reads one contrast, so
# these helpers are the seam: list the contrasts, take one out as an
# ordinary single-contrast bundle, and summarise all of them side by side.

#' The comparisons a diff bundle holds
#'
#' @param bundle An [`analysis_bundle`][is_analysis_bundle()] from
#'   [run_diff()].
#' @return Character vector of `comparison` labels, in the order they were
#'   fitted (for example `c("B_vs_A", "C_vs_A")`).
#' @export
#' @family diff
diff_comparisons <- function(bundle) {
  assert_diff_bundle(bundle)
  df <- bundle$results$diff_result_df
  from_params <- as.character(bundle$params$comparison %||% character(0))
  from_rows <- unique(as.character(stats::na.omit(df$comparison)))
  out <- c(intersect(from_params, from_rows), setdiff(from_rows, from_params))
  if (!length(out)) out <- from_params
  out
}

#' Take one contrast out of a multi-contrast diff bundle
#'
#' Returns a bundle shaped exactly like a single-contrast [run_diff()]
#' result -- `diff_result_df` holds the rows of that comparison only,
#' `params$comparison` and `params$case_group` name it -- so every function
#' that takes a diff bundle (plots, [run_enrichment()],
#' [run_integration()], the report) can be handed it unchanged.
#'
#' @param bundle An [`analysis_bundle`][is_analysis_bundle()] from
#'   [run_diff()].
#' @param comparison One of [diff_comparisons()]`(bundle)`. May be omitted
#'   when the bundle holds a single comparison.
#' @return An [`analysis_bundle`][is_analysis_bundle()].
#' @export
#' @family diff
select_comparison <- function(bundle, comparison = NULL) {
  assert_diff_bundle(bundle)
  assert_string(comparison, "comparison", allow_null = TRUE)
  available <- diff_comparisons(bundle)
  if (is.null(comparison)) {
    if (length(available) > 1L) {
      stop("The bundle holds ", length(available), " comparisons (",
           paste(available, collapse = ", "), "); name one.", call. = FALSE)
    }
    return(bundle)
  }
  if (!comparison %in% available) {
    stop("Comparison '", comparison, "' is not in the bundle. Available: ",
         paste(available, collapse = ", "), ".", call. = FALSE)
  }
  if (length(available) == 1L) return(bundle)

  idx <- match(comparison, available)
  out <- bundle
  df <- bundle$results$diff_result_df
  out$results$diff_result_df <- df[df$comparison %in% comparison, , drop = FALSE]
  rownames(out$results$diff_result_df) <- NULL
  raw <- bundle$results$diff_raw_df
  if (is.data.frame(raw) && "comparison" %in% names(raw)) {
    raw <- raw[raw$comparison %in% comparison, , drop = FALSE]
    raw$comparison <- NULL
    rownames(raw) <- NULL
    out$results$diff_raw_df <- raw
  }
  obj <- bundle$results$diff_object
  # A per-pair backend (t-test, lm) leaves one model per comparison; a
  # shared-fit backend leaves the one model all of them came from, which
  # stays as it is.
  if (is.list(obj) && !is.object(obj) && length(obj) == length(available)) {
    out$results$diff_object <- obj[[idx]]
  }
  out$params$comparison <- comparison
  cases <- bundle$params$case_group
  ct <- bundle$params$contrast_table
  if (is.data.frame(ct) && comparison %in% ct$comparison) {
    # Contrasts given as expressions: a pair still has a case and a
    # control; a weighted contrast has neither, and says so.
    row <- ct[match(comparison, ct$comparison), , drop = FALSE]
    out$params$contrasts <- row$spec
    out$params$all_contrasts <- ct$spec
    out$params$case_group <- if (is.na(row$case)) NULL else row$case
    out$params$control_group <- if (is.na(row$control)) NULL else row$control
  } else if (length(cases) == length(available)) {
    out$params$case_group <- cases[[idx]]
  }
  # What the contrast was fitted alongside. Repeating it on another layer
  # (as Integration does) needs the same groups in the model, not just
  # this pair.
  out$params$all_comparisons <- available
  out$params$all_case_groups <- cases
  out
}

#' Hit counts for every contrast in a diff bundle
#'
#' One row per comparison: how many features pass the thresholds, split
#' by direction. With several treatment groups against one control this
#' is the table to read first -- which treatment moved the most, and in
#' which direction -- before opening any single volcano.
#'
#' @param bundle An [`analysis_bundle`][is_analysis_bundle()] from
#'   [run_diff()].
#' @inheritParams filter_diff_results
#' @return A `data.frame` with `comparison`, `n_up`, `n_down`, `n_hits`,
#'   `n_tested`, and `n_shared` (hits also found in at least one other
#'   comparison).
#' @export
#' @family diff
summarize_diff_contrasts <- function(
  bundle,
  p_cutoff = 0.05,
  p_preference = c("adjusted", "raw"),
  effect_cutoff = NULL
) {
  assert_diff_bundle(bundle)
  p_preference <- match.arg(p_preference)
  df <- bundle$results$diff_result_df
  comps <- diff_comparisons(bundle)
  hits <- lapply(comps, function(cmp) {
    sub <- df[df$comparison %in% cmp, , drop = FALSE]
    list(
      tested = nrow(sub),
      sig = filter_diff_results(sub, p_cutoff = p_cutoff,
                                p_preference = p_preference,
                                effect_cutoff = effect_cutoff)
    )
  })
  ids <- lapply(hits, function(h) unique(h$sig$feature_id))
  shared <- vapply(seq_along(ids), function(i) {
    others <- unique(unlist(ids[-i]))
    sum(ids[[i]] %in% others)
  }, integer(1))
  out <- data.frame(
    comparison = comps,
    n_up = vapply(hits, function(h) sum(h$sig$direction %in% c("up", "positive")),
                  integer(1)),
    n_down = vapply(hits, function(h) sum(h$sig$direction %in% c("down", "negative")),
                    integer(1)),
    n_hits = vapply(hits, function(h) nrow(h$sig), integer(1)),
    n_tested = vapply(hits, function(h) as.integer(h$tested), integer(1)),
    n_shared = shared,
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  out
}

assert_diff_bundle <- function(bundle) {
  if (!is_analysis_bundle(bundle) || !identical(bundle$analysis_name, "run_diff")) {
    stop("`bundle` must be an analysis_bundle from run_diff().", call. = FALSE)
  }
  if (!is.data.frame(bundle$results$diff_result_df)) {
    stop("`bundle` is missing `results$diff_result_df`.", call. = FALSE)
  }
  invisible(bundle)
}

#' Up / down hit counts per contrast
#'
#' A diverging bar chart of [summarize_diff_contrasts()]: one row per
#' comparison, up-regulated hits to the right, down-regulated to the left.
#'
#' @inheritParams summarize_diff_contrasts
#' @return A `ggplot` object.
#' @export
#' @family diff
plot_diff_contrasts <- function(
  bundle,
  p_cutoff = 0.05,
  p_preference = c("adjusted", "raw"),
  effect_cutoff = NULL
) {
  p_preference <- match.arg(p_preference)
  s <- summarize_diff_contrasts(bundle, p_cutoff = p_cutoff,
                                p_preference = p_preference,
                                effect_cutoff = effect_cutoff)
  pretty <- gsub("_vs_", " vs ", s$comparison, fixed = TRUE)
  long <- data.frame(
    comparison = factor(rep(pretty, 2L), levels = rev(pretty)),
    direction = factor(rep(c("up", "down"), each = nrow(s)),
                       levels = c("up", "down")),
    n = c(s$n_up, -s$n_down),
    label = c(s$n_up, s$n_down),
    stringsAsFactors = FALSE
  )
  ggplot2::ggplot(long, ggplot2::aes(x = .data$n, y = .data$comparison,
                                     fill = .data$direction)) +
    ggplot2::geom_col(width = 0.65) +
    ggplot2::geom_vline(xintercept = 0, color = omics_colors$ns) +
    ggplot2::geom_text(
      data = long[long$label > 0, , drop = FALSE],
      ggplot2::aes(label = .data$label,
                   hjust = ifelse(.data$n >= 0, -0.2, 1.2)),
      size = 3.4, color = "#333333"
    ) +
    ggplot2::scale_fill_manual(values = c(up = omics_colors$up,
                                          down = omics_colors$down),
                               name = NULL) +
    ggplot2::scale_x_continuous(
      labels = function(x) abs(x),
      expand = ggplot2::expansion(mult = 0.15)
    ) +
    ggplot2::labs(
      title = "Hits per comparison",
      subtitle = sprintf("%s p < %s%s",
                         if (p_preference == "adjusted") "adjusted" else "raw",
                         format(p_cutoff),
                         if (is.null(effect_cutoff)) "" else
                           sprintf(", |effect| ≥ %s", format(effect_cutoff))),
      x = "features (down ← → up)", y = NULL
    ) +
    theme_omics_labelled()
}

#' The hits of every comparison, as sets
#'
#' @inheritParams summarize_diff_contrasts
#' @param direction `"any"` (default), `"up"` or `"down"`: which hits.
#' @return A named list of character vectors of `feature_id`, one per
#'   comparison, in [diff_comparisons()] order.
#' @export
#' @family diff
diff_hit_sets <- function(
  bundle,
  p_cutoff = 0.05,
  p_preference = c("adjusted", "raw"),
  effect_cutoff = NULL,
  direction = c("any", "up", "down")
) {
  assert_diff_bundle(bundle)
  p_preference <- match.arg(p_preference)
  direction <- match.arg(direction)
  df <- bundle$results$diff_result_df
  comps <- diff_comparisons(bundle)
  out <- lapply(comps, function(cmp) {
    sig <- filter_diff_results(df[df$comparison %in% cmp, , drop = FALSE],
                               p_cutoff = p_cutoff, p_preference = p_preference,
                               effect_cutoff = effect_cutoff)
    if (direction == "up") sig <- sig[sig$direction %in% c("up", "positive"), , drop = FALSE]
    if (direction == "down") sig <- sig[sig$direction %in% c("down", "negative"), , drop = FALSE]
    unique(as.character(sig$feature_id))
  })
  stats::setNames(out, comps)
}

#' Which comparisons share their hits (an UpSet plot)
#'
#' One column per combination of comparisons, its bar the number of
#' features that are hits in exactly those comparisons and no others; the
#' dots under it say which comparisons the combination is; each row label
#' carries that comparison's total. Read left to right: the tallest columns
#' are the patterns the data actually has -- "shared by every treatment",
#' "specific to TreatB" -- which a Venn diagram of more than three sets
#' cannot show legibly.
#'
#' @inheritParams diff_hit_sets
#' @param max_combinations Largest number of combinations drawn (the most
#'   populated ones).
#' @return A `patchwork` / `ggplot` object.
#' @export
#' @family diff
plot_diff_overlap <- function(
  bundle,
  p_cutoff = 0.05,
  p_preference = c("adjusted", "raw"),
  effect_cutoff = NULL,
  direction = c("any", "up", "down"),
  max_combinations = 20L
) {
  assert_count(max_combinations, "max_combinations")
  p_preference <- match.arg(p_preference)
  direction <- match.arg(direction)
  sets <- diff_hit_sets(bundle, p_cutoff = p_cutoff, p_preference = p_preference,
                        effect_cutoff = effect_cutoff, direction = direction)
  names(sets) <- gsub("_vs_", " vs ", names(sets), fixed = TRUE)
  empty <- function(msg) {
    ggplot2::ggplot() + ggplot2::theme_void() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = msg,
                        color = "#4D4D4D", size = 4) +
      ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1)
  }
  if (length(sets) < 2L) return(empty("Overlap needs two or more comparisons."))
  all_ids <- unique(unlist(sets))
  if (!length(all_ids)) return(empty("No comparison has hits at these thresholds."))

  member <- vapply(sets, function(s) all_ids %in% s, logical(length(all_ids)))
  member <- matrix(member, nrow = length(all_ids), dimnames = list(all_ids, names(sets)))
  key <- apply(member, 1L, function(r) paste(as.integer(r), collapse = ""))
  counts <- sort(table(key), decreasing = TRUE)
  counts <- utils::head(counts, max_combinations)
  combos <- names(counts)
  set_names <- names(sets)

  bars <- data.frame(combo = factor(combos, levels = combos),
                     n = as.integer(counts), stringsAsFactors = FALSE)
  grid <- expand.grid(combo = combos, set = set_names, stringsAsFactors = FALSE)
  grid$on <- mapply(function(cb, st) substr(cb, match(st, set_names),
                                            match(st, set_names)) == "1",
                    grid$combo, grid$set)
  grid$combo <- factor(grid$combo, levels = combos)
  grid$set <- factor(grid$set, levels = rev(set_names))
  lines <- do.call(rbind, lapply(combos, function(cb) {
    on <- grid[grid$combo == cb & grid$on, , drop = FALSE]
    if (nrow(on) < 2L) return(NULL)
    data.frame(combo = factor(cb, levels = combos),
               lo = min(as.integer(on$set)), hi = max(as.integer(on$set)))
  }))

  top <- ggplot2::ggplot(bars, ggplot2::aes(x = .data$combo, y = .data$n)) +
    ggplot2::geom_col(fill = omics_colors$fg_dark %||% "#333333", width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = .data$n), vjust = -0.4, size = 3.2) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
    ggplot2::labs(
      title = "Shared hits between comparisons",
      subtitle = sprintf("%s hits, %s p < %s%s",
                         switch(direction, any = "all", up = "up-regulated",
                                down = "down-regulated"),
                         if (p_preference == "adjusted") "adjusted" else "raw",
                         format(p_cutoff),
                         if (is.null(effect_cutoff)) "" else
                           sprintf(", |effect| ≥ %s", format(effect_cutoff))),
      x = NULL, y = "features in exactly\nthis combination") +
    theme_omics_labelled() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_blank())

  dots <- ggplot2::ggplot(grid, ggplot2::aes(x = .data$combo, y = .data$set))
  if (!is.null(lines) && nrow(lines)) {
    dots <- dots + ggplot2::geom_segment(
      data = lines,
      ggplot2::aes(x = .data$combo, xend = .data$combo, y = .data$lo, yend = .data$hi),
      inherit.aes = FALSE, color = "#333333", linewidth = 0.8)
  }
  dots <- dots +
    ggplot2::geom_point(ggplot2::aes(color = .data$on), size = 3.2,
                        show.legend = FALSE) +
    ggplot2::scale_color_manual(values = c(`TRUE` = "#333333", `FALSE` = "#DADADA")) +
    ggplot2::scale_y_discrete(labels = function(x) {
      sprintf("%s (%d)", x, lengths(sets)[x])
    }) +
    ggplot2::labs(x = NULL, y = NULL) +
    theme_omics_labelled() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank())

  patchwork::wrap_plots(top, dots, ncol = 1L,
                        heights = c(2, max(1, 0.35 * length(set_names))))
}
