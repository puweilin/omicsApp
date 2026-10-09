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
      expand = ggplot2::expansion(mult = 0.22)
    ) +
    ggplot2::scale_y_discrete(labels = function(x) wrap_comparison(x)) +
    ggplot2::labs(
      title = "Hits per comparison",
      subtitle = sprintf("%s p < %s%s",
                         if (p_preference == "adjusted") "adjusted" else "raw",
                         format(p_cutoff),
                         if (is.null(effect_cutoff)) "" else
                           sprintf("\n|%s| \u2265 %s", effect_label(bundle), format(effect_cutoff))),
      # Words, not arrows: the arrows are outside the PDF device's
      # encoding, so a PDF report warned (mbcsToSbcs) and dropped them.
      x = "features, down | up", y = NULL
    ) +
    theme_omics_labelled() +
    # Up/down above the bars rather than beside them: the panel shares
    # its width with the comparison names already, and on a phone the
    # side legend left the counts cut off at the right edge.
    ggplot2::theme(legend.position = "top", legend.justification = "left")
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
#' cannot show legibly. When more combinations have hits than are drawn,
#' the subtitle says so ("top 20 of 31 combinations").
#'
#' `compact = TRUE` lays the figure out for a phone-width panel: the six
#' most populated combinations, comparisons that all share a control
#' named by their treatment alone (the control is named once, in the
#' subtitle), other names cut to two lines, and smaller counts over the
#' bars.
#'
#' @inheritParams diff_hit_sets
#' @param max_combinations Largest number of combinations drawn (the most
#'   populated ones): 20, or 6 when `compact`.
#' @param compact Lay the figure out for a narrow (phone-width) panel.
#' @return A `patchwork` / `ggplot` object.
#' @export
#' @family diff
plot_diff_overlap <- function(
  bundle,
  p_cutoff = 0.05,
  p_preference = c("adjusted", "raw"),
  effect_cutoff = NULL,
  direction = c("any", "up", "down"),
  max_combinations = if (compact) 6L else 20L,
  compact = FALSE
) {
  assert_flag(compact, "compact")
  assert_count(max_combinations, "max_combinations")
  p_preference <- match.arg(p_preference)
  direction <- match.arg(direction)
  sets <- diff_hit_sets(bundle, p_cutoff = p_cutoff, p_preference = p_preference,
                        effect_cutoff = effect_cutoff, direction = direction)
  set_labels <- overlap_set_labels(names(sets), compact)
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
  n_combos <- length(counts)
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
    ggplot2::scale_x_discrete(limits = combos) +
    ggplot2::geom_col(fill = omics_colors$fg_dark %||% "#333333", width = 0.7) +
    # Counts over the bars, shortened ("1.2k") where the columns are
    # narrow -- many of them, or a phone's panel -- and four digits ran
    # into the next column's count.
    ggplot2::geom_text(ggplot2::aes(label = overlap_count_label(
                         .data$n, short = compact || length(combos) > 12L)),
                       vjust = -0.4, size = if (compact) 2.6 else 3.2) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
    ggplot2::labs(x = NULL, y = "features in exactly\nthis combination") +
    theme_omics_labelled() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_blank())

  # The combinations on an explicit axis, in count order, for both
  # panels. Left to ggplot, the dots panel's axis took its order from the
  # first layer to name a combination -- the connecting lines, which only
  # exist for multi-comparison combinations -- so whenever any hit was
  # shared, every column of dots sat under the wrong bar.
  dots <- ggplot2::ggplot(grid, ggplot2::aes(x = .data$combo, y = .data$set)) +
    ggplot2::scale_x_discrete(limits = combos)
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
      # The count goes after the wrapping, so a long name cut short
      # never takes its total with it.
      sprintf("%s (%d)", set_labels$labels[x], lengths(sets)[x])
    }) +
    ggplot2::labs(x = NULL, y = NULL) +
    theme_omics_labelled() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank())

  # The dot rows grow with their labels: a name wrapped over three lines
  # needs three lines of height, or it runs into its neighbours.
  n_lines <- sum(lengths(strsplit(set_labels$labels, "\n", fixed = TRUE)))
  # The title and subtitle are the figure's, over both panels, rather
  # than the bar panel's: the app drops a figure's title (its card has
  # one) and keeps the panels' own.
  patchwork::wrap_plots(top, dots, ncol = 1L,
                        heights = c(2, max(1, 0.13 * length(set_names) + 0.22 * n_lines))) +
    patchwork::plot_annotation(
      title = "Shared hits between comparisons",
      subtitle = overlap_subtitle(bundle, direction, p_preference, p_cutoff, effect_cutoff,
                                  shown = length(combos), total = n_combos,
                                  control = set_labels$control),
      theme = theme_omics_labelled())
}

# What the overlap plot counts, under its title: the hits and thresholds;
# how many of the populated combinations are drawn, when not all are;
# and the shared control, when the rows are named by treatment alone.
overlap_subtitle <- function(bundle, direction, p_preference, p_cutoff, effect_cutoff,
                             shown, total, control = NULL) {
  out <- sprintf("%s hits, %s p < %s%s",
                 switch(direction, any = "all", up = "up-regulated",
                        down = "down-regulated"),
                 if (p_preference == "adjusted") "adjusted" else "raw",
                 format(p_cutoff),
                 if (is.null(effect_cutoff)) "" else
                   sprintf("\n|%s| \u2265 %s", effect_label(bundle), format(effect_cutoff)))
  if (total > shown) out <- paste0(out, sprintf("\ntop %d of %d combinations", shown, total))
  if (!is.null(control)) {
    out <- paste0(out, "\n", wrap_label(paste("each vs", control), width = 40L, max_lines = 2L))
  }
  out
}

# Row labels of the overlap plot, named by comparison as the plot shows
# it ("A vs B"), and the control they share when it is left out of them.
# Full names, wrapped, by default. Compact (a phone): when every
# comparison has the same control, the treatment alone over at most two
# lines, as the comparison plot of enrichment names its columns
# (comparison_columns()); otherwise each side over at most two lines.
# One line a side was tried: groups named "Compound alpha ..." and
# "Compound beta ..." both came out "vs Compound...", and the rows could
# not be told apart.
overlap_set_labels <- function(cmps, compact = FALSE) {
  shown <- gsub("_vs_", " vs ", cmps, fixed = TRUE)
  if (!compact) {
    return(list(labels = stats::setNames(wrap_comparison(shown), shown), control = NULL))
  }
  width <- 16L
  parts <- strsplit(cmps, "_vs_", fixed = TRUE)
  control <- comparison_columns(cmps)$control
  labels <- if (!is.null(control)) {
    cases <- vapply(parts, `[`, character(1), 1L)
    short <- wrap_label(cases, width = width, max_lines = 2L)
    # Two treatments that only differ past the cut keep their full names.
    if (anyDuplicated(short)) wrap_label(cases, width = width, max_lines = 4L) else short
  } else {
    vapply(parts, function(p) {
      if (length(p) != 2L) return(wrap_label(paste(p, collapse = " vs "), width = width, max_lines = 2L))
      paste0(wrap_label(p[1L], width = width, max_lines = 2L), "\n",
             wrap_label(paste("vs", p[2L]), width = width, max_lines = 2L))
    }, character(1))
  }
  list(labels = stats::setNames(labels, shown), control = control)
}

# A bar's count; past 999 as "1.2k" when `short`.
overlap_count_label <- function(n, short = FALSE) {
  out <- as.character(n)
  if (!short) return(out)
  big <- n >= 1000
  out[big] <- paste0(trimws(formatC(n[big] / 1000, format = "fg", digits = 2)), "k")
  out
}

#' @rdname wrap_label
#' @details `wrap_comparison()` wraps a `"<case> vs <control>"` label with
#'   each side on its own lines -- the case over at most two, then
#'   `"vs <control>"` on one -- so a long case name cut short never takes
#'   the control, half of what the label says, with it.
#' @export
wrap_comparison <- function(x, width = 22L) {
  assert_count(width, "width")
  x <- gsub("_vs_", " vs ", as.character(x), fixed = TRUE)
  vapply(x, function(s) {
    if (is.na(s) || nchar(s, type = "width") <= width) return(s)
    sides <- strsplit(s, " vs ", fixed = TRUE)[[1L]]
    if (length(sides) != 2L) return(wrap_label(s, width = width, max_lines = 3L))
    paste0(wrap_label(sides[1L], width = width, max_lines = 2L), "\n",
           wrap_label(paste("vs", sides[2L]), width = width, max_lines = 1L))
  }, character(1), USE.NAMES = FALSE)
}

#' @rdname wrap_label
#' @details `comparison_label_lines()` gives the number of lines
#'   `wrap_comparison()` uses for each label, for sizing a plot to fit;
#'   with `compact = TRUE`, the lines of the row labels
#'   [plot_diff_overlap()] draws with `compact = TRUE` for the comparisons
#'   `x` (`"<case>_vs_<control>"` names, as [diff_comparisons()] gives
#'   them).
#' @param compact For `comparison_label_lines()`: count the lines of
#'   `plot_diff_overlap(compact = TRUE)`'s row labels instead.
#' @export
comparison_label_lines <- function(x, width = 22L, compact = FALSE) {
  assert_flag(compact, "compact")
  labels <- if (compact) overlap_set_labels(as.character(x), compact = TRUE)$labels
            else wrap_comparison(x, width = width)
  unname(lengths(strsplit(labels, "\n", fixed = TRUE)))
}
