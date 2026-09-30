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
  if (length(cases) == length(available)) {
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
