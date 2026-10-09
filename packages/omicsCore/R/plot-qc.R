#' QC visualizations
#'
#' Builds standard ggplot panels from a [run_qc()] bundle. They are drawn
#' from what [run_qc()] stored in it, so neither the input nor the cleaned
#' input is needed.
#'
#' @param bundle An `analysis_bundle` produced by [run_qc()].
#' @param view One of `"missing"`, `"pca"`, `"connectivity"`, `"imputation"`.
#' @param color_by Optional name of a column in the cleaned input's `meta_df`
#'   used to color samples in the `"pca"` view.
#' @param reference Optional group of `color_by` coloured first, as every
#'   figure colours the study design's reference group (see
#'   [group_palette()]); usually `study_design(input)$reference`. Ignored
#'   when it is not a group of `color_by`.
#' @param ... Reserved for future arguments.
#'
#' @return A `ggplot` object.
#' @export
#' @family qc
plot_qc <- function(bundle,
                    view = c("missing", "depth", "pca", "connectivity",
                             "imputation"),
                    color_by = NULL,
                    reference = NULL,
                    ...) {
  assert_string(color_by, "color_by", allow_null = TRUE)
  assert_label(reference, "reference", allow_null = TRUE)
  if (!is_analysis_bundle(bundle) || !identical(bundle$analysis_name, "run_qc")) {
    stop("`bundle` must be an analysis_bundle from run_qc().")
  }
  view <- match.arg(view)

  switch(view,
    missing      = plot_qc_missing(bundle),
    depth        = plot_qc_depth(bundle),
    pca          = plot_qc_pca(bundle, color_by = color_by, reference = reference),
    connectivity = plot_qc_connectivity(bundle),
    imputation   = plot_qc_imputation(bundle)
  )
}

# ---- views -------------------------------------------------------------

MISSING_FILL <- "#2C3E99"

# Samples and features above a cutoff -- the ones QC flags or removes.
# Amber rather than the up/down red or blue: being over a missingness
# cutoff is a warning, not a direction of change, and red here would
# read as "up". The same amber as omics_colors$conc_up_down, written out
# because plot-tokens.R is sourced after this file.
MISSING_OVER_FILL <- "#E0A030"

# The dashed cutoff line: dark enough to read over both fills.
MISSING_CUTOFF_COLOUR <- "#374151"

# ---- per-sample panels, shared by the missingness and depth views ----
#
# How many samples a per-sample panel draws as named bars. The app gives
# two stacked panels 360 px between them, and after titles and axes each
# has about 90 px for its bars: ten names at 7 pt is what fits there
# without the labels running into one another (24 overlapped at the old
# cap of 30, and 12 still touched). Past it, every sample is a point on
# a ranked curve and the ones that matter are named in the subtitle.
SAMPLE_MAX_NAMED_BARS <- 10L

# From this many bars the names drop from the theme's size to 7 pt, so
# that up to SAMPLE_MAX_NAMED_BARS fit.
SAMPLE_SMALL_LABEL_FROM <- 7L

# The names these had when only the missingness panel used them.
MISSING_MAX_NAMED_SAMPLES <- SAMPLE_MAX_NAMED_BARS
MISSING_SMALL_LABEL_FROM <- SAMPLE_SMALL_LABEL_FROM

# Sample names on a bar chart's y axis: smaller once there are enough
# bars to crowd, and cut at 24 characters.
sample_label_theme <- function(n) {
  if (n >= SAMPLE_SMALL_LABEL_FROM) {
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 7))
  }
}

# The y axis of a ranked curve: rank 1, the worst, at the top.
sample_rank_scale <- function() {
  ggplot2::scale_y_reverse(breaks = function(lim) {
    b <- scales::breaks_pretty(4)(lim)
    b[b >= 1]  # there is no rank 0
  })
}

# The feature histogram's bins. A missing rate can only be k / n -- with
# 24 samples there are 25 possible values -- so a bin per possible value
# is exact, and bins that ignore that split some values across two bars
# and leave others empty, which draws a comb that is not in the data.
# Past this many possible values, neighbouring values share a bin.
MISSING_MAX_BINS <- 40L

# Two panels, because samples and features are different questions:
#
# * per sample -- which sample is bad, and is any over the sample cutoff
#   (if one was set);
# * per feature -- how many features the missing-value filter removes,
#   which is the decision this panel exists to inform. The cutoff is
#   drawn on it and the bars past it are coloured, so the count in the
#   subtitle can be seen as well as read.
#
# Each panel has its own x axis. They used to share one, and since a
# feature missing in every sample is common, that shared axis ran to
# 100% and squeezed sample bars of 5-20% against its left edge.
plot_qc_missing <- function(bundle) {
  miss <- bundle$results$qc_summary$missingness
  rules <- missing_rules(bundle)
  n_samples <- nrow(miss$sample_metrics)
  patchwork::wrap_plots(
    missing_unaligned(
      plot_missing_by_sample(miss$sample_metrics, cutoff = rules$sample_cutoff)),
    missing_unaligned(
      plot_missing_by_feature(miss$feature_metrics,
                              cutoff = rules$feature_cutoff,
                              filter = rules$filter,
                              n_removed = rules$n_removed,
                              n_samples = n_samples)),
    ncol = 1,
    # Named bars need every pixel of their half; the ranked curve does
    # not, and the histogram's bars past the cutoff are short.
    heights = if (n_samples > SAMPLE_MAX_NAMED_BARS) c(1, 1.25) else c(1, 1)
  )
}

# The two panels' plotting areas are not lined up. Their x axes measure
# different things now, and lining the histogram up under the sample
# names gave a third of a phone-width panel to empty space. free() is
# patchwork >= 1.2; an older one lines them up, which is only wasteful.
missing_unaligned <- function(p) {
  if (utils::packageVersion("patchwork") >= "1.2.0") patchwork::free(p) else p
}

# The cutoffs and the outcome QC recorded, read from wherever the bundle
# keeps them. The missingness settings are what qc_missingness() actually
# used, so they come first; the run's params are the same values one
# level up, and an older saved bundle may carry only one of the two, or
# neither -- then the panel draws no line and states no count, rather
# than guess one.
missing_rules <- function(bundle) {
  miss <- bundle$results$qc_summary$missingness
  settings <- miss$settings
  params <- bundle$params
  cutoff_or_null <- function(x) {
    if (is.numeric(x) && length(x) == 1L && !is.na(x) && x >= 0 && x <= 1) x
  }
  filter <- settings$missing_filter %||% params$missing_filter %||% "global"
  if (!filter %in% c("global", "any_group", "all_groups")) filter <- "global"

  # The number removed is the number QC flagged -- every flagged feature
  # is dropped from the cleaned layer, and nothing else drops a feature,
  # so it equals n_features_in - n_features_out. The input counts are the
  # fallback for a bundle saved without the flagged list.
  n_removed <- if (!is.null(miss$flagged_features)) {
    length(unique(miss$flagged_features))
  } else {
    info <- bundle$input_info
    if (is.numeric(info$n_features_in) && is.numeric(info$n_features_out)) {
      as.integer(info$n_features_in - info$n_features_out)
    }
  }

  list(
    feature_cutoff = cutoff_or_null(settings$feature_missing_cutoff %||%
                                      params$missing_threshold),
    sample_cutoff = cutoff_or_null(settings$sample_missing_cutoff %||%
                                     params$sample_missing_threshold),
    filter = filter,
    n_removed = n_removed
  )
}

# Anchored at 0, so the reader can see where the floor is, and never
# past 1. The floor of 5% stops an all-but-complete dataset from being
# magnified into what looks like a problem. Used for the sample axis
# only: the feature panel always shows the whole 0-100%, because
# features missing in most samples are exactly what it is about.
missing_axis_upper <- function(rates) {
  rates <- rates[!is.na(rates)]
  if (!length(rates)) return(1)
  max(0.05, min(1, max(rates) * 1.15 + 0.01))
}

# Whether the sample cutoff is close enough to the samples to draw. A
# cutoff far past every sample (80% when the worst is at 15%) would
# stretch the axis back out and squeeze the bars the way the shared axis
# did; the subtitle then says none is over it, which is all the line
# would have shown.
missing_cutoff_in_view <- function(rates, cutoff) {
  if (is.null(cutoff)) return(FALSE)
  rates <- rates[!is.na(rates)]
  top <- if (length(rates)) max(rates) else 0
  cutoff <= max(0.05, 2 * top)
}

# How many of the worst samples to name when there are too many to name
# them all, at most.
SAMPLE_LABEL_WORST <- 3L

# How many characters of sample names the subtitle's second line holds.
# With its "Shallow: " in front and ", ..." after, 28 is what fits a
# phone-width panel (293 px at 9 pt). Short names fit three; long ones
# fewer, each whole where it can be -- a name cut short loses the part
# that tells samples apart ("Patient_005_pla...").
SAMPLE_NAME_CHARS <- 28L

# "S01, S02, S03". When the list is a set (the samples over a cutoff)
# and not all of it fits, it ends in an ellipsis: the line above it
# already says how many there are.
sample_name_list <- function(ids, n = SAMPLE_LABEL_WORST, more = FALSE) {
  ids <- truncate_pathway_name(as.character(ids), SAMPLE_NAME_CHARS)
  fits <- cumsum(nchar(ids) + 2L) - 2L <= SAMPLE_NAME_CHARS
  k <- max(1L, min(n, sum(fits)))
  shown <- paste(ids[seq_len(min(k, length(ids)))], collapse = ", ")
  if (more && length(ids) > k) paste0(shown, ", \u2026") else shown
}
missing_name_list <- sample_name_list

# "3 over the 25% cutoff", after a middle dot; nothing when no cutoff
# was set.
missing_sample_cutoff_text <- function(n_over, cutoff) {
  if (is.null(cutoff)) return("")
  pct <- format_missing_pct(cutoff)
  if (n_over == 0L) sprintf(" \u00b7 none over the %s cutoff", pct)
  else sprintf(" \u00b7 %d over the %s cutoff", n_over, pct)
}

# "8 samples", and says so when none of them misses a value: empty
# bars on their own look like a plot that failed to draw.
missing_sample_count <- function(df) {
  out <- sprintf("%d samples", nrow(df))
  if (nrow(df) && all(df$missing_rate == 0, na.rm = TRUE)) {
    out <- paste(out, "\u00b7 no missing values")
  }
  out
}

format_missing_pct <- function(x) scales::label_percent(accuracy = 1)(x)

# A bar per sample, worst first, while the names fit; past that a
# sorted curve: every sample is still a point, rank on y, and the shape
# of the curve is the distribution -- a flat line with a short tail
# reads very differently from a steady slope. The worst few (or the
# ones over the cutoff) keep their names, in the subtitle.
#
# Samples over the sample cutoff are amber, with the cutoff dashed.
plot_missing_by_sample <- function(sample_df, cutoff = NULL) {
  df <- sample_df[order(sample_df$missing_rate, decreasing = TRUE), ,
                  drop = FALSE]
  df$over <- if (is.null(cutoff)) {
    rep(FALSE, nrow(df))
  } else {
    !is.na(df$missing_rate) & df$missing_rate > cutoff
  }
  show_line <- missing_cutoff_in_view(df$missing_rate, cutoff)
  upper <- missing_axis_upper(c(df$missing_rate, if (show_line) cutoff))

  if (nrow(df) > SAMPLE_MAX_NAMED_BARS) {
    return(plot_missing_sample_curve(df, upper, cutoff, show_line))
  }
  # Reversed, because a discrete y axis is drawn bottom-up and the
  # worst sample belongs at the top.
  df$sample_id <- factor(df$sample_id, levels = rev(df$sample_id))

  ggplot2::ggplot(df, ggplot2::aes(x = .data$missing_rate,
                                   y = .data$sample_id,
                                   fill = .data$over)) +
    ggplot2::geom_col(width = 0.7) +
    missing_over_scale("fill") +
    missing_cutoff_line(cutoff, show_line) +
    ggplot2::scale_y_discrete(labels = function(x) truncate_pathway_name(x, 24L)) +
    ggplot2::scale_x_continuous(
      labels = scales::label_percent(),
      limits = c(0, upper),
      expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::labs(title = "Missing rate per sample",
                  subtitle = paste0(missing_sample_count(df),
                                    missing_sample_cutoff_text(sum(df$over), cutoff)),
                  x = NULL, y = NULL) +
    theme_omicsCore() +
    sample_label_theme(nrow(df))
}

# `df` is already sorted worst-first, with `over` set.
plot_missing_sample_curve <- function(df, upper = 1, cutoff = NULL,
                                      show_line = FALSE) {
  df$rank <- seq_len(nrow(df))
  # Named in the subtitle rather than beside their points: the worst few
  # sit at almost the same rank, so on the plot their labels land on top
  # of one another. With a cutoff, the samples over it are the ones to
  # name; otherwise the worst. On a second line, which a phone-width
  # panel needs for long names.
  # A complete matrix -- every counts layer -- has no worst sample to
  # name: "Worst: R01, R02, R03" over a column of zeros said there was.
  names <- if (any(df$over)) {
    paste0("\nOver it: ", sample_name_list(df$sample_id[df$over], more = TRUE))
  } else if (!all(df$missing_rate == 0, na.rm = TRUE)) {
    paste0("\nWorst: ", sample_name_list(df$sample_id))
  }
  # "ranked" is left to the axis title: the first line has to fit a
  # phone-width panel, about 42 characters in the app.
  subtitle <- paste0(missing_sample_count(df),
                     missing_sample_cutoff_text(sum(df$over), cutoff),
                     names)

  ggplot2::ggplot(df, ggplot2::aes(x = .data$missing_rate, y = .data$rank,
                                   colour = .data$over)) +
    ggplot2::geom_point(size = 1.3, alpha = 0.8) +
    missing_over_scale("colour") +
    missing_cutoff_line(cutoff, show_line) +
    sample_rank_scale() +
    ggplot2::scale_x_continuous(labels = scales::label_percent(),
                                limits = c(0, upper),
                                expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::labs(title = "Missing rate per sample", subtitle = subtitle,
                  x = NULL, y = "Rank") +
    theme_omicsCore()
}

missing_over_scale <- function(aesthetic) {
  ggplot2::scale_discrete_manual(
    aesthetics = aesthetic,
    values = c(`FALSE` = MISSING_FILL, `TRUE` = MISSING_OVER_FILL),
    guide = "none")
}

# A white line under the dashes: drawn on its own, the dark-grey dash
# all but vanished where it crossed the dark-blue bars, which is where
# the reader needs it -- on the sample nearest the cutoff.
missing_cutoff_line <- function(cutoff, show = !is.null(cutoff)) {
  if (!show || is.null(cutoff)) return(NULL)
  list(
    ggplot2::geom_vline(xintercept = cutoff, colour = "white", linewidth = 1.6),
    ggplot2::geom_vline(xintercept = cutoff, linetype = "dashed",
                        colour = MISSING_CUTOFF_COLOUR, linewidth = 0.5))
}

# What the feature filter compared with the cutoff, in words, for the
# axis -- and the rule for removal, for the subtitle. Worded like the
# app's own summary card ("filtered at 50% missing in every group").
missing_filter_words <- function(filter) {
  switch(filter,
    any_group = list(axis = "Lowest missing rate among groups",
                     rule = "more than %s missing in every group"),
    all_groups = list(axis = "Highest missing rate among groups",
                      rule = "more than %s missing in any group"),
    # No axis title for the overall rate: the panel's title already says
    # what x is, and the app's 360 px leave no height to spare.
    list(axis = NULL, rule = "missing in more than %s of samples"))
}

# The bins of the feature histogram, as edges on [0, 1]. Missing rates
# are fractions k / d, so the edges sit half-way between possible
# values: each bin holds whole values, never part of one. 0% -- complete
# features, usually the largest group -- gets a bin of its own; past
# MISSING_MAX_BINS possible values, neighbours are paired (or tripled,
# ...) after it.
missing_bins <- function(rate, n_samples = NA_integer_) {
  d <- missing_rate_denominator(rate, n_samples)
  if (is.na(d)) {
    edges <- seq(0, 1, length.out = MISSING_MAX_BINS + 1L)
  } else {
    g <- ceiling(d / MISSING_MAX_BINS)
    # In whole values first, then divided: every edge is half-way
    # between two possible values, so none can fall on one.
    edges <- c(-0.5, 0.5 + g * (0:ceiling(d / g))) / d
  }
  structure(edges, denominator = d)
}

# The smallest d for which every rate is a whole number of d-ths. For
# the overall rate that is the number of samples (or a divisor of it);
# for the group rates it is the least common multiple of the group
# sizes, which the bundle does not record but the rates themselves give
# away. NA when no small d fits -- then the bins are simply even.
missing_rate_denominator <- function(rate, n_samples = NA_integer_,
                                     max_d = 1000L) {
  r <- unique(round(rate[!is.na(rate)], 9))
  if (!length(r)) return(1L)
  cands <- seq_len(max_d)
  if (!is.na(n_samples) && n_samples >= 1L && n_samples <= max_d) {
    # Try the sample count first: for the overall rate it is the answer,
    # and a smaller d that happens to fit too (all rates multiples of a
    # quarter, with 24 samples) would draw coarser bins than the data has.
    cands <- c(n_samples, setdiff(cands, n_samples))
  }
  for (d in cands) {
    if (all(abs(r * d - round(r * d)) < 1e-6)) return(as.integer(d))
  }
  NA_integer_
}

# Counts of features per bin, split into kept and removed.
missing_feature_bins <- function(rate, cutoff, n_samples = NA_integer_) {
  edges <- missing_bins(rate, n_samples)
  bin <- findInterval(rate, edges, rightmost.closed = TRUE, all.inside = TRUE)
  removed <- if (is.null(cutoff)) rep(FALSE, length(rate)) else rate > cutoff
  counts <- as.data.frame(table(bin = factor(bin, levels = seq_len(length(edges) - 1L)),
                                removed = factor(removed, levels = c(FALSE, TRUE))),
                          responseName = "n", stringsAsFactors = FALSE)
  counts$bin <- as.integer(counts$bin)
  counts <- counts[counts$n > 0L, , drop = FALSE]
  # Each bar stands at the middle of the rates its bin can hold -- 0% at
  # 0%, 100% at 100% -- not at the middle of an interval that runs half
  # a value past either end.
  d <- attr(edges, "denominator")
  lo <- edges[counts$bin]
  hi <- edges[counts$bin + 1L]
  counts$x <- if (is.na(d)) {
    (lo + hi) / 2
  } else {
    (pmax(ceiling(lo * d), 0) + pmin(floor(hi * d), d)) / 2 / d
  }
  counts$width <- 0.85 * stats::median(diff(edges))
  counts$removed <- counts$removed == "TRUE"
  counts[order(counts$x, counts$removed), , drop = FALSE]
}

# A histogram of feature missing rates over the whole 0-100%, with the
# filter's cutoff dashed and the features it removes in amber. y is on a
# square-root scale: complete features usually outnumber those near the
# cutoff fifty to one, and on a linear axis the bars the filter removes
# -- the ones this panel is about -- are a pixel high.
#
# For the group rules the cutoff is not compared with the overall rate,
# so drawing that would put features on the wrong side of the line. The
# panel then plots the rate the filter used (the lowest or highest group
# rate) and says so on its axis.
plot_missing_by_feature <- function(feature_df, cutoff = NULL,
                                    filter = "global", n_removed = NULL,
                                    n_samples = NA_integer_) {
  group_rule <- !identical(filter, "global")
  has_group_rate <- group_rule && !is.null(feature_df$filter_missing_rate)
  rate <- if (has_group_rate) feature_df$filter_missing_rate
          else feature_df$missing_rate
  rate <- rate[!is.na(rate)]
  # A group rule with no group rates recorded: the overall rate is all
  # there is to draw, and a cutoff line over it would not be the line
  # the filter used.
  draw_cutoff <- !is.null(cutoff) && (!group_rule || has_group_rate)
  words <- missing_filter_words(if (has_group_rate) filter else "global")
  # The group rates are fractions of a group, not of all samples.
  if (has_group_rate) n_samples <- NA_integer_

  subtitle <- missing_feature_subtitle(length(rate), rate, cutoff, filter,
                                       n_removed)
  labs <- ggplot2::labs(title = "Missing rate per feature",
                        subtitle = subtitle,
                        x = words$axis, y = "Features (\u221a scale)")
  x_scale <- ggplot2::scale_x_continuous(
    labels = scales::label_percent(), breaks = seq(0, 1, 0.25))

  if (!length(rate)) {
    return(ggplot2::ggplot() + x_scale +
             ggplot2::coord_cartesian(xlim = c(0, 1)) + labs +
             theme_omicsCore())
  }

  bins <- missing_feature_bins(rate, if (draw_cutoff) cutoff, n_samples)
  half <- bins$width[1L] / 2
  xlim <- c(min(0, bins$x - half), max(1, bins$x + half))
  top <- max(tapply(bins$n, bins$x, sum))

  ggplot2::ggplot(bins, ggplot2::aes(x = .data$x, y = .data$n,
                                     fill = .data$removed)) +
    ggplot2::geom_col(width = bins$width[1L],
                      position = ggplot2::position_stack(reverse = TRUE)) +
    missing_over_scale("fill") +
    missing_cutoff_line(cutoff, draw_cutoff) +
    (if (draw_cutoff) missing_cutoff_label(cutoff, top)) +
    x_scale +
    ggplot2::scale_y_sqrt(breaks = missing_count_breaks(top),
                          labels = scales::label_comma()) +
    # Room above the tallest bar for the cutoff's label.
    ggplot2::coord_cartesian(xlim = xlim, ylim = c(0, top * 1.4),
                             expand = FALSE) +
    labs +
    theme_omicsCore() +
    # Smaller than the theme's, to fit the height of a panel the app
    # gives 200 px.
    ggplot2::theme(axis.title.y = ggplot2::element_text(size = ggplot2::rel(0.8)))
}

# The cutoff's value, written at the top of its line: on the side with
# more room, so a cutoff near 100% does not run off the panel.
missing_cutoff_label <- function(cutoff, top) {
  right <- cutoff <= 0.75
  ggplot2::annotate("text", x = cutoff, y = top * 1.35,
                    label = paste0(if (right) " " else "",
                                   format_missing_pct(cutoff), " cutoff",
                                   if (right) "" else " "),
                    hjust = if (right) 0 else 1, vjust = 1, size = 3,
                    colour = MISSING_CUTOFF_COLOUR)
}

# Round counts that spread out on a square-root axis: 0, 10, 100, 500,
# 1,000 ... rather than evenly spaced counts that bunch at the top.
missing_count_breaks <- function(top) {
  if (!is.finite(top) || top <= 0) return(0)
  steps <- c(1, 2, 5) * rep(10^(0:7), each = 3)
  cand <- c(0, steps[steps <= top])
  # Three or four labelled counts: 0, the largest round count under the
  # top, and one or two in between on the square-root scale.
  if (length(cand) <= 4L) return(cand)
  hi <- cand[length(cand)]
  mids <- vapply(c(1 / 9, 4 / 9), function(f) {
    cand[which.min(abs(sqrt(cand) - sqrt(hi * f)))]
  }, numeric(1))
  unique(c(0, mids, hi))
}

# "5,000 features, 312 removed" over "(missing in more than 50% of
# samples)" -- on two lines, because one does not fit a phone-width
# panel, and the second line is the rule, not the result.
missing_feature_subtitle <- function(n, rate, cutoff, filter, n_removed) {
  head <- sprintf("%s features", format(n, big.mark = ","))
  if (n > 0L && all(rate == 0)) {
    head <- paste(head, "\u00b7 no missing values")
  }
  if (is.null(cutoff) || is.null(n_removed)) return(head)
  result <- if (n_removed == 0L) "none removed"
            else sprintf("%s removed", format(n_removed, big.mark = ","))
  sprintf("%s \u00b7 %s\n(%s)", head, result,
          sprintf(missing_filter_words(filter)$rule, format_missing_pct(cutoff)))
}

plot_qc_pca <- function(bundle, color_by = NULL, reference = NULL) {
  pd <- qc_pca_data(bundle)
  if (!is.null(pd$error)) stop(pd$error, call. = FALSE)
  scores <- as.data.frame(pd$scores)
  scores$sample_id <- rownames(scores)
  rownames(scores) <- NULL
  meta <- pd$meta_df

  if (!is.null(color_by)) {
    if (!color_by %in% colnames(meta)) {
      stop("`color_by` not found in the cleaned input's meta_df: ", color_by)
    }
    scores[[color_by]] <- meta[scores$sample_id, color_by]
  }

  var_pct <- pd$var_pct

  group_vals <- if (is.null(color_by)) NULL else scores[[color_by]]
  mapping <- if (is.null(color_by)) {
    ggplot2::aes(x = .data$PC1, y = .data$PC2)
  } else if (use_group_shape(group_vals)) {
    ggplot2::aes(x = .data$PC1, y = .data$PC2,
                 color = .data[[color_by]], shape = .data[[color_by]])
  } else {
    ggplot2::aes(x = .data$PC1, y = .data$PC2,
                 color = .data[[color_by]])
  }

  # Said on the plot rather than left implicit. A counts matrix carries
  # every annotated gene, so a fifth of them being zero in every sample
  # is ordinary -- but "PCA of 63,241 features" and "PCA of the 49,000
  # that vary" are different claims, and only one of them is true.
  dropped <- pd$n_dropped %||% 0L
  subtitle <- if (dropped > 0L) {
    sprintf("%s features; %s constant across all samples, excluded",
            format(pd$n_features, big.mark = ","),
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
    group_legend_scales(group_vals, reference = reference) +
    theme_omicsCore()
}

plot_qc_connectivity <- function(bundle) {
  stats_df <- qc_connectivity_from_outliers(bundle$results$qc_summary$outliers) %||%
    bundle$results$plot_data$connectivity
  if (is.null(stats_df)) {
    # outlier_method was something else, in a result saved before
    # run_qc() worked this out itself: recompute it on the cleaned input
    # so users always see this view.
    cleaned <- bundle$results$cleaned_input
    stats_df <- qc_outliers_connectivity(qc_log_scale(cleaned)$mat,
                                         sd_threshold = bundle$params$outlier_sd_threshold)$stats
  }
  if (!is.data.frame(stats_df)) stop(stats_df$error, call. = FALSE)
  stats_df <- stats_df[order(stats_df$mean_correlation), , drop = FALSE]
  stats_df$sample_id <- factor(stats_df$sample_id, levels = stats_df$sample_id)

  ggplot2::ggplot(stats_df,
                  ggplot2::aes(x = .data$sample_id,
                               y = .data$mean_correlation,
                               fill = .data$is_outlier)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = c(`TRUE` = "#C0392B", `FALSE` = "#2C3E99")) +
    ggplot2::scale_x_discrete(labels = function(x) truncate_pathway_name(x, 20L)) +
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
  imp <- bundle$results$qc_summary$imputation
  cleaned <- bundle$results$cleaned_input
  pd <- bundle$results$plot_data$imputation
  if (is.null(imp) && is.null(cleaned$raw_mat)) {
    stop("This bundle has no imputation step (run_qc was called with impute_method='none').")
  }
  if (!is.null(pd)) {
    # Observed against imputed values, both from the matrix that goes
    # downstream, on the log scale. The observed side is a summary of
    # every observed value (see qc_plot_data()), smoothed as the full set
    # would have been.
    imputed <- bundle$results$cleaning$imputed_values
    at <- bundle$results$cleaning$assay_type
    if (!isTRUE(at %in% LOG_SCALE_ASSAY_TYPES)) imputed <- log2(imputed)
    curves <- rbind(
      density_curve(pd$observed_quantiles, "observed", bw = pd$observed_bw),
      density_curve(imputed, "imputed"))
  } else {
    if (!is.null(imp)) {
      # `raw_mat` was the wrong "before": after normalize_omics() it holds
      # the linear values, so the two curves were on different scales.
      mat <- cleaned$expr_mat
      if (!cleaned$assay_type %in% LOG_SCALE_ASSAY_TYPES) mat <- log2(mat)
      before <- as.numeric(mat[-imp$imputed_cells])
      after <- as.numeric(mat[imp$imputed_cells])
      first <- "observed"
    } else {
      # Bundles from before the imputation record.
      before <- as.numeric(cleaned$raw_mat)
      after <- as.numeric(cleaned$expr_mat)
      first <- "raw"
    }
    curves <- rbind(density_curve(before, first), density_curve(after, "imputed"))
  }

  ggplot2::ggplot(curves,
                  ggplot2::aes(x = .data$value, y = .data$density,
                               fill = .data$type, color = .data$type)) +
    ggplot2::geom_area(alpha = 0.35, position = "identity") +
    ggplot2::scale_fill_manual(values = c(raw = "#9AA3AE", observed = "#9AA3AE",
                                          imputed = "#1FBF9E")) +
    ggplot2::scale_color_manual(values = c(raw = "#9AA3AE", observed = "#9AA3AE",
                                           imputed = "#1FBF9E")) +
    ggplot2::labs(
      title = "Imputation effect on intensity distribution",
      x = "Value",
      y = "Density",
      fill = NULL, color = NULL
    ) +
    theme_omicsCore()
}

# A kernel density as a data frame, the way geom_density() would have drawn
# it. Computed here so a curve can be drawn from a summary of the values,
# with the bandwidth the full set had. Fewer than two finite values have no
# density, and draw nothing.
density_curve <- function(x, type, bw = NULL) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) {
    return(data.frame(value = numeric(0), density = numeric(0),
                      type = character(0), stringsAsFactors = FALSE))
  }
  if (is.null(bw) || !is.finite(bw) || bw <= 0) bw <- "nrd0"
  d <- tryCatch(stats::density(x, bw = bw), error = function(e) NULL)
  if (is.null(d)) {
    return(data.frame(value = numeric(0), density = numeric(0),
                      type = character(0), stringsAsFactors = FALSE))
  }
  data.frame(value = d$x, density = d$y, type = type, stringsAsFactors = FALSE)
}

# What the PCA panel draws: stored by run_qc(), or -- for results saved
# before it stored it -- worked out from the cleaned input they carry.
qc_pca_data <- function(bundle) {
  pd <- bundle$results$plot_data
  if (!is.null(pd$pca)) return(c(pd$pca, list(meta_df = pd$meta_df)))
  cleaned <- bundle$results$cleaned_input
  if (is.null(cleaned)) {
    stop("This QC result carries nothing to draw a PCA from; run run_qc() again.",
         call. = FALSE)
  }
  # On the scale qc_outliers() used, so the plot shows what was tested:
  # a PCA of raw counts is a plot of library size and a few huge genes.
  mat <- mean_impute_rows(qc_log_scale(cleaned)$mat)
  if (ncol(mat) < 2L) {
    return(list(error = "Need at least 2 samples to draw a PCA scatter."))
  }
  pca <- pca_over_samples(mat)
  list(scores = pca$x[, 1:2, drop = FALSE],
       var_pct = (pca$sdev^2) / sum(pca$sdev^2) * 100,
       n_features = nrow(mat),
       n_dropped = attr(pca, "n_dropped") %||% 0L,
       meta_df = cleaned$meta_df)
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
  assert_number(base_size, "base_size", lower = 0)
  assert_string(base_family, "base_family", allow_empty = TRUE)
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.05)),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "#E5E7EB"),
      axis.line = ggplot2::element_line(color = "#9AA3AE"),
      strip.text = ggplot2::element_text(face = "bold"),
      # Titles from the plot's left edge, not the panel's: with long
      # axis labels the panel starts half-way across and the title ran
      # off the right side ("Hits per compar").
      plot.title.position = "plot"
    ) +
    legend_key_spacing()
}

# Room between legend entries, so a label wrapped over three lines does
# not run into the next one. The element is ggplot2 >= 3.5; older
# versions do not know it and would refuse the theme.
legend_key_spacing <- function() {
  if (!"legend.key.spacing.y" %in% names(ggplot2::get_element_tree())) {
    return(ggplot2::theme())
  }
  ggplot2::theme(legend.key.spacing.y = ggplot2::unit(5, "pt"))
}

# ---- depth (RNA-seq) ---------------------------------------------------
# What replaces the missingness panels for a counts matrix, where
# missingness is 0% for every feature and the panel says nothing.
#
# Two questions, drawn together because the pair is what distinguishes
# the two problems: a shallow library drops both total counts and genes
# detected, while a degraded sample drops detection with the total
# holding up.
#
# Drawn by the same rules as the missingness panel's samples: named bars
# while the names fit, worst at the top, a ranked curve past
# SAMPLE_MAX_NAMED_BARS; flagged samples amber, their cutoff dashed.

# A library below this fraction of the median is flagged as shallow --
# the cutoff qc_depth_outliers() uses by default, and so the one the
# app's caption under this panel reports. run_qc() does not record a
# depth cutoff of its own.
DEPTH_LOW_RATIO <- 0.3

plot_qc_depth <- function(bundle) {
  depth <- bundle$results$qc_summary$depth
  if (is.null(depth) || nrow(depth) == 0L) {
    stop("This QC bundle carries no depth summary.", call. = FALSE)
  }
  # Flagged once, on library size, and shown amber in both panels: that
  # a shallow library also detects fewer genes is the expected half of
  # the pair, and a low detection bar that is *not* amber is the
  # degraded sample worth a second look.
  depth$.low <- depth$sample_id %in% qc_depth_outliers(depth, DEPTH_LOW_RATIO)
  words <- depth_words(bundle$input_info$omics_type)
  patchwork::wrap_plots(
    missing_unaligned(plot_depth_library(depth, words)),
    missing_unaligned(plot_depth_detection(depth, words)),
    ncol = 1
  )
}

# What the two numbers are called. A proteomics run has no library and
# is not sequenced: its total is the summed intensity, and a sample far
# below the rest was under-loaded or poorly ionised, not "shallow".
DEPTH_WORDS <- list(
  library = list(total = "Library size per sample",
                 detected = "Features detected per sample",
                 low = "shallow library", low_prefix = "Shallow:"),
  intensity = list(total = "Total intensity per sample",
                   detected = "Features quantified per sample",
                   low = "low total intensity", low_prefix = "Low:")
)

depth_words <- function(omics_type) {
  if (identical(omics_type, "proteomics")) DEPTH_WORDS$intensity else DEPTH_WORDS$library
}

# Ordered worst-first for the same reason the missingness panel is: the
# question is which sample is bad, and sorting is what answers it
# without reading every label. Lowest first; the bar chart reverses its
# levels so the lowest is drawn at the top.
depth_ordered <- function(depth, col) {
  depth <- depth[order(depth[[col]]), , drop = FALSE]
  if (is.null(depth$.low)) depth$.low <- rep(FALSE, nrow(depth))
  depth
}

plot_depth_library <- function(depth, words = DEPTH_WORDS$library) {
  if (is.null(depth$.low)) {
    depth$.low <- depth$sample_id %in% qc_depth_outliers(depth, DEPTH_LOW_RATIO)
  }
  d <- depth_ordered(depth, "library_size")
  med <- stats::median(d$library_size[is.finite(d$library_size)])
  cutoff <- if (is.finite(med) && med > 0) DEPTH_LOW_RATIO * med
  n_low <- sum(d$.low)
  # Short enough for a phone-width panel (about 42 characters there).
  head <- sprintf("%d samples \u00b7 %s below %s of median",
                  nrow(d), if (n_low == 0L) "none" else as.character(n_low),
                  format_missing_pct(DEPTH_LOW_RATIO))
  if (is.null(cutoff)) head <- sprintf("%d samples", nrow(d))
  names <- if (n_low > 0L) {
    paste(words$low_prefix, sample_name_list(d$sample_id[d$.low], more = TRUE))
  } else {
    paste("Lowest:", sample_name_list(d$sample_id))
  }
  depth_sample_panel(d, "library_size", title = words$total,
                     head = head, names = names, cutoff = cutoff)
}

plot_depth_detection <- function(depth, words = DEPTH_WORDS$library) {
  d <- depth_ordered(depth, "n_detected")
  n_feat <- suppressWarnings(
    round(max(d$n_detected) / max(d$detection_rate, na.rm = TRUE)))
  head <- if (is.finite(n_feat)) {
    sprintf("of %s", format(n_feat, big.mark = ","))
  } else {
    "with any signal"
  }
  if (any(d$.low)) head <- paste(head, "\u00b7 amber:", words$low)
  depth_sample_panel(d, "n_detected", title = words$detected,
                     head = head,
                     names = paste("Fewest:", sample_name_list(d$sample_id)))
}

# One per-sample depth panel. `d` is sorted lowest first, with `.low`
# set. `head` is the subtitle; `names` its second line, used only by
# the ranked curve, where the samples are not named on the axis.
depth_sample_panel <- function(d, col, title, head, names, cutoff = NULL) {
  x_scale <- ggplot2::scale_x_continuous(
    labels = depth_axis_labels,
    limits = c(0, NA),
    expand = ggplot2::expansion(mult = c(0, 0.03)))
  line <- missing_cutoff_line(cutoff)

  if (nrow(d) > SAMPLE_MAX_NAMED_BARS) {
    d$rank <- seq_len(nrow(d))
    return(
      ggplot2::ggplot(d, ggplot2::aes(x = .data[[col]], y = .data$rank,
                                      colour = .data$.low)) +
        line +
        ggplot2::geom_point(size = 1.3, alpha = 0.8) +
        missing_over_scale("colour") +
        sample_rank_scale() +
        x_scale +
        ggplot2::labs(title = title, subtitle = paste0(head, "\n", names),
                      x = NULL, y = "Rank") +
        theme_omicsCore()
    )
  }

  # Reversed, because a discrete y axis is drawn bottom-up and the
  # lowest sample belongs at the top.
  d$sample_id <- factor(d$sample_id, levels = rev(d$sample_id))
  ggplot2::ggplot(d, ggplot2::aes(x = .data[[col]], y = .data$sample_id,
                                  fill = .data$.low)) +
    ggplot2::geom_col(width = 0.7) +
    line +
    missing_over_scale("fill") +
    ggplot2::scale_y_discrete(labels = function(x) truncate_pathway_name(x, 24L)) +
    x_scale +
    # No x title: the panel title says what is counted, and the app's
    # 360 px for two panels leave no height for one.
    ggplot2::labs(title = title, subtitle = head, x = NULL, y = NULL) +
    theme_omicsCore() +
    sample_label_theme(nrow(d))
}

# 0, 50k, 2.5M. Each value formatted on its own: format() over the
# whole vector gave every label the notation the largest needed, so the
# axis started at "0e+00".
depth_axis_labels <- function(x) {
  vapply(x, function(v) {
    if (is.na(v)) return("")
    a <- abs(v)
    # Proteomics totals run to 1e11 and past: "1e+05M" is not a label.
    if (a >= 1e12) paste0(signif(v / 1e12, 3), "T")
    else if (a >= 1e9) paste0(signif(v / 1e9, 3), "G")
    else if (a >= 1e6) paste0(signif(v / 1e6, 3), "M")
    else if (a >= 1e3) paste0(signif(v / 1e3, 3), "k")
    else format(v, big.mark = ",", scientific = FALSE, trim = TRUE)
  }, character(1))
}
