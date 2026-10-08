# The missingness panel: a panel per sample (which sample is bad, and is
# any over the sample cutoff) and a histogram per feature (how many
# features the missing-value filter removes). A figure review found the
# previous design hard to read: one x axis shared by both panels, which
# a single feature missing everywhere stretched to 100%, squeezing
# sample bars of 5-20% against the left edge; sample names that
# overlapped past ~20 samples; and a density curve of feature rates
# that never showed the cutoff the user set or what it removed.

missing_bundle <- function(n_features = 40L, n_samples = 8L, frac = 0.05,
                           all_missing = 0L, ...) {
  set.seed(11L)
  m <- matrix(rnorm(n_features * n_samples, 20, 2),
              nrow = n_features,
              dimnames = list(paste0("f", seq_len(n_features)),
                              sprintf("S%02d", seq_len(n_samples))))
  if (frac > 0) m[sample.int(length(m), ceiling(frac * length(m)))] <- NA_real_
  if (all_missing > 0L) m[seq_len(all_missing), ] <- NA_real_
  meta <- data.frame(group = rep(c("A", "B"), length.out = n_samples),
                     row.names = colnames(m))
  input <- omics_input(m, meta,
                       data.frame(feature_id = rownames(m)),
                       omics_type = "proteomics",
                       assay_type = "normalized_intensity")
  args <- utils::modifyList(list(missing_threshold = 0.9,
                                 outlier_method = "iqr"), list(...))
  do.call(run_qc, c(list(input), args))
}

# The two panels of plot_qc(view = "missing").
missing_panels <- function(bundle) {
  p <- plot_qc(bundle, view = "missing")
  list(sample = p[[1]], feature = p[[2]])
}

x_range <- function(p) ggplot2::ggplot_build(p)$layout$panel_params[[1]]$x.range

vline_x <- function(p) {
  for (l in p$layers) {
    if (inherits(l$geom, "GeomVline")) return(l$data$xintercept %||% l$aes_params$xintercept)
  }
  NULL
}

test_that("the panel is still a ggplot, so every caller keeps working", {
  # patchwork objects subclass ggplot; the app renders this through
  # renderPlot and the report through knitr, and both only need that.
  p <- plot_qc(missing_bundle(), view = "missing")
  expect_s3_class(p, "ggplot")
  expect_s3_class(p, "patchwork")
})

test_that("the panel survives the app's phone-width theme", {
  # fit_to_width() applies its theme to every patch with `&`.
  p <- plot_qc(missing_bundle(), view = "missing") &
    ggplot2::theme(text = ggplot2::element_text(size = 9))
  f <- tempfile(fileext = ".png")
  on.exit(unlink(f))
  grDevices::png(f, width = 300, height = 360, res = 96)
  expect_no_error(print(p))
  grDevices::dev.off()
})

# ---- the sample panel -------------------------------------------------

test_that("samples are ordered worst first", {
  b <- missing_bundle()
  sdf <- b$results$qc_summary$missingness$sample_metrics
  p <- plot_missing_by_sample(sdf)
  # The factor levels are reversed, because a discrete y axis draws
  # bottom-up and the worst sample belongs at the top.
  lv <- levels(p$data$sample_id)
  expect_identical(
    rev(lv),
    sdf$sample_id[order(sdf$missing_rate, decreasing = TRUE)]
  )
})

test_that("every sample gets a bar when there are few of them", {
  sdf <- missing_bundle(n_samples = 8L)$results$qc_summary$missingness$sample_metrics
  p <- plot_missing_by_sample(sdf)
  expect_equal(nrow(p$data), 8L)
  expect_match(p$labels$subtitle, "^8 samples$")
})

big_cohort <- function(n = 100L) {
  data.frame(sample_id = sprintf("S%03d", seq_len(n)),
             missing_rate = seq(0.5, 0.01, length.out = n))
}

test_that("a large cohort keeps every sample, as a ranked curve", {
  # Bars past a few dozen are thinner than their own labels. Truncating
  # to the worst N answered "which sample is bad" but threw away "how
  # bad is this cohort", which is the other half of the question.
  p <- plot_missing_by_sample(big_cohort(100L))
  expect_equal(nrow(p$data), 100L)
  expect_true("rank" %in% names(p$data))
})

test_that("the ranked curve names the worst samples in its subtitle", {
  # On the plot their labels would land on top of each other: the worst
  # few sit at almost the same rank.
  p <- plot_missing_by_sample(big_cohort(100L))
  expect_match(p$labels$subtitle, "100 samples, ranked", fixed = TRUE)
  expect_match(p$labels$subtitle, "Worst: S001", fixed = TRUE)
  expect_false(grepl("S100", p$labels$subtitle, fixed = TRUE))  # the best
})

test_that("long sample names are listed whole, fewer of them", {
  # Cut short, "Patient_005_plasma_rep1" lost the part that tells
  # samples apart; the list stays short enough for a phone-width panel.
  ids <- sprintf("Patient_%03d_plasma_rep1", 1:20)
  expect_identical(missing_name_list(ids), "Patient_001_plasma_rep1")
  expect_identical(missing_name_list(ids, more = TRUE),
                   "Patient_001_plasma_rep1 and 19 more")
  expect_identical(missing_name_list(c("P01", "P02", "P03", "P04")),
                   "P01, P02, P03")
})

test_that("rank 1 is the worst sample", {
  p <- plot_missing_by_sample(big_cohort(100L))
  worst <- p$data[p$data$rank == 1L, ]
  expect_equal(worst$missing_rate, max(p$data$missing_rate))
})

test_that("the switch happens at the documented threshold, not near it", {
  expect_false("rank" %in% names(
    plot_missing_by_sample(big_cohort(MISSING_MAX_NAMED_SAMPLES))$data))
  expect_true("rank" %in% names(
    plot_missing_by_sample(big_cohort(MISSING_MAX_NAMED_SAMPLES + 1L))$data))
})

test_that("sample names stop overlapping: named bars only while they fit", {
  # 24 samples drew 24 names into ~90 px of the app's 360 px card. The
  # cap is set by what fits there at 7 pt, and the names shrink to 7 pt
  # before it is reached.
  expect_lte(MISSING_MAX_NAMED_SAMPLES, 10L)
  expect_lte(MISSING_SMALL_LABEL_FROM, MISSING_MAX_NAMED_SAMPLES)
  p <- plot_missing_by_sample(big_cohort(MISSING_MAX_NAMED_SAMPLES))
  expect_equal(p$theme$axis.text.y$size, 7)
  # Few samples keep the theme's size.
  p <- plot_missing_by_sample(big_cohort(4L))
  expect_null(p$theme$axis.text.y$size)
})

test_that("the sample axis is not stretched by features missing everywhere", {
  # The regression: one feature missing in every sample took the shared
  # axis to 100%, and sample bars of a few percent became slivers.
  b <- missing_bundle(n_features = 60L, frac = 0.05, all_missing = 3L)
  fm <- b$results$qc_summary$missingness$feature_metrics
  sm <- b$results$qc_summary$missingness$sample_metrics
  expect_equal(max(fm$missing_rate), 1)
  panels <- missing_panels(b)
  upper <- x_range(panels$sample)[2]
  expect_gte(upper, max(sm$missing_rate))
  expect_lt(upper, 0.5)
  # ...while the feature panel shows the whole 0-100%.
  fr <- x_range(panels$feature)
  expect_lte(fr[1], 0)
  expect_gte(fr[2], 1)
})

test_that("samples over the sample cutoff are coloured and counted", {
  b <- missing_bundle(n_samples = 8L, frac = 0.1,
                      sample_missing_threshold = 0.1)
  miss <- b$results$qc_summary$missingness
  expect_gt(length(miss$flagged_samples), 0L)
  expect_lt(length(miss$flagged_samples), 8L)
  sp <- missing_panels(b)$sample
  over <- as.character(sp$data$sample_id[sp$data$over])
  expect_setequal(over, miss$flagged_samples)
  expect_equal(vline_x(sp), 0.1)
  expect_match(sp$labels$subtitle,
               sprintf("%d over the 10%% cutoff", length(miss$flagged_samples)),
               fixed = TRUE)
  # Amber for over, the usual blue for the rest -- not the up/down red.
  built <- ggplot2::ggplot_build(sp)$data[[1]]
  expect_setequal(unique(built$fill), c(MISSING_FILL, MISSING_OVER_FILL))
  expect_identical(MISSING_OVER_FILL, omics_colors$conc_up_down)
  expect_false(MISSING_OVER_FILL %in% c(omics_colors$up, omics_colors$down))
})

test_that("the ranked curve names the samples over the cutoff", {
  df <- big_cohort(40L)
  p <- plot_missing_by_sample(df, cutoff = 0.45)
  n_over <- sum(df$missing_rate > 0.45)
  expect_equal(sum(p$data$over), n_over)
  expect_match(p$labels$subtitle, sprintf("%d over the 45%% cutoff", n_over),
               fixed = TRUE)
  expect_match(p$labels$subtitle, "Over it: S001", fixed = TRUE)
})

test_that("a sample cutoff far past every sample does not stretch the axis", {
  # The line would take the axis back to 80% to show that nothing is
  # near it; the subtitle says so instead.
  sdf <- data.frame(sample_id = paste0("S", 1:6),
                    missing_rate = c(0.05, 0.08, 0.1, 0.12, 0.14, 0.15))
  p <- plot_missing_by_sample(sdf, cutoff = 0.8)
  expect_null(vline_x(p))
  expect_lt(x_range(p)[2], 0.5)
  expect_match(p$labels$subtitle, "none over the 80% cutoff", fixed = TRUE)
  # A cutoff near the samples is drawn, and in range.
  p <- plot_missing_by_sample(sdf, cutoff = 0.2)
  expect_equal(vline_x(p), 0.2)
  expect_gt(x_range(p)[2], 0.2)
})

# ---- the feature panel ------------------------------------------------

test_that("features are a histogram of counts over the whole 0-100%", {
  fdf <- missing_bundle()$results$qc_summary$missingness$feature_metrics
  p <- plot_missing_by_feature(fdf, cutoff = 0.5, n_removed = 0L, n_samples = 8L)
  expect_true(inherits(p$layers[[1]]$geom, "GeomCol"))
  expect_equal(sum(p$data$n), nrow(fdf))
  r <- x_range(p)
  expect_lte(r[1], 0)
  expect_gte(r[2], 1)
})

test_that("bins hold whole missing-rate values, one per possible k / n", {
  # With n samples a rate can only be k / n. Bins that ignore that split
  # some values across two bars and leave others empty.
  rate <- (0:24) / 24
  bins <- missing_feature_bins(rep(rate, 3), cutoff = 0.5, n_samples = 24L)
  expect_equal(nrow(bins), 25L)
  expect_equal(sort(bins$x), rate)
  expect_true(all(bins$n == 3L))
  # No bin edge falls on a possible value.
  edges <- missing_bins(rate, 24L)
  expect_false(any(abs(outer(edges, rate, "-")) < 1e-9))
})

test_that("many samples share bins, still whole values, capped", {
  rate <- (0:60) / 60
  edges <- missing_bins(rate, 60L)
  expect_lte(length(edges) - 1L, MISSING_MAX_BINS + 1L)
  expect_false(any(abs(outer(edges, rate, "-")) < 1e-9))
  bins <- missing_feature_bins(rate, cutoff = 0.5, n_samples = 60L)
  # 0% alone, then pairs -- so every bin is on one side of the cutoff.
  expect_equal(bins$n[bins$x == 0], 1L)
  expect_false(any(duplicated(bins$x)))
})

test_that("the denominator is found from the rates themselves", {
  # Group rates are fractions of each group: 5 and 6 samples give k/30.
  expect_equal(missing_rate_denominator(c(0, 1 / 5, 1 / 6, 1)), 30L)
  # The sample count is preferred over a smaller d that also fits.
  expect_equal(missing_rate_denominator(c(0, 0.25, 0.5), n_samples = 24L), 24L)
  expect_equal(missing_rate_denominator(numeric(0)), 1L)
  expect_true(is.na(missing_rate_denominator(c(0, 1 / 7919), max_d = 100L)))
})

test_that("the cutoff is drawn at the threshold the run used", {
  b <- missing_bundle(n_features = 80L, frac = 0.3, missing_threshold = 0.3)
  fp <- missing_panels(b)$feature
  expect_equal(vline_x(fp), 0.3)
  labels <- vapply(fp$layers, function(l) {
    if (inherits(l$geom, "GeomText")) as.character(l$aes_params$label %||% l$data$label)
    else NA_character_
  }, character(1))
  expect_true("30% cutoff" %in% trimws(labels))
})

test_that("the removed count is the number QC removed, and the bars agree", {
  b <- missing_bundle(n_features = 80L, frac = 0.3, missing_threshold = 0.3)
  miss <- b$results$qc_summary$missingness
  n_flagged <- length(miss$flagged_features)
  expect_gt(n_flagged, 0L)
  # flagged features are exactly the ones dropped from the cleaned layer.
  expect_equal(n_flagged,
               b$input_info$n_features_in - b$input_info$n_features_out)
  fp <- missing_panels(b)$feature
  expect_match(fp$labels$subtitle,
               sprintf("80 features · %d removed\n(missing in more than 30%% of samples)",
                       n_flagged), fixed = TRUE)
  # The amber bars hold exactly those features.
  expect_equal(sum(fp$data$n[fp$data$removed]), n_flagged)
  expect_true(all(fp$data$x[fp$data$removed] > 0.3))
  expect_true(all(fp$data$x[!fp$data$removed] <= 0.3))
})

test_that("nothing removed says so", {
  b <- missing_bundle(frac = 0.05, missing_threshold = 0.9)
  fp <- missing_panels(b)$feature
  expect_match(fp$labels$subtitle, "none removed", fixed = TRUE)
  expect_false(any(fp$data$removed))
})

group_bundle <- function(filter) {
  set.seed(5L)
  grp <- rep(c("ctrl", "trt"), c(5L, 6L))
  m <- matrix(rnorm(60L * 11L, 20, 2), 60L,
              dimnames = list(sprintf("P%02d", 1:60), sprintf("S%02d", 1:11)))
  m[1:15, grp == "ctrl"] <- NA        # absent in one group only
  m[16:20, ] <- NA                    # absent everywhere
  m[21:30, c(1, 2, 3, 6, 7, 8, 9)] <- NA  # >50% missing in both groups
  m[sample.int(length(m), 40L)] <- NA
  meta <- data.frame(group = grp, row.names = colnames(m))
  input <- omics_input(m, meta, data.frame(feature_id = rownames(m)),
                       omics_type = "proteomics",
                       assay_type = "normalized_intensity")
  run_qc(input, missing_threshold = 0.5, missing_filter = filter,
         group_col = "group", outlier_method = "iqr")
}

test_that("a group filter plots the rate it compared with the cutoff", {
  # The overall rate of a protein absent from the controls is ~45%, under
  # the line -- but the any_group rule keeps it for a different reason,
  # and the all_groups rule removes it although its overall rate is
  # under 50%. Only the group rate puts every feature on its true side.
  for (filter in c("any_group", "all_groups")) {
    b <- group_bundle(filter)
    miss <- b$results$qc_summary$missingness
    fm <- miss$feature_metrics
    fp <- missing_panels(b)$feature
    expect_equal(sum(fp$data$n[fp$data$removed]),
                 sum(fm$filter_missing_rate > 0.5), info = filter)
    expect_equal(sum(fp$data$n[fp$data$removed]),
                 length(miss$flagged_features), info = filter)
    # The overall rate would have told a different story.
    expect_false(sum(fm$missing_rate > 0.5) == length(miss$flagged_features),
                 info = filter)
    expect_equal(vline_x(fp), 0.5)
    expect_match(fp$labels$subtitle,
                 sprintf("%d removed", length(miss$flagged_features)),
                 fixed = TRUE)
  }
  any_p <- missing_panels(group_bundle("any_group"))$feature
  expect_identical(any_p$labels$x, "Lowest missing rate among groups")
  expect_match(any_p$labels$subtitle, "missing in every group", fixed = TRUE)
  all_p <- missing_panels(group_bundle("all_groups"))$feature
  expect_identical(all_p$labels$x, "Highest missing rate among groups")
  expect_match(all_p$labels$subtitle, "missing in at least one group",
               fixed = TRUE)
})

test_that("a group filter without its group rates draws no false line", {
  # An older bundle with the rule but not the rates: the overall rate is
  # all there is to draw, and the cutoff over it would not be the line
  # the filter used. The count is still the one QC removed.
  b <- group_bundle("any_group")
  b$results$qc_summary$missingness$feature_metrics$filter_missing_rate <- NULL
  fp <- missing_panels(b)$feature
  expect_null(vline_x(fp))
  expect_false(any(fp$data$removed))
  expect_match(fp$labels$subtitle,
               sprintf("%d removed",
                       length(b$results$qc_summary$missingness$flagged_features)),
               fixed = TRUE)
})

test_that("an old bundle without settings or flags draws no line or count", {
  b <- missing_bundle(n_features = 80L, frac = 0.3, missing_threshold = 0.3)
  b$results$qc_summary$missingness$settings <- NULL
  b$results$qc_summary$missingness$flagged_features <- NULL
  b$params <- list()
  b$input_info$n_features_out <- NULL
  panels <- missing_panels(b)
  expect_null(vline_x(panels$feature))
  expect_null(vline_x(panels$sample))
  expect_identical(panels$feature$labels$subtitle, "80 features")
})

test_that("the removed count falls back to the input counts", {
  b <- missing_bundle(n_features = 80L, frac = 0.3, missing_threshold = 0.3)
  n <- length(b$results$qc_summary$missingness$flagged_features)
  b$results$qc_summary$missingness$flagged_features <- NULL
  expect_identical(missing_rules(b)$n_removed, as.integer(n))
})

test_that("all-complete data still draws, and says so", {
  # The good case -- nothing missing. The density it replaces had no
  # spread to estimate from and needed a special case.
  b <- missing_bundle(frac = 0)
  panels <- missing_panels(b)
  fp <- panels$feature
  expect_s3_class(fp, "ggplot")
  expect_equal(nrow(fp$data), 1L)
  expect_equal(fp$data$x, 0)
  expect_equal(fp$data$n, 40L)
  expect_match(fp$labels$subtitle, "40 features · no missing values",
               fixed = TRUE)
  expect_match(fp$labels$subtitle, "none removed", fixed = TRUE)
  expect_equal(vline_x(fp), 0.9)
  expect_match(panels$sample$labels$subtitle, "no missing values", fixed = TRUE)
  f <- tempfile(fileext = ".png")
  on.exit(unlink(f))
  grDevices::png(f, width = 620, height = 360, res = 96)
  expect_no_error(print(plot_qc(b, view = "missing")))
  grDevices::dev.off()
})

test_that("a single feature, or none, does not error", {
  p <- plot_missing_by_feature(data.frame(feature_id = "f1", missing_rate = 0.2))
  expect_s3_class(p, "ggplot")
  p <- plot_missing_by_feature(data.frame(feature_id = character(0),
                                          missing_rate = numeric(0)))
  expect_s3_class(p, "ggplot")
})

test_that("count breaks spread out on the square-root axis", {
  expect_equal(missing_count_breaks(2450), c(0, 200, 1000, 2000))
  expect_equal(missing_count_breaks(3), c(0, 1, 2))
  expect_equal(missing_count_breaks(0), 0)
})

# ---- the sample axis --------------------------------------------------

test_that("the axis is anchored at zero and never runs past 100%", {
  expect_equal(missing_axis_upper(c(0, 0.02)), 0.05)      # floor
  expect_equal(missing_axis_upper(c(0.9, 1)), 1)          # ceiling
  expect_gt(missing_axis_upper(c(0.1, 0.3)), 0.3)         # headroom
  expect_lte(missing_axis_upper(c(0.1, 0.3)), 1)
})

test_that("an empty or all-NA set of rates still gives a usable range", {
  expect_equal(missing_axis_upper(numeric(0)), 1)
  expect_equal(missing_axis_upper(c(NA_real_, NA_real_)), 1)
})
