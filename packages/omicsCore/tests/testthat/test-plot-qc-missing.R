# The missingness panel used to be one 30-bin histogram faceted over
# sample and feature. With a dozen samples most bins were empty and the
# rest read as unexplained spikes, and a histogram throws away the one
# thing the sample panel exists to answer: which sample is bad.

missing_bundle <- function(n_features = 40L, n_samples = 8L, frac = 0.05) {
  set.seed(11L)
  m <- matrix(rnorm(n_features * n_samples, 20, 2),
              nrow = n_features,
              dimnames = list(paste0("f", seq_len(n_features)),
                              sprintf("S%02d", seq_len(n_samples))))
  if (frac > 0) m[sample.int(length(m), ceiling(frac * length(m)))] <- NA_real_
  meta <- data.frame(group = rep(c("A", "B"), length.out = n_samples),
                     row.names = colnames(m))
  input <- omics_input(m, meta,
                       data.frame(feature_id = rownames(m)),
                       omics_type = "proteomics",
                       assay_type = "normalized_intensity")
  run_qc(input, missing_threshold = 0.9, outlier_method = "iqr")
}

test_that("the panel is still a ggplot, so every caller keeps working", {
  # patchwork objects subclass ggplot; the app renders this through
  # renderPlot and the report through knitr, and both only need that.
  p <- plot_qc(missing_bundle(), view = "missing")
  expect_s3_class(p, "ggplot")
  expect_s3_class(p, "patchwork")
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

test_that("a large cohort is truncated, and says so", {
  # Silently dropping rows would be worse than a cramped panel; the
  # subtitle is what makes the truncation honest.
  sdf <- data.frame(sample_id = sprintf("S%03d", 1:100),
                    missing_rate = seq(0.5, 0.01, length.out = 100))
  p <- plot_missing_by_sample(sdf)
  expect_equal(nrow(p$data), MISSING_MAX_SAMPLE_BARS)
  expect_match(p$labels$subtitle, "worst 30 of 100 samples", fixed = TRUE)
  # And it is the worst 30 that survive, not the first 30.
  expect_true(all(p$data$missing_rate >= sort(sdf$missing_rate,
                                              decreasing = TRUE)[30]))
})

# ---- the feature panel ------------------------------------------------

test_that("features are drawn as a density over the possible range", {
  fdf <- missing_bundle()$results$qc_summary$missingness$feature_metrics
  p <- plot_missing_by_feature(fdf)
  # The curve is evaluated on [0, 1]: a missing rate cannot be negative,
  # and an unbounded smooth would draw one.
  expect_gte(min(p$data$x), 0)
  expect_lte(max(p$data$x), 1)
  expect_true(all(p$data$y >= 0))
})

test_that("all-complete data does not error on a zero bandwidth", {
  # The good case -- nothing missing -- gives density() no spread to
  # work from. It has to say so rather than fail.
  fdf <- data.frame(feature_id = paste0("f", 1:10), missing_rate = rep(0, 10))
  p <- plot_missing_by_feature(fdf)
  expect_s3_class(p, "ggplot")
  expect_match(p$labels$subtitle, "all at 0%", fixed = TRUE)
})

test_that("a single feature does not error", {
  p <- plot_missing_by_feature(
    data.frame(feature_id = "f1", missing_rate = 0.2))
  expect_s3_class(p, "ggplot")
})

# ---- the shared axis --------------------------------------------------

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

test_that("both panels are given the same range", {
  # They measure the same quantity. facet_wrap(scales = "free") gave
  # them separate ones, so two panels that looked alike could be an
  # order of magnitude apart.
  b <- missing_bundle()
  miss <- b$results$qc_summary$missingness
  upper <- missing_axis_upper(c(miss$sample_metrics$missing_rate,
                                miss$feature_metrics$missing_rate))
  expect_gte(upper, max(miss$sample_metrics$missing_rate))
  expect_gte(upper, max(miss$feature_metrics$missing_rate))
})
