# The leave-one-out outlier test: the one method that can flag a sample
# in a study of ten or fewer, where a z-score over the cohort cannot
# reach 3. Its thresholds were set by simulation (the rates are in the
# change that introduced it); the fast version below guards them.

# Log2 intensities shaped like a small proteomics study: a dynamic range
# of features, two groups with a tenth of the features shifted, a noise
# level that differs a little between samples and grows at low
# abundance, and left-censored missing values.
sim_loo_matrix <- function(n, p = 500L, groups = 2L) {
  mu <- stats::rnorm(p, 23, 2.5)
  grp <- rep(seq_len(groups), length.out = n)
  eff <- matrix(0, p, groups)
  for (g in seq_len(groups)[-1L]) {
    de <- sample(p, p / 10)
    eff[de, g] <- stats::rnorm(length(de), 0, 1)
  }
  sig_s <- 0.3 * exp(stats::rnorm(n, 0, 0.2))
  sig_f <- 1 + 0.5 * pmax(0, (23 - mu) / 2.5)
  x <- vapply(seq_len(n), function(s) {
    mu + eff[, grp[s]] + stats::rnorm(p, 0, sig_s[s] * sig_f)
  }, numeric(p))
  lod <- stats::quantile(x, 0.1)
  x[x < lod + stats::rnorm(length(x), 0, 0.5)] <- NA
  dimnames(x) <- list(paste0("f", seq_len(p)), paste0("s", seq_len(n)))
  x
}

loo_input <- function(x, group = rep(c("A", "B"), length.out = ncol(x))) {
  omics_input(x, data.frame(group = group, row.names = colnames(x)),
              data.frame(feature_id = rownames(x), row.names = rownames(x)),
              omics_type = "proteomics", assay_type = "normalized_intensity")
}

test_that("clean small studies flag nothing; broken samples are flagged", {
  clean_flags <- 0L
  scrambled_hits <- 0L
  swapped_hits <- 0L
  runs <- 0L
  for (n in 4:10) for (seed in 1:12) {
    set.seed(seed * 100 + n)
    x <- sim_loo_matrix(n)
    runs <- runs + 1L
    clean_flags <- clean_flags + (length(qc_outliers_loo(x, 3)$flagged_samples) > 0L)
    # A failed run: every intensity off by 1.5 log2, up or down.
    broken <- x
    broken[, 1] <- broken[, 1] + sample(c(-1.5, 1.5), nrow(x), replace = TRUE)
    scrambled_hits <- scrambled_hits + ("s1" %in% qc_outliers_loo(broken, 3)$flagged_samples)
    # A sample of something else: 40% of the proteome at other levels.
    swapped <- x
    idx <- sample(nrow(x), 0.4 * nrow(x))
    swapped[idx, 1] <- swapped[idx, 1] + stats::rnorm(length(idx), 0, 1.5)
    swapped_hits <- swapped_hits + ("s1" %in% qc_outliers_loo(swapped, 3)$flagged_samples)
  }
  expect_lte(clean_flags / runs, 0.05)
  expect_gte(scrambled_hits / runs, 0.95)
  expect_gte(swapped_hits / runs, 0.85)
})

test_that("group structure is not mistaken for an outlier", {
  # Three vs three with a strong difference: every sample has a close
  # partner in its own group.
  set.seed(7)
  mu <- stats::rnorm(400, 23, 2.5)
  shift <- c(rep(3, 120), rep(0, 280))
  x <- cbind(replicate(3, mu + stats::rnorm(400, 0, 0.3)),
             replicate(3, mu + shift + stats::rnorm(400, 0, 0.3)))
  dimnames(x) <- list(paste0("f", 1:400), paste0("s", 1:6))
  res <- qc_outliers(loo_input(x, rep(c("A", "B"), each = 3)), method = "loo")
  expect_identical(res$flagged_samples, character(0))
})

test_that("loo reports what it compared, as correlations", {
  set.seed(11)
  x <- sim_loo_matrix(6)
  x[, 2] <- x[, 2] + sample(c(-1.5, 1.5), nrow(x), replace = TRUE)
  res <- qc_outliers(loo_input(x), method = "loo")
  expect_identical(res$method, "loo")
  expect_named(res$stats, c("sample_id", "nearest_sample", "nearest_correlation",
                            "reference_correlation", "z_score", "is_outlier"))
  expect_identical(res$flagged_samples, "s2")
  s2 <- res$stats[res$stats$sample_id == "s2", ]
  expect_lt(s2$nearest_correlation, s2$reference_correlation)
  expect_gt(s2$z_score, 3)
  # Its threshold is the shared one.
  expect_identical(qc_outliers(loo_input(x), method = "loo",
                               sd_threshold = 1e6)$flagged_samples, character(0))
})

test_that("loo needs four samples and says so", {
  set.seed(3)
  x <- sim_loo_matrix(3)
  res <- qc_outliers(loo_input(x), method = "loo")
  expect_identical(res$flagged_samples, character(0))
  expect_true(all(is.na(res$stats$z_score)))
  expect_match(res$note, "at least 4 samples", all = FALSE)
})

test_that("the default runs all four methods and leaves the other three unchanged", {
  set.seed(5)
  x <- sim_loo_matrix(8)
  x[, 1] <- x[, 1] + sample(c(-1.5, 1.5), nrow(x), replace = TRUE)
  inp <- loo_input(x)
  all4 <- qc_outliers(inp)
  expect_identical(all4$method, c("pca", "connectivity", "iqr", "loo"))
  expect_named(all4$by_method, c("pca", "connectivity", "iqr", "loo"))
  for (m in c("pca", "connectivity", "iqr")) {
    alone <- qc_outliers(inp, method = m)
    expect_identical(all4$by_method[[m]]$stats, alone$stats)
    expect_identical(all4$by_method[[m]]$flagged_samples, alone$flagged_samples)
  }
  three <- qc_outliers(inp, method = c("pca", "connectivity", "iqr"))
  expect_identical(three$by_method, all4$by_method[c("pca", "connectivity", "iqr")])
  # With eight samples the z-score tests cannot reach 3; loo can.
  expect_true("s1" %in% all4$flagged_samples)
  expect_identical(sort(unique(all4$stats$method)),
                   sort(c("pca", "connectivity", "iqr", "loo")))
  expect_match(all4$note, "leave-one-out test still applies", all = FALSE)
})

test_that("run_qc accepts loo and keeps its per-modality defaults", {
  set.seed(9)
  x <- sim_loo_matrix(6)
  x[, 3] <- x[, 3] + sample(c(-1.5, 1.5), nrow(x), replace = TRUE)
  inp <- loo_input(x)
  b <- run_qc(inp, outlier_method = c("pca", "connectivity", "iqr", "loo"),
              impute_method = "none")
  expect_true("s3" %in% b$results$qc_summary$outliers$flagged_samples)
  expect_true("s3" %in% colnames(qc_cleaned_input(b, inp)$expr_mat))
  expect_identical(run_qc(inp, impute_method = "none")$params$outlier_method, "pca")
})
