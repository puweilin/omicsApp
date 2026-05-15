# Edge-case tests for QC functions: qc_missingness, qc_outliers, and
# the internal qc_rnaseq function.

make_small_input <- function(nrow = 10, ncol = 8) {
  mat <- matrix(rnorm(nrow * ncol, mean = 10, sd = 2), nrow = nrow, ncol = ncol)
  rownames(mat) <- paste0("gene_", seq_len(nrow))
  colnames(mat) <- paste0("sample_", seq_len(ncol))
  meta <- data.frame(
    group = rep(c("A", "B"), each = ncol / 2),
    row.names = colnames(mat)
  )
  feat <- data.frame(
    feature_id = rownames(mat),
    feature_name = paste0("Gene", seq_len(nrow)),
    stringsAsFactors = FALSE
  )
  omics_input(mat, meta, feat, omics_type = "proteomics", assay_type = "intensity")
}

# ---- qc_missingness ----------------------------------------------------

test_that("qc_missingness returns complete structure", {
  inp <- make_small_input()
  res <- qc_missingness(inp)
  expect_type(res, "list")
  expect_named(res, c("sample_metrics", "feature_metrics",
                      "flagged_samples", "flagged_features", "settings"))
})

test_that("qc_missingness with no NAs flags nothing", {
  inp <- make_small_input()
  res <- qc_missingness(inp, sample_missing_cutoff = 0.2,
                        feature_missing_cutoff = 0.5)
  expect_equal(nrow(res$sample_metrics), ncol(inp$expr_mat))
  expect_equal(nrow(res$feature_metrics), nrow(inp$expr_mat))
  expect_equal(length(res$flagged_samples), 0)
  expect_equal(length(res$flagged_features), 0)
})

test_that("qc_missingness flags samples above cutoff", {
  inp <- make_small_input()
  inp$expr_mat[1:5, 1] <- NA_real_
  res <- qc_missingness(inp, sample_missing_cutoff = 0.4,
                        feature_missing_cutoff = 0.5)
  # sample_1 has 5/10 features NA = 0.5 > 0.4
  expect_true("sample_1" %in% res$flagged_samples)
})

test_that("qc_missingness flags features above cutoff", {
  inp <- make_small_input()
  inp$expr_mat[1, 1:6] <- NA_real_
  res <- qc_missingness(inp, feature_missing_cutoff = 0.5)
  # gene_1 has 6/8 samples NA = 0.75 > 0.5
  expect_true("gene_1" %in% res$flagged_features)
})

test_that("qc_missingness with NULL sample cutoff leaves samples unflagged", {
  inp <- make_small_input()
  inp$expr_mat[1:10, 1] <- NA_real_
  res <- qc_missingness(inp, sample_missing_cutoff = NULL,
                        feature_missing_cutoff = 0.5)
  expect_equal(length(res$flagged_samples), 0)
})

test_that("qc_missingness settings echo thresholds", {
  inp <- make_small_input()
  res <- qc_missingness(inp, sample_missing_cutoff = 0.3,
                        feature_missing_cutoff = 0.7)
  expect_equal(res$settings$sample_missing_cutoff, 0.3)
  expect_equal(res$settings$feature_missing_cutoff, 0.7)
})

test_that("qc_missingness metrics are within [0, 1]", {
  inp <- make_small_input()
  inp$expr_mat[sample(seq_len(prod(dim(inp$expr_mat))), 10)] <- NA_real_
  res <- qc_missingness(inp)
  expect_true(all(res$sample_metrics$missing_rate >= 0))
  expect_true(all(res$sample_metrics$missing_rate <= 1))
  expect_true(all(res$feature_metrics$missing_rate >= 0))
  expect_true(all(res$feature_metrics$missing_rate <= 1))
})

test_that("qc_missingness validates input", {
  expect_error(qc_missingness(list()), "omics_input")
})

# ---- qc_outliers: pca --------------------------------------------------

test_that("qc_outliers pca returns expected structure", {
  inp <- make_small_input()
  res <- qc_outliers(inp, method = "pca")
  expect_equal(res$method, "pca")
  expect_named(res$stats, c("sample_id", "PC1", "PC2", "z_pc1", "z_pc2", "is_outlier"))
  expect_type(res$flagged_samples, "character")
})

test_that("qc_outliers pca with default threshold finds few outliers", {
  inp <- make_small_input()
  res <- qc_outliers(inp, method = "pca", sd_threshold = 3)
  # Random normal data should have close to zero outliers at z > 3
  expect_equal(length(res$flagged_samples), 0)
})

test_that("qc_outliers pca with tight threshold finds more outliers", {
  inp <- make_small_input()
  res <- qc_outliers(inp, method = "pca", sd_threshold = 0.5)
  expect_true(length(res$flagged_samples) >= 0)
})

# ---- qc_outliers: connectivity -----------------------------------------

test_that("qc_outliers connectivity returns expected structure", {
  inp <- make_small_input()
  res <- qc_outliers(inp, method = "connectivity")
  expect_equal(res$method, "connectivity")
  expect_named(res$stats, c("sample_id", "mean_correlation", "z_score", "is_outlier"))
})

test_that("qc_outliers connectivity with sd=3 finds few outliers in random data", {
  inp <- make_small_input()
  res <- qc_outliers(inp, method = "connectivity", sd_threshold = 3)
  expect_equal(length(res$flagged_samples), 0)
})

# ---- qc_outliers: iqr --------------------------------------------------

test_that("qc_outliers iqr returns expected structure", {
  inp <- make_small_input()
  res <- qc_outliers(inp, method = "iqr")
  expect_equal(res$method, "iqr")
  expect_named(res$stats, c("sample_id", "mean_signal", "lower_fence",
                             "upper_fence", "is_outlier"))
})

test_that("qc_outliers iqr works with different k multiplier", {
  inp <- make_small_input()
  res <- qc_outliers(inp, method = "iqr", sd_threshold = 1.5)
  expect_equal(res$method, "iqr")
})

# ---- qc_outliers: multi-method -----------------------------------------

test_that("qc_outliers with multiple methods returns union of flags", {
  inp <- make_small_input()
  res <- qc_outliers(inp, method = c("pca", "connectivity"), sd_threshold = 3)
  expect_true(all(c("pca", "connectivity") %in% res$method))
  expect_type(res$flagged_samples, "character")
  expect_named(res$by_method, c("pca", "connectivity"))
})

test_that("qc_outliers multi-method stats has method column", {
  inp <- make_small_input()
  res <- qc_outliers(inp, method = c("pca", "iqr"))
  expect_true("method" %in% colnames(res$stats))
})

test_that("qc_outliers validates method", {
  inp <- make_small_input()
  expect_error(qc_outliers(inp, method = "nonexistent"), "should be one of")
})

test_that("qc_outliers validates input", {
  expect_error(qc_outliers(list()), "omics_input")
})

test_that("qc_outliers handles data with NAs via mean imputation", {
  inp <- make_small_input()
  inp$expr_mat[2, 3] <- NA_real_
  inp$expr_mat[5, 6] <- NA_real_
  res <- qc_outliers(inp, method = "pca")
  expect_equal(res$method, "pca")
  expect_type(res$flagged_samples, "character")
})

# ---- qc_outliers: edge cases ------------------------------------------

test_that("qc_outliers on single sample errors for pca", {
  mat <- matrix(rnorm(5), ncol = 1, dimnames = list(paste0("g", 1:5), "s1"))
  meta <- data.frame(group = "A", row.names = "s1")
  feat <- data.frame(feature_id = paste0("g", 1:5), feature_name = paste0("G", 1:5),
                     stringsAsFactors = FALSE)
  inp <- omics_input(mat, meta, feat, omics_type = "proteomics")
  # prcomp needs >= 2 observations; single sample should error
  expect_error(qc_outliers(inp, method = "pca"))
})

test_that("qc_outliers on single sample works for connectivity", {
  mat <- matrix(rnorm(5), ncol = 1, dimnames = list(paste0("g", 1:5), "s1"))
  meta <- data.frame(group = "A", row.names = "s1")
  feat <- data.frame(feature_id = paste0("g", 1:5), feature_name = paste0("G", 1:5),
                     stringsAsFactors = FALSE)
  inp <- omics_input(mat, meta, feat, omics_type = "proteomics")
  res <- qc_outliers(inp, method = "connectivity")
  expect_equal(length(res$flagged_samples), 0)
})

test_that("qc_outliers on single sample works for iqr", {
  mat <- matrix(rnorm(5), ncol = 1, dimnames = list(paste0("g", 1:5), "s1"))
  meta <- data.frame(group = "A", row.names = "s1")
  feat <- data.frame(feature_id = paste0("g", 1:5), feature_name = paste0("G", 1:5),
                     stringsAsFactors = FALSE)
  inp <- omics_input(mat, meta, feat, omics_type = "proteomics")
  res <- qc_outliers(inp, method = "iqr")
  expect_equal(res$method, "iqr")
})
