# Edge-case tests for impute_matrix and winsorize_counts.
# Covers all backends, degenerate inputs, and error paths.

make_test_mat <- function(nrow = 10, ncol = 6, seed = 123) {
  set.seed(seed)
  mat <- matrix(rnorm(nrow * ncol, mean = 10, sd = 2), nrow = nrow, ncol = ncol)
  rownames(mat) <- paste0("gene_", seq_len(nrow))
  colnames(mat) <- paste0("sample_", seq_len(ncol))
  mat
}

# ---- impute_matrix: no missing data -----------------------------------

test_that("impute_matrix returns unchanged when method='none'", {
  mat <- make_test_mat()
  res <- impute_matrix(mat, method = "none")
  expect_equal(res, mat)
})

test_that("impute_matrix returns unchanged when no NAs present", {
  mat <- make_test_mat()
  res <- impute_matrix(mat, method = "min")
  expect_equal(res, mat)
})

# ---- impute_matrix: min method -----------------------------------------

test_that("impute_matrix min replaces NAs with per-feature minimum", {
  mat <- make_test_mat()
  mat[1, c(1, 3)] <- NA_real_
  mat[2, c(2, 5)] <- NA_real_
  res <- impute_matrix(mat, method = "min")
  expect_false(anyNA(res))
  r1_min <- min(mat[1, ], na.rm = TRUE)
  expect_equal(res[1, 1], r1_min)
  expect_equal(res[1, 3], r1_min)
  r2_min <- min(mat[2, ], na.rm = TRUE)
  expect_equal(res[2, 2], r2_min)
  expect_equal(res[2, 5], r2_min)
  expect_equal(dim(res), dim(mat))
  expect_equal(rownames(res), rownames(mat))
  expect_equal(colnames(res), colnames(mat))
})

test_that("MinDet imputes at the bottom of the sample, not of the feature", {
  # The substantive difference from the old per-feature `half_min`: a
  # detection limit is a property of the run, not of the protein. So the
  # value comes from a low quantile of the *column* it sits in, and two
  # NAs in the same feature but different samples get different values.
  mat <- make_test_mat()
  mat[3, c(4, 6)] <- NA_real_
  res <- impute_matrix(mat, method = "MinDet")

  expect_false(anyNA(res))
  for (j in c(4L, 6L)) {
    obs <- mat[, j][!is.na(mat[, j])]
    expect_lte(res[3, j], stats::median(obs))
    expect_lte(abs(res[3, j] - min(obs)), 0.1 * diff(range(obs)))
  }
  # Deterministic, unlike MinProb -- same input, same answer.
  expect_identical(res, impute_matrix(mat, method = "MinDet"))
})

test_that("impute_matrix min handles features with no non-NA values", {
  mat <- make_test_mat()
  mat[5, ] <- NA_real_
  res <- impute_matrix(mat, method = "min")
  expect_false(anyNA(res))
  expect_equal(as.vector(res[5, ]), rep(0, ncol(mat)))
})

# ---- MAR vs MNAR --------------------------------------------------------

test_that("MAR methods stay inside the observed range, MNAR go below it", {
  # This is the whole reason the control groups them. knn infers the
  # missing value from samples where the protein *was* seen; MinProb
  # assumes it is missing because it was too low to see. Picking the
  # wrong one is not a rounding difference.
  set.seed(7)
  mat <- matrix(stats::rnorm(120, 20, 2), nrow = 15,
                dimnames = list(paste0("g", 1:15), paste0("s", 1:8)))
  mat[sample(length(mat), 25)] <- NA_real_
  obs_min <- min(mat, na.rm = TRUE)

  mnar <- impute_matrix(mat, method = "MinProb")
  mar  <- impute_matrix(mat, method = "knn")
  filled <- is.na(mat)

  expect_lt(min(mnar[filled]), obs_min)
  expect_gte(min(mar[filled]), obs_min)
})

test_that("every method fills every NA and keeps the shape", {
  set.seed(8)
  mat <- matrix(stats::rnorm(120, 20, 2), nrow = 15,
                dimnames = list(paste0("g", 1:15), paste0("s", 1:8)))
  mat[sample(length(mat), 25)] <- NA_real_
  for (m in setdiff(IMPUTE_METHODS, "none")) {
    res <- impute_matrix(mat, method = m)
    expect_false(anyNA(res), info = m)
    expect_identical(dimnames(res), dimnames(mat), info = m)
  }
})

test_that("draws are reproducible, so a report and its script agree", {
  set.seed(9)
  mat <- matrix(stats::rnorm(80, 20, 2), nrow = 10,
                dimnames = list(paste0("g", 1:10), paste0("s", 1:8)))
  mat[sample(length(mat), 15)] <- NA_real_
  for (m in c("MinProb", "QRILC", "man")) {
    expect_identical(impute_matrix(mat, method = m),
                     impute_matrix(mat, method = m), info = m)
  }
})

# ---- impute_matrix: error for missing backend packages -----------------

# ---- impute_matrix: edge cases -----------------------------------------

test_that("impute_matrix works on single-row matrix", {
  mat <- matrix(1:4, nrow = 1, dimnames = list("g1", paste0("s", 1:4)))
  mat[1, 2] <- NA_real_
  res <- impute_matrix(mat, method = "min")
  expect_false(anyNA(res))
  expect_equal(dim(res), c(1, 4))
})

test_that("impute_matrix works on single-column matrix", {
  mat <- matrix(1:8, ncol = 1, dimnames = list(paste0("g", 1:8), "s1"))
  mat[3, 1] <- NA_real_
  res <- impute_matrix(mat, method = "min")
  expect_false(anyNA(res))
  expect_equal(dim(res), c(8, 1))
})

test_that("impute_matrix preserves integer-like values with min method", {
  mat <- matrix(c(5L, NA, 10L, 2L, NA, 8L), nrow = 2, byrow = TRUE,
                dimnames = list(c("g1", "g2"), c("s1", "s2", "s3")))
  res <- impute_matrix(mat, method = "min")
  expect_false(anyNA(res))
  # g1 min = 5 (column 1)
  expect_equal(res[1, 2], min(mat[1, ], na.rm = TRUE))
})

test_that("impute_matrix validates method argument", {
  mat <- make_test_mat()
  expect_error(impute_matrix(mat, method = "xyz"), "should be one of")
})

test_that("impute_matrix converts data.frame to matrix", {
  mat <- make_test_mat()
  df <- as.data.frame(mat)
  res <- impute_matrix(df, method = "none")
  expect_true(is.matrix(res))
})

# ---- winsorize_counts --------------------------------------------------

test_that("winsorize_counts returns list with expected components", {
  mat <- make_test_mat()
  res <- winsorize_counts(mat, k = 20)
  expect_type(res, "list")
  expect_true(all(c("count_mat", "stats", "n_clipped", "n_genes_affected", "k") %in% names(res)))
  expect_equal(res$k, 20)
  expect_equal(dim(res$count_mat), dim(mat))
  expect_equal(rownames(res$count_mat), rownames(mat))
  expect_equal(colnames(res$count_mat), colnames(mat))
})

test_that("winsorize_counts clips extreme values", {
  mat <- make_test_mat()
  # Introduce an extreme outlier in row 1
  mat[1, 1] <- 1e6
  res <- winsorize_counts(mat, k = 3)
  expect_true(max(res$count_mat[1, ]) < 1e6)
  expect_true(res$n_clipped > 0)
})

test_that("winsorize_counts handles NAs transparently", {
  mat <- make_test_mat()
  mat[1, 2] <- NA_real_
  mat[3, 5] <- NA_real_
  res <- winsorize_counts(mat, k = 20)
  expect_equal(which(is.na(res$count_mat)), which(is.na(mat)))
})

test_that("winsorize_counts on constant row does not error", {
  mat <- matrix(rep(5, 20), nrow = 2, dimnames = list(c("g1", "g2"), paste0("s", 1:10)))
  res <- winsorize_counts(mat, k = 20)
  expect_equal(dim(res$count_mat), dim(mat))
  expect_equal(res$n_clipped, 0)
})

test_that("winsorize_counts with low k clips more", {
  mat <- make_test_mat()
  mat[1, 1] <- 1e6
  res_low <- winsorize_counts(mat, k = 3)
  res_high <- winsorize_counts(mat, k = 50)
  expect_true(res_low$n_clipped >= res_high$n_clipped)
})

# ---- winsorize_counts: edge cases --------------------------------------

test_that("winsorize_counts works on single-row matrix", {
  mat <- matrix(1:6, nrow = 1, dimnames = list("g1", paste0("s", 1:6)))
  res <- winsorize_counts(mat, k = 20)
  expect_equal(dim(res$count_mat), c(1, 6))
})

test_that("winsorize_counts works on single-column matrix", {
  mat <- matrix(1:10, ncol = 1, dimnames = list(paste0("g", 1:10), "s1"))
  res <- winsorize_counts(mat, k = 20)
  expect_equal(dim(res$count_mat), c(10, 1))
})

test_that("winsorize_counts preserves dimnames", {
  mat <- make_test_mat()
  res <- winsorize_counts(mat, k = 20)
  expect_equal(rownames(res$count_mat), rownames(mat))
  expect_equal(colnames(res$count_mat), colnames(mat))
})

test_that("winsorize_counts stats df has expected columns", {
  mat <- make_test_mat()
  res <- winsorize_counts(mat, k = 20)
  expect_s3_class(res$stats, "data.frame")
  expect_true(all(c("feature_id", "q1", "q3", "iqr", "threshold", "n_clipped") %in% colnames(res$stats)))
})
