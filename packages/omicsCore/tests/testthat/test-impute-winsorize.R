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

# ---- winsorize_counts: expressed genes only, whole numbers -------------

test_that("winsorize_counts leaves an on/off gene's expressed samples alone", {
  # Expressed in 4 of 20 samples, zero in the rest. Over all counts its
  # quartiles are both 0, and the old rule clipped every expressed sample
  # to zero.
  m <- rbind(on_off = c(rep(0, 16), 480, 510, 495, 530),
             steady = c(100, 105, 98, 102, 99, 101, 103, 97, 100, 104,
                        96, 100, 99, 101, 102, 98, 100, 103, 97, 101))
  colnames(m) <- paste0("s", 1:20)
  storage.mode(m) <- "integer"
  res <- winsorize_counts(m, k = 3)
  expect_identical(res$count_mat["on_off", ], m["on_off", ])
  expect_identical(res$stats$n_clipped[res$stats$feature_id == "on_off"], 0L)
  expect_true(res$stats$winsorized[res$stats$feature_id == "on_off"])

  legacy <- winsorize_counts(m, k = 3, legacy = TRUE)
  expect_true(all(legacy$count_mat["on_off", ] == 0))
})

test_that("winsorize_counts computes the bound from the non-zero counts", {
  set.seed(3)
  x <- c(rep(0, 20), rpois(19, 200), 50000)
  m <- matrix(as.integer(x), nrow = 1, dimnames = list("g", paste0("s", seq_along(x))))
  res <- winsorize_counts(m, k = 5)
  nz <- x[x > 0]
  q <- stats::quantile(nz, c(0.25, 0.75), names = FALSE)
  bound <- q[2] + 5 * diff(q)
  expect_equal(res$stats$threshold, round(bound, 2))
  expect_identical(res$n_clipped, 1L)
  expect_identical(res$count_mat[1, 40], as.integer(floor(bound)))
  # Zeros stay zero
  expect_true(all(res$count_mat[1, 1:20] == 0L))
})

test_that("winsorize_counts keeps counts whole and integer", {
  set.seed(4)
  m <- matrix(rpois(200, 40), nrow = 10,
              dimnames = list(paste0("g", 1:10), paste0("s", 1:20)))
  m[1, 1] <- 100000L
  m[2, 5] <- 77777L
  res <- winsorize_counts(m, k = 3)
  expect_true(is.integer(res$count_mat))
  expect_gte(res$n_clipped, 2L)
  # A double matrix of counts gets whole numbers back too
  md <- m
  storage.mode(md) <- "double"
  resd <- winsorize_counts(md, k = 3)
  expect_true(all(resd$count_mat == round(resd$count_mat)))
  expect_equal(resd$count_mat, res$count_mat + 0)
  # Every clipped value lies at or under its bound
  st <- res$stats
  for (i in which(st$n_clipped > 0)) {
    expect_lte(max(res$count_mat[i, ]), st$threshold[i])
  }
})

test_that("winsorize_counts skips genes expressed in too few samples", {
  m <- rbind(rare = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 9000),
             flat = c(10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 900))
  colnames(m) <- paste0("s", 1:11)
  res <- winsorize_counts(m, k = 3, min_expressed = 3)
  expect_identical(res$count_mat, m)
  expect_identical(res$stats$winsorized, c(FALSE, FALSE))
  expect_true(all(is.na(res$stats$threshold)))
  expect_identical(res$stats$n_nonzero, c(2L, 11L))
  # A lower floor lets the rare gene through; the flat gene has no spread
  # to measure an outlier against and stays untouched.
  res1 <- winsorize_counts(m, k = 3, min_expressed = 2)
  expect_identical(res1$stats$winsorized, c(TRUE, FALSE))
})

test_that("winsorize_counts legacy = TRUE is the old rule, unrounded", {
  set.seed(5)
  m <- matrix(rpois(120, 30), nrow = 6,
              dimnames = list(paste0("g", 1:6), paste0("s", 1:20)))
  m[1, 1] <- 5000L
  m[2, 1:16] <- 0L
  res <- winsorize_counts(m, k = 2, legacy = TRUE)
  q <- apply(m, 1, stats::quantile, probs = c(0.25, 0.75))
  thr <- q[2, ] + 2 * (q[2, ] - q[1, ])
  expect_equal(res$stats$threshold, unname(round(thr, 2)))
  expected <- m + 0
  for (i in seq_len(nrow(m))) expected[i, m[i, ] > thr[i]] <- thr[i]
  expect_equal(res$count_mat, expected)
  expect_true(res$legacy)
  # The on/off gene: zero in 16 of 20, so the old bound is 0
  expect_true(all(res$count_mat[2, ] == 0))
})

test_that("winsorize_counts validates its new arguments", {
  m <- matrix(1:6, nrow = 1, dimnames = list("g1", paste0("s", 1:6)))
  expect_error(winsorize_counts(m, min_expressed = 0), "min_expressed")
  expect_error(winsorize_counts(m, legacy = NA), "legacy")
})
