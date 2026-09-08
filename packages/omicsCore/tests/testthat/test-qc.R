make_qc_input <- function(n_features = 8, n_samples = 8, missing_frac = 0.1,
                          omics_type = "proteomics") {
  set.seed(2024)
  feat_ids <- paste0("g", seq_len(n_features))
  samp_ids <- paste0("s", seq_len(n_samples))
  expr <- matrix(
    rnorm(n_features * n_samples, mean = 5, sd = 1.5),
    nrow = n_features,
    dimnames = list(feat_ids, samp_ids)
  )
  n_na <- floor(length(expr) * missing_frac)
  na_idx <- sample(length(expr), n_na)
  expr[na_idx] <- NA
  meta <- data.frame(
    group = rep(c("A", "B"), length.out = n_samples),
    row.names = samp_ids,
    stringsAsFactors = FALSE
  )
  feat <- data.frame(
    feature_id = feat_ids,
    row.names = feat_ids,
    stringsAsFactors = FALSE
  )
  omics_input(expr, meta, feat, omics_type = omics_type,
              assay_type = "normalized_intensity")
}

test_that("qc_missingness returns matched sample / feature frames", {
  x <- make_qc_input()
  m <- qc_missingness(x, feature_missing_cutoff = 0.5)
  expect_named(m, c("sample_metrics", "feature_metrics", "flagged_samples",
                    "flagged_features", "settings"))
  expect_equal(nrow(m$sample_metrics), 8L)
  expect_equal(nrow(m$feature_metrics), 8L)
  expect_true(all(m$sample_metrics$missing_rate >= 0))
  expect_true(all(m$sample_metrics$missing_rate <= 1))
})

test_that("qc_missingness flags features above the cutoff", {
  x <- make_qc_input()
  x$expr_mat[1, ] <- NA  # 100% missing
  m <- qc_missingness(x, feature_missing_cutoff = 0.5)
  expect_true("g1" %in% m$flagged_features)
})

test_that("qc_outliers iqr method flags extreme samples", {
  x <- make_qc_input(missing_frac = 0)
  x$expr_mat[, 1] <- x$expr_mat[, 1] + 50   # far above the rest
  res <- qc_outliers(x, method = "iqr", sd_threshold = 1.5)
  expect_equal(res$method, "iqr")
  expect_true("s1" %in% res$flagged_samples)
})

test_that("qc_outliers pca runs and returns z-scores", {
  x <- make_qc_input(missing_frac = 0)
  res <- qc_outliers(x, method = "pca", sd_threshold = 3)
  expect_named(res, c("method", "stats", "flagged_samples"))
  expect_true(all(c("PC1", "PC2", "z_pc1", "z_pc2", "is_outlier") %in%
                  colnames(res$stats)))
})

test_that("qc_outliers connectivity returns correlation-based stats", {
  x <- make_qc_input(missing_frac = 0)
  res <- qc_outliers(x, method = "connectivity", sd_threshold = 3)
  expect_true("mean_correlation" %in% colnames(res$stats))
})

test_that("qc_outliers with multiple methods unions the flags", {
  x <- make_qc_input(missing_frac = 0)
  x$expr_mat[, 1] <- x$expr_mat[, 1] + 100
  res <- qc_outliers(x, method = c("pca", "iqr"), sd_threshold = 2)
  expect_named(res, c("method", "stats", "flagged_samples", "by_method"))
  expect_true("s1" %in% res$flagged_samples)
})

test_that("impute_matrix passes through when method='none'", {
  x <- make_qc_input()
  out <- impute_matrix(x$expr_mat, method = "none")
  expect_equal(sum(is.na(out)), sum(is.na(x$expr_mat)))
})

test_that("impute_matrix min replaces NA with row min", {
  # matrix() fills column-major: row 'a' = c(1, NA, 5); row 'b' = c(2, 4, NA)
  m <- matrix(c(1, 2, NA, 4, 5, NA), nrow = 2,
              dimnames = list(c("a", "b"), c("s1", "s2", "s3")))
  out <- impute_matrix(m, method = "min")
  expect_false(anyNA(out))
  expect_equal(out["a", "s2"], 1)   # row 'a' min = 1
  expect_equal(out["b", "s3"], 2)   # row 'b' min = 2
})

test_that("impute_matrix min fills with the row minimum", {
  # row 'a' = c(2, NA); row 'b' = c(8, 4)
  m <- matrix(c(2, 8, NA, 4), nrow = 2,
              dimnames = list(c("a", "b"), c("s1", "s2")))
  out <- impute_matrix(m, method = "min")
  expect_equal(out["a", "s2"], 2)   # row 'a' min
  expect_false(anyNA(out))
})

test_that("winsorize_counts clips per-gene extremes", {
  set.seed(1)
  m <- matrix(rpois(50, 50), nrow = 5,
              dimnames = list(paste0("g", 1:5), paste0("s", 1:10)))
  m[1, 1] <- 1e6  # one extreme value
  res <- winsorize_counts(m, k = 5)
  expect_equal(dim(res$count_mat), dim(m))
  expect_true(res$n_clipped >= 1L)
  expect_lt(res$count_mat[1, 1], 1e6)
})

# ---- low-expression filter --------------------------------------------

make_counts_input <- function(n_features = 10, n_samples = 6, counts = NULL) {
  ids <- paste0("g", seq_len(n_features))
  sids <- paste0("s", seq_len(n_samples))
  m <- counts %||% matrix(100L, n_features, n_samples,
                          dimnames = list(ids, sids))
  dimnames(m) <- list(ids, sids)
  omics_input(m,
              data.frame(group = rep(c("a", "b"), length.out = n_samples),
                         row.names = sids),
              data.frame(feature_id = ids, row.names = ids),
              omics_type = "rnaseq", assay_type = "raw_count")
}

test_that("low_expression_features keeps a feature seen often enough", {
  m <- rbind(
    plenty  = c(50L, 50L, 50L, 50L),
    exactly = c(10L, 10L, 10L,  0L),   # 3 samples reach 10
    thin    = c(10L, 10L,  0L,  0L),   # only 2 do
    spike   = c(999L, 0L,  0L,  0L),   # a big total, one sample
    silent  = c(0L,  0L,  0L,  0L)
  )
  colnames(m) <- paste0("s", 1:4)
  dropped <- low_expression_features(m, min_count = 10, min_samples = 3)
  expect_setequal(dropped, c("thin", "spike", "silent"))

  # A row of NA is not a row of detections
  m2 <- rbind(m, allna = rep(NA_integer_, 4))
  expect_true("allna" %in%
                low_expression_features(m2, min_count = 10, min_samples = 3))

  # Off when either number is zero
  expect_length(low_expression_features(m, 0, 3), 0L)
  expect_length(low_expression_features(m, 10, 0), 0L)
})

test_that("run_qc filters low-expression features on raw counts by default", {
  set.seed(4)
  m <- matrix(100L, 10, 6, dimnames = list(paste0("g", 1:10), paste0("s", 1:6)))
  m[9:10, ] <- 0L                       # two features nobody detected
  x <- make_counts_input(counts = m)

  b <- run_qc(x)
  kept <- rownames(b$results$cleaned_input$expr_mat)
  expect_false(any(c("g9", "g10") %in% kept))
  expect_length(kept, 8L)

  # and the numbers used are recorded, not just their effect
  expect_equal(b$params$min_count, 10)
  expect_equal(b$params$min_samples, 3)

  # switched off, nothing is dropped for being quiet
  b0 <- run_qc(x, min_count = 0, min_samples = 0)
  expect_length(rownames(b0$results$cleaned_input$expr_mat), 10L)
})

# 10 reads is a low number on raw counts and a high one on TPM or logCPM,
# where it would drop most of the matrix. The filter is not adapted to those
# scales -- it is withheld from them, and the control is not offered either.
test_that("run_qc leaves non-count layers alone", {
  set.seed(9)
  # Varying, not constant: a matrix of one repeated value makes the
  # correlation in qc_outliers() undefined, and the warning has nothing to
  # do with what is being tested here.
  m <- matrix(rnorm(60, mean = 1), 10, 6,
              dimnames = list(paste0("g", 1:10), paste0("s", 1:6)))
  ids <- rownames(m)
  x <- omics_input(m,
                   data.frame(group = rep(c("a", "b"), 3), row.names = colnames(m)),
                   data.frame(feature_id = ids, row.names = ids),
                   omics_type = "rnaseq", assay_type = "logcpm")
  b <- run_qc(x)
  expect_equal(b$params$min_count, 0)
  expect_length(rownames(b$results$cleaned_input$expr_mat), 10L)
})

test_that("run_qc returns an analysis_bundle with the expected slots", {
  x <- make_qc_input()
  b <- run_qc(x, missing_threshold = 0.5,
              outlier_method = "iqr", outlier_sd_threshold = 3)
  expect_true(is_analysis_bundle(b))
  expect_equal(b$analysis_name, "run_qc")
  expect_named(b$results, c("qc_summary", "cleaned_input"))
  expect_true(is_omics_input(b$results$cleaned_input))
})

test_that("run_qc removes flagged samples and features", {
  x <- make_qc_input(missing_frac = 0)
  x$expr_mat[1, ] <- NA              # feature -> 100% missing
  x$expr_mat[, 1] <- x$expr_mat[, 1] + 100  # sample -> extreme
  b <- run_qc(x, missing_threshold = 0.5,
              outlier_method = "iqr", outlier_sd_threshold = 1.5)
  cleaned <- b$results$cleaned_input
  expect_false("g1" %in% rownames(cleaned$expr_mat))
  expect_false("s1" %in% colnames(cleaned$expr_mat))
})

test_that("run_qc with impute_method='min' fills NAs in cleaned input", {
  x <- make_qc_input(missing_frac = 0.05)
  b <- run_qc(x, missing_threshold = 1,             # keep all features
              outlier_method = "none",
              impute_method = "min")
  expect_false(anyNA(b$results$cleaned_input$expr_mat))
  expect_true(!is.null(b$results$cleaned_input$raw_mat))
})

test_that("run_qc errors when filters wipe out the input", {
  x <- make_qc_input()
  # One missing value in every feature, and a threshold that flags any:
  # nothing survives. (A negative threshold used to be the way to say
  # this; it is now refused as not a probability.)
  x$expr_mat[cbind(seq_len(nrow(x$expr_mat)), 1L)] <- NA
  expect_error(
    run_qc(x, missing_threshold = 0),
    "loosen the thresholds"
  )
})

test_that("plot_qc returns ggplot objects for each view", {
  x <- make_qc_input()
  b <- run_qc(x, missing_threshold = 1, outlier_method = "connectivity",
              impute_method = "min")
  for (v in c("missing", "pca", "connectivity", "imputation")) {
    p <- plot_qc(b, view = v)
    expect_s3_class(p, "ggplot")
  }
})

test_that("plot_qc 'imputation' errors when no imputation was run", {
  x <- make_qc_input()
  b <- run_qc(x, missing_threshold = 1, outlier_method = "none",
              impute_method = "none")
  expect_error(plot_qc(b, view = "imputation"), "no imputation step")
})
