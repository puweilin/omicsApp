# QC defects found in the 2026-10 audit.

aq_counts <- function(seed = 3L) {
  set.seed(seed)
  ids <- paste0("S", 1:12)
  mu <- exp(stats::rnorm(400, log(100), 1.5))
  m <- matrix(stats::rnbinom(400 * 12, mu = mu, size = 30), 400,
              dimnames = list(paste0("G", 1:400), ids))
  # One library sequenced four times deeper: the same sample, scaled.
  m[, 1] <- m[, 1] * 4L
  omics_input(m, data.frame(group = rep(c("A", "B"), each = 6), row.names = ids),
              data.frame(feature_id = rownames(m), row.names = rownames(m)),
              omics_type = "rnaseq", assay_type = "raw_count")
}

aq_prot <- function(linear = FALSE, seed = 4L) {
  set.seed(seed)
  ids <- paste0("S", 1:8)
  m <- matrix(stats::rnorm(200 * 8, 22, 1), 200, dimnames = list(paste0("P", 1:200), ids))
  m[sample(length(m), 150)] <- NA
  if (linear) m <- 2^m
  omics_input(m, data.frame(group = rep(c("A", "B"), each = 4), row.names = ids),
              data.frame(feature_id = rownames(m), row.names = rownames(m)),
              omics_type = "proteomics",
              assay_type = if (linear) "raw_intensity" else "normalized_intensity")
}

test_that("a deeper library is not an outlier: counts are tested as log2-CPM", {
  inp <- aq_counts()
  res <- qc_outliers(inp, method = c("pca", "connectivity"), sd_threshold = 3)
  expect_false("S1" %in% res$flagged_samples)
  expect_match(res$note, "log2-CPM", all = FALSE)
})

test_that("flagged samples are kept unless removal is asked for", {
  inp <- aq_prot()
  inp$expr_mat[, 1] <- inp$expr_mat[, 1] + 10
  b <- run_qc(inp, outlier_method = "iqr", outlier_sd_threshold = 1.5,
              impute_method = "none")
  expect_identical(b$results$qc_summary$recommended_filters$remove_samples, "S1")
  expect_true("S1" %in% colnames(qc_cleaned_input(b, inp)$expr_mat))
  expect_match(b$warnings, "kept", all = FALSE)

  b2 <- run_qc(inp, outlier_method = "iqr", outlier_sd_threshold = 1.5,
               impute_method = "none", remove_outliers = TRUE)
  expect_false("S1" %in% colnames(qc_cleaned_input(b2, inp)$expr_mat))
})

test_that("a z-score threshold no sample size can reach says so", {
  inp <- aq_prot()
  res <- qc_outliers(inp, method = "pca", sd_threshold = 3)
  expect_match(res$note, "cannot exceed", all = FALSE)
})

test_that("imputing linear intensities gives positive values on the same scale", {
  skip_if_not_installed("imputeLCMD")
  inp <- aq_prot(linear = TRUE)
  na <- is.na(inp$expr_mat)
  b <- run_qc(inp, outlier_method = "none", impute_method = "MinProb",
              missing_threshold = 1)
  out <- qc_cleaned_input(b, inp)$expr_mat
  expect_false(anyNA(out))
  expect_true(all(out[na] > 0))
  # Left-censored: below what was observed, but on the same scale.
  expect_lt(stats::median(out[na]), stats::median(out[!na]))
  expect_gt(stats::median(out[na]), stats::median(out[!na]) / 100)
  expect_identical(qc_cleaned_input(b, inp)$assay_type, "raw_intensity")
  expect_match(b$warnings, "log2 scale", all = FALSE)
})

test_that("imputed log intensities are labelled, and the plot compares observed with imputed", {
  inp <- aq_prot()
  b <- run_qc(inp, outlier_method = "none", impute_method = "min",
              missing_threshold = 1)
  expect_identical(qc_cleaned_input(b, inp)$assay_type, "imputed_intensity")
  imp <- b$results$qc_summary$imputation
  expect_identical(imp$n_imputed, 150L)
  p <- plot_qc(b, view = "imputation")
  expect_setequal(unique(p$data$type), c("observed", "imputed"))
})

test_that("knn imputation takes a neighbour count", {
  skip_if_not_installed("imputeLCMD")
  m <- aq_prot()$expr_mat
  expect_false(anyNA(impute_matrix(m, method = "knn", K = 5)))
  expect_false(anyNA(impute_matrix(m, method = "knn", k = 5)))
})

test_that("normalize_omics() names a missing assay type instead of crashing", {
  inp <- aq_prot(linear = TRUE)
  inp$assay_type <- NULL
  expect_error(normalize_omics(inp, method = "log2"), "assay_type")
})
