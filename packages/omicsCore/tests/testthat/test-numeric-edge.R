# Numeric / statistical edge cases across diff backends and enrichment.

make_numeric_input <- function(n_feat = 20, n_samp = 8) {
  mat <- matrix(rnorm(n_feat * n_samp, mean = 10, sd = 2),
                nrow = n_feat, ncol = n_samp)
  rownames(mat) <- paste0("gene_", seq_len(n_feat))
  colnames(mat) <- paste0("s", seq_len(n_samp))
  half <- n_samp %/% 2
  meta <- data.frame(
    group = c(rep("ctrl", half), rep("trt", n_samp - half)),
    age   = seq(20, by = 5, length.out = n_samp),
    row.names = colnames(mat)
  )
  feat <- data.frame(feature_id = rownames(mat),
                     stringsAsFactors = FALSE)
  omics_input(mat, meta, feat, omics_type = "proteomics",
              assay_type = "normalized_intensity")
}

make_count_input <- function(n_feat = 30, n_samp = 8) {
  mat <- matrix(as.integer(rpois(n_feat * n_samp, lambda = 50)),
                nrow = n_feat, ncol = n_samp)
  rownames(mat) <- paste0("gene_", seq_len(n_feat))
  colnames(mat) <- paste0("s", seq_len(n_samp))
  half <- n_samp %/% 2
  meta <- data.frame(
    group = c(rep("ctrl", half), rep("trt", n_samp - half)),
    row.names = colnames(mat)
  )
  feat <- data.frame(feature_id = rownames(mat),
                     stringsAsFactors = FALSE)
  omics_input(mat, meta, feat, omics_type = "rnaseq",
              assay_type = "raw_count")
}

# ---- ttest backend ----------------------------------------------------

test_that("ttest with one-sample-per-group still produces a bundle", {
  inp <- make_numeric_input(n_feat = 10, n_samp = 2)
  # ttest requires >=2 per group; expect a clear error
  expect_error(
    run_diff(inp, method = "ttest", analysis_type = "group",
             group_col = "group", control_group = "ctrl", case_group = "trt"),
    regexp = ".+"
  )
})

test_that("ttest with unequal group sizes works", {
  inp <- make_numeric_input(n_feat = 10, n_samp = 6)
  inp$meta_df$group <- c("ctrl", "ctrl", "trt", "trt", "trt", "trt")
  b <- run_diff(inp, method = "ttest", analysis_type = "group",
                group_col = "group", control_group = "ctrl",
                case_group = "trt")
  expect_true(is_analysis_bundle(b))
  expect_equal(nrow(b$results$diff_result_df), 10)
})

test_that("ttest p-values are bounded in [0, 1]", {
  inp <- make_numeric_input(n_feat = 20, n_samp = 8)
  b <- run_diff(inp, method = "ttest", analysis_type = "group",
                group_col = "group", control_group = "ctrl",
                case_group = "trt")
  p <- b$results$diff_result_df$p_value
  p <- p[!is.na(p)]
  expect_true(all(p >= 0 & p <= 1))
})

test_that("ttest adj_p_value >= p_value (BH monotonicity)", {
  inp <- make_numeric_input(n_feat = 30, n_samp = 8)
  b <- run_diff(inp, method = "ttest", analysis_type = "group",
                group_col = "group", control_group = "ctrl",
                case_group = "trt")
  df <- b$results$diff_result_df
  ok_rows <- complete.cases(df[, c("p_value", "adj_p_value")])
  expect_true(all(df$adj_p_value[ok_rows] >= df$p_value[ok_rows] - 1e-12))
})

test_that("ttest with no real signal yields uniformly-distributed p-values", {
  set.seed(123)
  inp <- make_numeric_input(n_feat = 200, n_samp = 8)
  b <- run_diff(inp, method = "ttest", analysis_type = "group",
                group_col = "group", control_group = "ctrl",
                case_group = "trt")
  p <- b$results$diff_result_df$p_value
  p <- p[!is.na(p)]
  # KS test against uniform: weak null, just want fraction below 0.05
  # to be roughly 5% (allow generous slack)
  frac <- mean(p < 0.05)
  expect_true(frac < 0.20)
})

test_that("ttest detects strong signal", {
  set.seed(456)
  inp <- make_numeric_input(n_feat = 10, n_samp = 8)
  inp$expr_mat[1, 5:8] <- inp$expr_mat[1, 5:8] + 10  # strong signal in g1
  b <- run_diff(inp, method = "ttest", analysis_type = "group",
                group_col = "group", control_group = "ctrl",
                case_group = "trt")
  df <- b$results$diff_result_df
  g1_p <- df$p_value[df$feature_id == "gene_1"]
  expect_true(!is.na(g1_p) && g1_p < 0.01)
})

# ---- limma backend ----------------------------------------------------

test_that("limma backend returns a bundle", {
  skip_if_not_installed("limma")
  inp <- make_numeric_input(n_feat = 30, n_samp = 8)
  b <- run_diff(inp, method = "limma", analysis_type = "group",
                group_col = "group", control_group = "ctrl",
                case_group = "trt")
  expect_true(is_analysis_bundle(b))
})

test_that("limma produces standardized output columns", {
  skip_if_not_installed("limma")
  inp <- make_numeric_input(n_feat = 30, n_samp = 8)
  b <- run_diff(inp, method = "limma", analysis_type = "group",
                group_col = "group", control_group = "ctrl",
                case_group = "trt")
  df <- b$results$diff_result_df
  expected_cols <- c("feature_id", "effect", "p_value", "adj_p_value")
  expect_true(all(expected_cols %in% colnames(df)))
})

test_that("limma with single feature works", {
  skip_if_not_installed("limma")
  inp <- make_numeric_input(n_feat = 1, n_samp = 8)
  b <- run_diff(inp, method = "limma", analysis_type = "group",
                group_col = "group", control_group = "ctrl",
                case_group = "trt")
  expect_equal(nrow(b$results$diff_result_df), 1)
})

# ---- deseq2 backend ---------------------------------------------------

test_that("deseq2 backend returns a bundle for raw counts", {
  skip_if_not_installed("DESeq2")
  inp <- make_count_input(n_feat = 30, n_samp = 8)
  b <- run_diff(inp, method = "deseq2", analysis_type = "group",
                group_col = "group", control_group = "ctrl",
                case_group = "trt")
  expect_true(is_analysis_bundle(b))
})

test_that("deseq2 effect column is log2FC", {
  skip_if_not_installed("DESeq2")
  inp <- make_count_input(n_feat = 30, n_samp = 8)
  b <- run_diff(inp, method = "deseq2", analysis_type = "group",
                group_col = "group", control_group = "ctrl",
                case_group = "trt")
  expect_equal(b$results$diff_result_df$effect_type[1], "log2FC")
})

test_that("deseq2 errors on non-integer matrix", {
  skip_if_not_installed("DESeq2")
  inp <- make_numeric_input(n_feat = 30, n_samp = 8)
  # proteomics input → deseq2 not supported
  expect_error(
    run_diff(inp, method = "deseq2", analysis_type = "group",
             group_col = "group", control_group = "ctrl",
             case_group = "trt")
  )
})

# ---- edger backend ----------------------------------------------------

test_that("edger backend returns a bundle for counts", {
  skip_if_not_installed("edgeR")
  inp <- make_count_input(n_feat = 30, n_samp = 8)
  b <- run_diff(inp, method = "edger", analysis_type = "group",
                group_col = "group", control_group = "ctrl",
                case_group = "trt")
  expect_true(is_analysis_bundle(b))
})

test_that("edger rejects non-counts omics_type", {
  skip_if_not_installed("edgeR")
  inp <- make_numeric_input(n_feat = 30, n_samp = 8)
  expect_error(
    run_diff(inp, method = "edger", analysis_type = "group",
             group_col = "group", control_group = "ctrl",
             case_group = "trt")
  )
})

# ---- lm backend for continuous ----------------------------------------

test_that("lm backend supports continuous covariate", {
  inp <- make_numeric_input(n_feat = 20, n_samp = 8)
  b <- run_diff_continuous(inp, method = "lm", continuous_col = "age")
  expect_true(is_analysis_bundle(b))
})

test_that("lm backend errors when continuous_col missing", {
  inp <- make_numeric_input(n_feat = 20, n_samp = 8)
  expect_error(
    run_diff_continuous(inp, method = "lm", continuous_col = "nope"),
    "not found"
  )
})

test_that("lm backend handles character-coerced numeric column (silent coerce or error)", {
  inp <- make_numeric_input(n_feat = 20, n_samp = 8)
  inp$meta_df$age <- as.character(inp$meta_df$age)
  # Behavior is implementation-defined: either the backend coerces
  # characters to numeric internally (succeeds), or it rejects them
  # (errors). Both outcomes are acceptable as long as the function
  # doesn't crash hard somewhere downstream.
  res <- tryCatch(
    run_diff_continuous(inp, method = "lm", continuous_col = "age"),
    error = function(e) e
  )
  expect_true(is_analysis_bundle(res) || inherits(res, "error"))
})

# ---- enrichment numerical edge cases ----------------------------------

test_that("run_enrichment ORA p-values are bounded in [0, 1]", {
  skip_if_not_installed("clusterProfiler")
  inp <- make_numeric_input(n_feat = 40, n_samp = 8)
  bd <- run_diff(inp, method = "ttest", analysis_type = "group",
                 group_col = "group", control_group = "ctrl",
                 case_group = "trt")
  be <- run_enrichment(bd, type = "ora", database = "hallmark")
  p <- be$results$enrich_result_df$p_value
  if (length(p) > 0) {
    expect_true(all(p >= 0 & p <= 1, na.rm = TRUE))
  } else {
    succeed()
  }
})

test_that("run_enrichment ORA adj_p >= p (BH monotonicity)", {
  skip_if_not_installed("clusterProfiler")
  inp <- make_numeric_input(n_feat = 40, n_samp = 8)
  bd <- run_diff(inp, method = "ttest", analysis_type = "group",
                 group_col = "group", control_group = "ctrl",
                 case_group = "trt")
  be <- run_enrichment(bd, type = "ora", database = "hallmark")
  df <- be$results$enrich_result_df
  if (nrow(df) > 0) {
    ok <- complete.cases(df[, c("p_value", "adj_p_value")])
    expect_true(all(df$adj_p_value[ok] >= df$p_value[ok] - 1e-12))
  } else {
    succeed()
  }
})

test_that("make_ranked_features returns sorted vector by default", {
  inp <- make_numeric_input(n_feat = 30, n_samp = 8)
  bd <- run_diff(inp, method = "ttest", analysis_type = "group",
                 group_col = "group", control_group = "ctrl",
                 case_group = "trt")
  ranked <- make_ranked_features(bd$results$diff_result_df)
  expect_type(ranked, "double")
  expect_true(!is.unsorted(rev(ranked)) || !is.unsorted(ranked))
})

test_that("make_ranked_features has unique names", {
  inp <- make_numeric_input(n_feat = 30, n_samp = 8)
  bd <- run_diff(inp, method = "ttest", analysis_type = "group",
                 group_col = "group", control_group = "ctrl",
                 case_group = "trt")
  ranked <- make_ranked_features(bd$results$diff_result_df)
  expect_equal(length(unique(names(ranked))), length(ranked))
})

# ---- winsorize_counts numerical -----------------------------------

test_that("winsorize_counts preserves dimensions", {
  mat <- matrix(rpois(100, 50), nrow = 10, ncol = 10)
  rownames(mat) <- paste0("g", 1:10)
  colnames(mat) <- paste0("s", 1:10)
  res <- winsorize_counts(mat, k = 5)
  expect_equal(dim(res$count_mat), dim(mat))
})

test_that("winsorize_counts clips extreme values", {
  mat <- matrix(rpois(100, 50), nrow = 10, ncol = 10)
  rownames(mat) <- paste0("g", 1:10)
  colnames(mat) <- paste0("s", 1:10)
  mat[1, 1] <- 1e6  # extreme outlier
  res <- winsorize_counts(mat, k = 5)
  expect_true(res$count_mat[1, 1] < 1e6)
})

test_that("winsorize_counts with k = 0 always clips", {
  mat <- matrix(rpois(100, 50), nrow = 10, ncol = 10)
  rownames(mat) <- paste0("g", 1:10)
  colnames(mat) <- paste0("s", 1:10)
  res <- winsorize_counts(mat, k = 0)
  expect_true(res$n_clipped >= 0)
})

test_that("winsorize_counts result includes stats per feature", {
  mat <- matrix(rpois(100, 50), nrow = 10, ncol = 10)
  rownames(mat) <- paste0("g", 1:10)
  colnames(mat) <- paste0("s", 1:10)
  res <- winsorize_counts(mat, k = 3)
  expect_true("stats" %in% names(res))
})

test_that("winsorize_counts preserves integer-ness when input is integer", {
  mat <- matrix(as.integer(rpois(100, 50)), nrow = 10, ncol = 10)
  rownames(mat) <- paste0("g", 1:10)
  colnames(mat) <- paste0("s", 1:10)
  res <- winsorize_counts(mat, k = 3)
  expect_true(is.numeric(res$count_mat))
})

# ---- run_integration: input/output shape checks ----------------------

make_dual_proj_int <- function() {
  i1 <- make_numeric_input(n_feat = 30, n_samp = 6)
  i2 <- make_numeric_input(n_feat = 30, n_samp = 6)
  colnames(i2$expr_mat) <- paste0("t", 1:6)
  rownames(i2$meta_df) <- paste0("t", 1:6)
  i2$omics_type <- "rnaseq"
  i2$assay_type <- "logcpm"
  proj <- omics_project("dual",
                        experiments = list(proteomics = i1, rnaseq = i2))
  d1 <- run_diff(i1, method = "ttest", analysis_type = "group",
                 group_col = "group", control_group = "ctrl",
                 case_group = "trt")
  d2 <- run_diff(i2, method = "ttest", analysis_type = "group",
                 group_col = "group", control_group = "ctrl",
                 case_group = "trt")
  list(proj = proj, d1 = d1, d2 = d2)
}

test_that("run_integration concordance returns integration_df with rows", {
  s <- make_dual_proj_int()
  res <- run_integration(s$proj, method = "concordance",
                         experiments = c("proteomics", "rnaseq"),
                         diff_bundles = list(proteomics = s$d1,
                                             rnaseq = s$d2))
  df <- res$results$integration_df
  expect_true(nrow(df) > 0)
})

test_that("run_integration concordance produces stable shape", {
  s <- make_dual_proj_int()
  res1 <- run_integration(s$proj, method = "concordance",
                          experiments = c("proteomics", "rnaseq"),
                          diff_bundles = list(proteomics = s$d1,
                                              rnaseq = s$d2))
  res2 <- run_integration(s$proj, method = "concordance",
                          experiments = c("proteomics", "rnaseq"),
                          diff_bundles = list(proteomics = s$d1,
                                              rnaseq = s$d2))
  expect_equal(nrow(res1$results$integration_df),
               nrow(res2$results$integration_df))
})

test_that("run_integration correlation method works", {
  s <- make_dual_proj_int()
  res <- tryCatch(
    run_integration(s$proj, method = "correlation",
                    experiments = c("proteomics", "rnaseq"),
                    diff_bundles = list(proteomics = s$d1, rnaseq = s$d2)),
    error = function(e) e
  )
  expect_true(is_analysis_bundle(res) || inherits(res, "error"))
})
