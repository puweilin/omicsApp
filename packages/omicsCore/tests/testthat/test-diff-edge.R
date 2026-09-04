# Edge-case and robustness tests for run_diff and run_diff_continuous.
# These tests use ttest/lm backends (pure R) to avoid Bioconductor Suggests.

make_diff_input <- function(n_feat = 20, n_samp = 8) {
  mat <- matrix(rnorm(n_feat * n_samp, mean = 10, sd = 2),
                nrow = n_feat, ncol = n_samp)
  rownames(mat) <- paste0("gene_", seq_len(n_feat))
  colnames(mat) <- paste0("sample_", seq_len(n_samp))
  meta <- data.frame(
    group = rep(c("control", "treatment"), each = n_samp / 2),
    batch = rep(c("B1", "B2"), n_samp / 2),
    row.names = colnames(mat)
  )
  feat <- data.frame(
    feature_id = rownames(mat),
    feature_name = paste0("Gene", seq_len(n_feat)),
    stringsAsFactors = FALSE
  )
  omics_input(mat, meta, feat, omics_type = "proteomics", assay_type = "normalized_intensity")
}

# ---- run_diff: basic execution -----------------------------------------

test_that("run_diff with ttest backend returns an analysis_bundle", {
  inp <- make_diff_input()
  res <- run_diff(inp, method = "ttest", analysis_type = "group",
                  group_col = "group",
                  control_group = "control", case_group = "treatment")
  expect_true(is_analysis_bundle(res))
})

test_that("run_diff returns standardized result columns", {
  inp <- make_diff_input()
  res <- run_diff(inp, method = "ttest", analysis_type = "group",
                  group_col = "group",
                  control_group = "control", case_group = "treatment")
  df <- res$results$diff_result_df
  expect_true("feature_id" %in% colnames(df))
  expect_true("effect" %in% colnames(df))
  expect_true("adj_p_value" %in% colnames(df))
  expect_equal(nrow(df), nrow(inp$expr_mat))
})

test_that("run_diff includes input_info and params", {
  inp <- make_diff_input()
  res <- run_diff(inp, method = "ttest", analysis_type = "group",
                  group_col = "group",
                  control_group = "control", case_group = "treatment")
  expect_equal(res$input_info$omics_type, "proteomics")
  expect_equal(res$params$method, "ttest")
  expect_equal(res$params$control_group, "control")
  expect_equal(res$params$case_group, "treatment")
})

test_that("run_diff with lm backend returns an analysis_bundle", {
  inp <- make_diff_input()
  res <- run_diff(inp, method = "lm", analysis_type = "group",
                  group_col = "group",
                  control_group = "control", case_group = "treatment")
  expect_true(is_analysis_bundle(res))
})

test_that("run_diff with auto method selects a backend", {
  inp <- make_diff_input()
  res <- run_diff(inp, method = "auto", analysis_type = "group",
                  group_col = "group",
                  control_group = "control", case_group = "treatment")
  expect_true(is_analysis_bundle(res))
})

# ---- run_diff: covariates -----------------------------------------------

test_that("run_diff with covariates runs successfully", {
  inp <- make_diff_input()
  res <- run_diff(inp, method = "lm", analysis_type = "group",
                  group_col = "group",
                  control_group = "control", case_group = "treatment",
                  covariates = "batch")
  expect_true(is_analysis_bundle(res))
})

# ---- run_diff: error paths ---------------------------------------------

test_that("run_diff errors on unsupported method for analysis_type", {
  inp <- make_diff_input()
  # edgeR only supports group analysis
  skip_if_not(requireNamespace("edgeR", quietly = TRUE),
              "edgeR not installed; can't test this path")
  expect_error(
    run_diff(inp, method = "edger", analysis_type = "continuous",
             continuous_col = "batch"),
    "edgeR.*only supports"
  )
})

test_that("run_diff errors when group_col missing from meta", {
  inp <- make_diff_input()
  expect_error(
    run_diff(inp, method = "ttest", analysis_type = "group",
             group_col = "nonexistent",
             control_group = "control", case_group = "treatment"),
    "group_col"
  )
})

test_that("run_diff errors when control_group not found in meta", {
  inp <- make_diff_input()
  expect_error(
    run_diff(inp, method = "ttest", analysis_type = "group",
             group_col = "group",
             control_group = "not_a_group", case_group = "treatment"),
    "Each group must have at least"
  )
})

test_that("run_diff errors on invalid input type", {
  expect_error(
    run_diff(list(), method = "ttest", analysis_type = "group",
             group_col = "x", control_group = "a", case_group = "b"),
    "omics_input"
  )
})

test_that("run_diff errors when control equals case", {
  inp <- make_diff_input()
  expect_error(
    run_diff(inp, method = "ttest", analysis_type = "group",
             group_col = "group",
             control_group = "control", case_group = "control"),
    "duplicated|identical|same"
  )
})

test_that("run_diff errors on invalid analysis_type", {
  inp <- make_diff_input()
  expect_error(
    run_diff(inp, method = "ttest", analysis_type = "invalid",
             group_col = "group",
             control_group = "control", case_group = "treatment"),
    "should be one of"
  )
})

# ---- run_diff_continuous: basic execution -------------------------------

test_that("run_diff_continuous with lm backend returns an analysis_bundle", {
  inp <- make_diff_input()
  inp$meta_df$age <- c(25, 30, 35, 40, 45, 50, 55, 60)
  res <- run_diff_continuous(inp, method = "lm", continuous_col = "age")
  expect_true(is_analysis_bundle(res))
})

test_that("run_diff_continuous errors without continuous_col", {
  inp <- make_diff_input()
  expect_error(
    run_diff_continuous(inp, method = "lm", continuous_col = "nonexistent"),
    "continuous_col"
  )
})

# ---- filter_diff_results -----------------------------------------------

test_that("filter_diff_results filters by effect and adj_p", {
  inp <- make_diff_input()
  bundle <- run_diff(inp, method = "ttest", analysis_type = "group",
                     group_col = "group",
                     control_group = "control", case_group = "treatment")
  df <- bundle$results$diff_result_df
  res <- filter_diff_results(df, p_cutoff = 1, p_preference = "adjusted",
                              effect_cutoff = 0)
  expect_equal(nrow(res), nrow(inp$expr_mat))
})

test_that("filter_diff_results with strict cutoff returns fewer features", {
  inp <- make_diff_input()
  bundle <- run_diff(inp, method = "ttest", analysis_type = "group",
                     group_col = "group",
                     control_group = "control", case_group = "treatment")
  df <- bundle$results$diff_result_df
  res_strict <- filter_diff_results(df, p_cutoff = 1e-100,
                                     effect_cutoff = 100)
  expect_equal(nrow(res_strict), 0)
})

test_that("filter_diff_results validates result_df schema", {
  expect_error(filter_diff_results(list()), "Missing required")
})

# ---- is_analysis_bundle -------------------------------------------------

test_that("is_analysis_bundle returns TRUE for valid bundles", {
  inp <- make_diff_input()
  bundle <- run_diff(inp, method = "ttest", analysis_type = "group",
                     group_col = "group",
                     control_group = "control", case_group = "treatment")
  expect_true(is_analysis_bundle(bundle))
})

test_that("is_analysis_bundle returns FALSE for other objects", {
  expect_false(is_analysis_bundle(list()))
  expect_false(is_analysis_bundle(NULL))
  expect_false(is_analysis_bundle("string"))
})

# ---- new_analysis_bundle -----------------------------------------------

test_that("new_analysis_bundle creates a valid bundle", {
  b <- new_analysis_bundle("diff", list(x = 1),
                           params = list(method = "test"))
  expect_true(is_analysis_bundle(b))
  expect_equal(b$analysis_name, "diff")
  expect_equal(b$params$method, "test")
})

# ---- check_diff_result_schema ------------------------------------------

test_that("check_diff_result_schema passes for valid result df", {
  inp <- make_diff_input()
  bundle <- run_diff(inp, method = "ttest", analysis_type = "group",
                     group_col = "group",
                     control_group = "control", case_group = "treatment")
  df <- bundle$results$diff_result_df
  expect_true(check_diff_result_schema(df))
})

test_that("check_diff_result_schema errors when columns missing", {
  expect_error(check_diff_result_schema(data.frame(x = 1)), "Missing required")
})
