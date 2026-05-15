# Edge-case tests for run_integration, integration utilities,
# and integration plot functions.

make_dual_project <- function() {
  mat1 <- matrix(rnorm(100, mean = 10, sd = 2), nrow = 20, ncol = 5)
  rownames(mat1) <- paste0("gene_", 1:20)
  colnames(mat1) <- paste0("s", 1:5)
  meta1 <- data.frame(group = c("control", "control", "treatment", "treatment", "treatment"),
                      row.names = paste0("s", 1:5))
  feat1 <- data.frame(feature_id = paste0("gene_", 1:20),
                      feature_name = paste0("Gene", 1:20),
                      stringsAsFactors = FALSE)
  inp1 <- omics_input(mat1, meta1, feat1, omics_type = "proteomics",
                      assay_type = "intensity")

  mat2 <- matrix(rnorm(100, mean = 10, sd = 2), nrow = 20, ncol = 5)
  rownames(mat2) <- paste0("gene_", 1:20)
  colnames(mat2) <- paste0("t", 1:5)
  meta2 <- data.frame(group = c("control", "control", "treatment", "treatment", "treatment"),
                      row.names = paste0("t", 1:5))
  feat2 <- data.frame(feature_id = paste0("gene_", 1:20),
                      feature_name = paste0("Gene", 1:20),
                      stringsAsFactors = FALSE)
  inp2 <- omics_input(mat2, meta2, feat2, omics_type = "rnaseq",
                      assay_type = "normalized_count")

  omics_project("dual", experiments = list(proteomics = inp1, rnaseq = inp2))
}

# ---- run_integration: concordance --------------------------------------

test_that("run_integration concordance returns an analysis_bundle", {
  proj <- make_dual_project()
  diff1 <- run_diff(proj$experiments$proteomics, method = "ttest",
                    analysis_type = "group", group_col = "group",
                    control_group = "control", case_group = "treatment")
  diff2 <- run_diff(proj$experiments$rnaseq, method = "ttest",
                    analysis_type = "group", group_col = "group",
                    control_group = "control", case_group = "treatment")
  bundles <- list(proteomics = diff1, rnaseq = diff2)
  res <- run_integration(proj, method = "concordance",
                         experiments = c("proteomics", "rnaseq"),
                         diff_bundles = bundles)
  expect_true(is_analysis_bundle(res))
  expect_equal(res$analysis_name, "run_integration")
})

test_that("run_integration concordance result has integration_df", {
  proj <- make_dual_project()
  diff1 <- run_diff(proj$experiments$proteomics, method = "ttest",
                    analysis_type = "group", group_col = "group",
                    control_group = "control", case_group = "treatment")
  diff2 <- run_diff(proj$experiments$rnaseq, method = "ttest",
                    analysis_type = "group", group_col = "group",
                    control_group = "control", case_group = "treatment")
  bundles <- list(proteomics = diff1, rnaseq = diff2)
  res <- run_integration(proj, method = "concordance",
                         experiments = c("proteomics", "rnaseq"),
                         diff_bundles = bundles)
  df <- res$results$integration_df
  expect_true("effect" %in% colnames(df) || "feature_id" %in% colnames(df))
})

test_that("run_integration concordance with feature_link works", {
  proj <- make_dual_project()
  diff1 <- run_diff(proj$experiments$proteomics, method = "ttest",
                    analysis_type = "group", group_col = "group",
                    control_group = "control", case_group = "treatment")
  diff2 <- run_diff(proj$experiments$rnaseq, method = "ttest",
                    analysis_type = "group", group_col = "group",
                    control_group = "control", case_group = "treatment")

  fl <- data.frame(
    proteomics = paste0("gene_", 1:20),
    rnaseq = paste0("gene_", 1:20),
    stringsAsFactors = FALSE
  )
  proj$feature_link <- fl
  bundles <- list(proteomics = diff1, rnaseq = diff2)
  res <- run_integration(proj, method = "concordance",
                         experiments = c("proteomics", "rnaseq"),
                         diff_bundles = bundles)
  expect_true(is_analysis_bundle(res))
})

# ---- run_integration: error paths --------------------------------------

test_that("run_integration errors on invalid project", {
  proj <- make_dual_project()
  diff1 <- run_diff(proj$experiments$proteomics, method = "ttest",
                    analysis_type = "group", group_col = "group",
                    control_group = "control", case_group = "treatment")
  expect_error(
    run_integration(list(), method = "concordance",
                    experiments = "x", diff_bundles = list(x = diff1)),
    "omics_project"
  )
})

test_that("run_integration errors when too few experiments", {
  proj <- make_dual_project()
  diff1 <- run_diff(proj$experiments$proteomics, method = "ttest",
                    analysis_type = "group", group_col = "group",
                    control_group = "control", case_group = "treatment")
  expect_error(
    run_integration(proj, method = "concordance",
                    experiments = "proteomics",
                    diff_bundles = list(proteomics = diff1)),
    "length-2"
  )
})

test_that("run_integration errors on invalid method", {
  proj <- make_dual_project()
  expect_error(
    run_integration(proj, method = "invalid",
                    experiments = c("proteomics", "rnaseq"),
                    diff_bundles = list()),
    "should be one of"
  )
})

test_that("run_integration errors when experiment tag not found", {
  proj <- make_dual_project()
  expect_error(
    run_integration(proj, method = "concordance",
                    experiments = c("proteomics", "nonexistent"),
                    diff_bundles = list()),
    "not found in project"
  )
})

# ---- plot_integration --------------------------------------------------

test_that("plot_integration renders a ggplot for concordance result", {
  proj <- make_dual_project()
  diff1 <- run_diff(proj$experiments$proteomics, method = "ttest",
                    analysis_type = "group", group_col = "group",
                    control_group = "control", case_group = "treatment")
  diff2 <- run_diff(proj$experiments$rnaseq, method = "ttest",
                    analysis_type = "group", group_col = "group",
                    control_group = "control", case_group = "treatment")
  bundles <- list(proteomics = diff1, rnaseq = diff2)
  res <- run_integration(proj, method = "concordance",
                         experiments = c("proteomics", "rnaseq"),
                         diff_bundles = bundles)
  p <- plot_integration(res)
  expect_s3_class(p, "ggplot")
})

# ---- check_integration_result_schema -----------------------------------

test_that("check_integration_result_schema passes for valid result df", {
  proj <- make_dual_project()
  diff1 <- run_diff(proj$experiments$proteomics, method = "ttest",
                    analysis_type = "group", group_col = "group",
                    control_group = "control", case_group = "treatment")
  diff2 <- run_diff(proj$experiments$rnaseq, method = "ttest",
                    analysis_type = "group", group_col = "group",
                    control_group = "control", case_group = "treatment")
  bundles <- list(proteomics = diff1, rnaseq = diff2)
  res <- run_integration(proj, method = "concordance",
                         experiments = c("proteomics", "rnaseq"),
                         diff_bundles = bundles)
  df <- res$results$integration_df
  expect_true(check_integration_result_schema(df))
})

test_that("check_integration_result_schema errors when columns missing", {
  expect_error(check_integration_result_schema(data.frame(x = 1)),
               "Missing required")
})
