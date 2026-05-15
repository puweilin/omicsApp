# Extended edge-case tests for example data fixtures.

# ---- example_input -----------------------------------------------------

test_that("example_input proteomics creates valid omics_input", {
  inp <- example_input("proteomics")
  expect_s3_class(inp, "omics_input")
  expect_equal(inp$omics_type, "proteomics")
  expect_equal(ncol(inp$expr_mat), 12)
  expect_equal(nrow(inp$expr_mat), 50)
  expect_true(inherits(inp$meta_df, "data.frame"))
  expect_true(inherits(inp$feature_df, "data.frame"))
})

test_that("example_input rnaseq creates valid omics_input", {
  inp <- example_input("rnaseq")
  expect_s3_class(inp, "omics_input")
  expect_equal(inp$omics_type, "rnaseq")
  expect_equal(ncol(inp$expr_mat), 12)
  expect_equal(nrow(inp$expr_mat), 60)
})

test_that("example_input has all metadata columns", {
  inp <- example_input("proteomics")
  expect_true(all(c("group", "age", "sex", "donor_id", "batch") %in%
                  colnames(inp$meta_df)))
})

test_that("example_input proteomics has expression values", {
  inp <- example_input("proteomics")
  expect_true(is.numeric(inp$expr_mat))
  expect_true(all(is.finite(inp$expr_mat) | is.na(inp$expr_mat)))
})

test_that("example_input is deterministic with seed", {
  a1 <- example_input("proteomics")$expr_mat[1, 1]
  a2 <- example_input("proteomics")$expr_mat[1, 1]
  expect_equal(a1, a2)
})

test_that("example_input validates omics_type", {
  expect_error(example_input("invalid"), "should be one of")
})

# ---- example_diff_bundle -----------------------------------------------

test_that("example_diff_bundle returns an analysis_bundle", {
  bundle <- example_diff_bundle()
  expect_true(is_analysis_bundle(bundle))
  expect_equal(bundle$analysis_name, "run_diff")
  df <- bundle$results$diff_result_df
  expect_true("feature_id" %in% colnames(df))
  expect_true("effect" %in% colnames(df))
  expect_true("adj_p_value" %in% colnames(df))
})

test_that("example_diff_bundle has some significant features", {
  bundle <- example_diff_bundle()
  df <- bundle$results$diff_result_df
  sig <- df[df$adj_p_value < 0.05 & abs(df$effect) > 0, ]
  expect_gt(nrow(sig), 0)
})

test_that("example_diff_bundle is deterministic", {
  sig1 <- sum(example_diff_bundle()$results$diff_result_df$adj_p_value < 0.05)
  sig2 <- sum(example_diff_bundle()$results$diff_result_df$adj_p_value < 0.05)
  expect_equal(sig1, sig2)
})

# ---- example_enrich_table ----------------------------------------------

test_that("example_enrich_table returns a data.frame", {
  df <- example_enrich_table()
  expect_s3_class(df, "data.frame")
  expect_gt(nrow(df), 0)
})

test_that("example_enrich_table has required columns", {
  df <- example_enrich_table()
  expect_true("pathway_name" %in% colnames(df))
  expect_true("adj_p_value" %in% colnames(df))
  expect_true("effect" %in% colnames(df))
})

test_that("example_enrich_table is deterministic", {
  n1 <- nrow(example_enrich_table())
  n2 <- nrow(example_enrich_table())
  expect_equal(n1, n2)
})

# ---- example_integration_tables ----------------------------------------

test_that("example_integration_tables returns expected named list", {
  tables <- example_integration_tables()
  expect_type(tables, "list")
  expect_true("concordance_df" %in% names(tables))
  expect_true("active_pathways_df" %in% names(tables))
})

test_that("example_integration_tables concordance_df has paired features", {
  tables <- example_integration_tables()
  df <- tables$concordance_df
  expect_s3_class(df, "data.frame")
  expect_gt(nrow(df), 0)
  expect_true("effect" %in% colnames(df) || "effect_diff" %in% colnames(df))
})

test_that("example_integration_tables active_pathways_df has pathways", {
  tables <- example_integration_tables()
  ap <- tables$active_pathways_df
  expect_s3_class(ap, "data.frame")
  expect_gt(nrow(ap), 0)
  expect_true("pathway_name" %in% colnames(ap))
})

test_that("example_integration_tables is deterministic", {
  t1 <- example_integration_tables()
  t2 <- example_integration_tables()
  expect_equal(nrow(t1$concordance_df), nrow(t2$concordance_df))
  expect_equal(nrow(t1$active_pathways_df), nrow(t2$active_pathways_df))
})
