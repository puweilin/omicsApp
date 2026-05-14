# Tests for the built-in example fixtures used to populate the UI
# views before any user upload exists.

test_that("example_input('proteomics') returns a valid omics_input", {
  x <- example_input("proteomics")
  expect_true(omicsCore::is_omics_input(x))
  expect_identical(x$omics_type, "proteomics")
  expect_identical(x$assay_type, "normalized_intensity")
  expect_equal(dim(x$expr_mat), c(50L, 12L))
  expect_named(x$meta_df, c("group", "age", "sex", "donor_id", "batch"))
  expect_true(all(c("feature_id", "feature_symbol", "description")
                  %in% colnames(x$feature_df)))
  expect_setequal(unique(x$meta_df$group), c("G1", "G2"))
})

test_that("example_input('rnaseq') returns a valid omics_input", {
  x <- example_input("rnaseq")
  expect_true(omicsCore::is_omics_input(x))
  expect_identical(x$omics_type, "rnaseq")
  expect_identical(x$assay_type, "raw_count")
  expect_equal(dim(x$expr_mat), c(60L, 12L))
  # ENSG-style feature IDs.
  expect_true(all(grepl("^ENSG[0-9]+$", rownames(x$expr_mat))))
})

test_that("example_input() builds are deterministic", {
  a <- example_input("proteomics")
  b <- example_input("proteomics")
  expect_identical(a$expr_mat, b$expr_mat)
  expect_identical(a$meta_df,  b$meta_df)
})

test_that("example_diff_bundle() returns a populated analysis_bundle", {
  b <- example_diff_bundle()
  expect_true(omicsCore::is_analysis_bundle(b))
  expect_identical(b$analysis_name, "run_diff")
  expect_identical(b$params$method, "limma")
  expect_true("diff_result_df" %in% names(b$results))
  expect_s3_class(b$results$diff_result_df, "data.frame")
  expect_gt(nrow(b$results$diff_result_df), 0L)
})

test_that("example_diff_bundle() top-10 |effect| features are the seeded ones", {
  # The fixture injects up/down effects on features P001-P010. If limma
  # is healthy and the signal stays strong, every one of the top 10 by
  # |effect| should be a seeded feature.
  df <- example_diff_bundle()$results$diff_result_df
  top10 <- df[order(-abs(df$effect)), ][1:10, ]
  seeded <- sprintf("P%03d", 1:10)
  expect_setequal(top10$feature_id, seeded)
  # All seeded features should clear adj_p < 0.1 with this signal size.
  expect_true(all(top10$adj_p_value < 0.1))
})

test_that("example_project() carries both proteomics and rnaseq experiments", {
  p <- example_project()
  expect_true(omicsCore::is_omics_project(p))
  expect_named(p$experiments, c("proteomics", "rnaseq"))
  expect_identical(p$experiments$proteomics$omics_type, "proteomics")
  expect_identical(p$experiments$rnaseq$omics_type,     "rnaseq")
})

test_that("example_qc_bundle() returns a populated QC analysis_bundle", {
  b <- example_qc_bundle()
  expect_true(omicsCore::is_analysis_bundle(b))
  expect_identical(b$analysis_name, "run_qc")
  # The injected ~5% NAs must surface in the per-feature missing-rate
  # table so the missingness panel has something to draw.
  miss <- b$results$qc_summary$missingness$feature_metrics
  expect_s3_class(miss, "data.frame")
  expect_true(any(miss$missing_rate > 0))
})

test_that("example_qc_bundle() is deterministic", {
  a <- example_qc_bundle()
  b <- example_qc_bundle()
  expect_identical(
    a$results$qc_summary$missingness$feature_metrics$missing_rate,
    b$results$qc_summary$missingness$feature_metrics$missing_rate
  )
})
