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

# ---- slice 2E fixtures ------------------------------------------------

test_that("example_enrich_table() matches the enrich_result_df schema", {
  df <- example_enrich_table()
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 15L)
  required <- c("database", "result_type", "comparison", "pathway_id",
                "pathway_name", "effect", "effect_type", "direction",
                "p_value", "adj_p_value", "q_value", "gene_set_size",
                "overlap_size", "overlap_features", "leading_features",
                "source_label")
  expect_true(all(required %in% colnames(df)))
  expect_equal(unique(df$database), "hallmark")
  expect_equal(unique(df$result_type), "gsea")
  expect_equal(unique(df$effect_type), "nes")
  expect_setequal(unique(df$direction), c("up", "down"))
  # All rows clear a strict significance bar so the dotplot has signal.
  expect_true(all(df$adj_p_value < 0.01))
  # NES sign agrees with the direction column.
  expect_true(all(ifelse(df$effect >= 0, "up", "down") == df$direction))
})

test_that("example_enrich_table() is deterministic", {
  a <- example_enrich_table()
  b <- example_enrich_table()
  expect_identical(a, b)
})

test_that("example_integration_tables() returns the expected list", {
  tt <- example_integration_tables()
  expect_named(tt, c("concordance_df", "active_pathways_df"))

  conc <- tt$concordance_df
  expect_s3_class(conc, "data.frame")
  expect_equal(nrow(conc), 60L)
  expect_true(all(c("feature_id", "feature_symbol", "effect_a", "effect_b",
                    "effect_diff", "p_value", "adj_p_value", "quadrant")
                  %in% colnames(conc)))
  expect_true(all(conc$quadrant %in%
                  c("up_up", "down_down", "up_down", "down_up", "ns")))
  # The seeded signal in features 1-15 / 16-30 must surface as the
  # dominant concordant quadrants — the dual-volcano / scatter rely on
  # this for visual interpretability.
  expect_gte(sum(conc$quadrant == "up_up"),     5L)
  expect_gte(sum(conc$quadrant == "down_down"), 5L)

  ap <- tt$active_pathways_df
  expect_s3_class(ap, "data.frame")
  expect_equal(nrow(ap), 6L)
  expect_true(all(c("pathway_id", "pathway_name",
                    "p_a", "p_b", "p_combined", "direction")
                  %in% colnames(ap)))
  expect_true(all(ap$p_combined > 0 & ap$p_combined < 1))
  # Brown's combined p must be at least as significant as the smaller
  # per-omics p for each row (a basic sanity check on the fixture).
  expect_true(all(ap$p_combined <= pmin(ap$p_a, ap$p_b)))
})

test_that("example_integration_tables() is deterministic", {
  a <- example_integration_tables()
  b <- example_integration_tables()
  expect_identical(a, b)
})
