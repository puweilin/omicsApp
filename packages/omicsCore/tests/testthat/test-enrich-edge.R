# Edge-case and robustness tests for run_enrichment, enrich utilities,
# and filter_enrich_results.

make_diff_bundle_for_enrich <- function() {
  n_feat <- 30
  n_samp <- 8
  mat <- matrix(rnorm(n_feat * n_samp, mean = 10, sd = 2),
                nrow = n_feat, ncol = n_samp)
  rownames(mat) <- paste0("gene_", seq_len(n_feat))
  colnames(mat) <- paste0("sample_", seq_len(n_samp))
  meta <- data.frame(
    group = rep(c("control", "treatment"), each = n_samp / 2),
    row.names = colnames(mat)
  )
  feat <- data.frame(
    feature_id = rownames(mat),
    feature_name = paste0("Gene", seq_len(n_feat)),
    stringsAsFactors = FALSE
  )
  inp <- omics_input(mat, meta, feat, omics_type = "proteomics",
                     assay_type = "normalized_intensity")
  run_diff(inp, method = "ttest", analysis_type = "group",
           group_col = "group",
           control_group = "control", case_group = "treatment")
}

# ---- run_enrichment: basic execution -----------------------------------

test_that("run_enrichment ORA with hallmark returns an analysis_bundle", {
  skip_if_not_installed("clusterProfiler")
  bundle <- make_diff_bundle_for_enrich()
  res <- run_enrichment(bundle, type = "ora", database = "hallmark")
  expect_true(is_analysis_bundle(res))
  expect_equal(res$analysis_name, "run_enrichment")
})

test_that("run_enrichment returns standardized result columns", {
  skip_if_not_installed("clusterProfiler")
  bundle <- make_diff_bundle_for_enrich()
  res <- run_enrichment(bundle, type = "ora", database = "hallmark")
  df <- res$results$enrich_result_df
  expect_true("pathway_name" %in% colnames(df))
  expect_true("adj_p_value" %in% colnames(df))
})

test_that("run_enrichment ORA with 'up' direction works", {
  skip_if_not_installed("clusterProfiler")
  bundle <- make_diff_bundle_for_enrich()
  res <- run_enrichment(bundle, type = "ora", database = "hallmark",
                        direction = "up")
  expect_true(is_analysis_bundle(res))
})

test_that("run_enrichment ORA with 'down' direction works", {
  skip_if_not_installed("clusterProfiler")
  bundle <- make_diff_bundle_for_enrich()
  res <- run_enrichment(bundle, type = "ora", database = "hallmark",
                        direction = "down")
  expect_true(is_analysis_bundle(res))
})

test_that("run_enrichment ORA with different databases works", {
  skip_if_not_installed("clusterProfiler")
  bundle <- make_diff_bundle_for_enrich()
  for (db in c("kegg", "reactome", "go_bp")) {
    has_sets <- tryCatch(
      nrow(msigdbr::msigdbr("Homo sapiens", db)) > 0,
      error = function(e) FALSE
    )
    skip_if_not(has_sets, sprintf("No %s gene sets available", db))
    res <- run_enrichment(bundle, type = "ora", database = db)
    expect_true(is_analysis_bundle(res))
  }
})

# ---- run_enrichment: GSEA ----------------------------------------------

test_that("run_enrichment GSEA with hallmark works", {
  skip_if_not_installed("clusterProfiler")
  # Real symbols, or no set overlaps the list and GSEA fails inside a
  # tryCatch while this test still passes on the class of the bundle.
  # test-enrich-realistic.R holds the full assertions.
  bundle <- realistic_diff_bundle()
  res <- suppressWarnings(run_enrichment(bundle, type = "gsea", database = "hallmark"))
  expect_true(is_analysis_bundle(res))
  expect_gt(nrow(res$results$enrich_result_df), 0L)
})

# ---- run_enrichment: error paths ---------------------------------------

test_that("run_enrichment errors on invalid bundle input", {
  expect_error(run_enrichment(list(), type = "ora", database = "hallmark"),
               "analysis_bundle")
})

test_that("run_enrichment errors on invalid type", {
  skip_if_not_installed("clusterProfiler")
  bundle <- make_diff_bundle_for_enrich()
  expect_error(
    run_enrichment(bundle, type = "invalid", database = "hallmark"),
    "should be one of"
  )
})

test_that("run_enrichment errors on invalid database", {
  skip_if_not_installed("clusterProfiler")
  bundle <- make_diff_bundle_for_enrich()
  expect_error(
    run_enrichment(bundle, type = "ora", database = "nonexistent_db"),
    "database"
  )
})

test_that("run_enrichment errors on invalid direction", {
  skip_if_not_installed("clusterProfiler")
  bundle <- make_diff_bundle_for_enrich()
  expect_error(
    run_enrichment(bundle, type = "ora", database = "hallmark",
                   direction = "sideways"),
    "should be one of"
  )
})

# ---- filter_enrich_results ---------------------------------------------

test_that("filter_enrich_results filters by adj_p", {
  skip_if_not_installed("clusterProfiler")
  bundle <- make_diff_bundle_for_enrich()
  enrich <- run_enrichment(bundle, type = "ora", database = "hallmark")
  df <- enrich$results$enrich_result_df
  res <- filter_enrich_results(df, p_cutoff = 1, p_preference = "adjusted")
  expect_s3_class(res, "data.frame")
})

test_that("filter_enrich_results with strict cutoff returns empty", {
  skip_if_not_installed("clusterProfiler")
  bundle <- make_diff_bundle_for_enrich()
  enrich <- run_enrichment(bundle, type = "ora", database = "hallmark")
  df <- enrich$results$enrich_result_df
  res <- filter_enrich_results(df, p_cutoff = 1e-100)
  expect_equal(nrow(res), 0)
})

test_that("filter_enrich_results validates input", {
  expect_error(filter_enrich_results(list()), "enrich_df.*must be")
})

# ---- make_ranked_features ----------------------------------------------

test_that("make_ranked_features returns a named numeric vector", {
  bundle <- make_diff_bundle_for_enrich()
  df <- bundle$results$diff_result_df
  ranked <- make_ranked_features(df)
  expect_type(ranked, "double")
  expect_true(length(names(ranked)) > 0)
  expect_true(all(!is.na(ranked)))
})

test_that("make_ranked_features validates input", {
  expect_error(make_ranked_features(list()), "Missing required")
})

# ---- list_gene_sets ----------------------------------------------------

test_that("list_gene_sets returns a non-empty data.frame", {
  skip_if_not_installed("clusterProfiler")
  res <- list_gene_sets("hallmark")
  expect_s3_class(res, "data.frame")
  expect_gt(nrow(res), 0)
})

# ---- check_enrich_result_schema ----------------------------------------

test_that("check_enrich_result_schema passes for valid result df", {
  skip_if_not_installed("clusterProfiler")
  bundle <- make_diff_bundle_for_enrich()
  enrich <- run_enrichment(bundle, type = "ora", database = "hallmark")
  df <- enrich$results$enrich_result_df
  expect_true(check_enrich_result_schema(df))
})

test_that("check_enrich_result_schema errors when columns missing", {
  expect_error(check_enrich_result_schema(data.frame(x = 1)),
               "Missing required")
})
