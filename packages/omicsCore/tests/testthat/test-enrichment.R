# Tests for the enrichment slice (ORA, GSEA, GSVA, plots). Heavy
# Bioconductor backends (clusterProfiler, msigdbr, GSVA, ComplexHeatmap) are
# gated by skip_if_not_installed() so the suite stays green on a fresh
# install.

make_proteo_enrich <- function(n_features = 80L, n_per_group = 6L) {
  set.seed(2025)
  n_samples <- 2L * n_per_group
  feat_ids <- paste0("f", seq_len(n_features))
  samp_ids <- paste0("s", seq_len(n_samples))
  groups <- rep(c("ctrl", "case"), each = n_per_group)
  expr <- matrix(
    rnorm(n_features * n_samples, mean = 5, sd = 1.0),
    nrow = n_features,
    dimnames = list(feat_ids, samp_ids)
  )
  expr[1:10, groups == "case"] <- expr[1:10, groups == "case"] + 3
  meta <- data.frame(
    group = groups,
    row.names = samp_ids,
    stringsAsFactors = FALSE
  )
  symbols <- c(
    "TP53", "EGFR", "MYC", "AKT1", "BRCA1", "BRCA2", "PTEN", "RB1", "KRAS",
    "PIK3CA", "STAT3", "JUN", "FOS", "VEGFA", "TNF", "IL6", "CDKN1A",
    "MAPK1", "MAPK3", "BAX", "BCL2", "CASP3", "CASP8", "ATM", "ATR",
    "CHEK1", "CHEK2", "MDM2", "RAD51", "XRCC1"
  )
  symbols <- c(symbols, paste0("GENE", seq_len(max(0L, n_features - length(symbols)))))
  feat <- data.frame(
    feature_id = feat_ids,
    feature_symbol = symbols[seq_len(n_features)],
    row.names = feat_ids,
    stringsAsFactors = FALSE
  )
  omics_input(expr, meta, feat, omics_type = "proteomics",
              assay_type = "normalized_intensity")
}

build_test_diff_bundle <- function() {
  x <- make_proteo_enrich()
  run_diff(x, method = "ttest", analysis_type = "group",
           group_col = "group", control_group = "ctrl", case_group = "case")
}

# ---- list_gene_sets ----------------------------------------------------

test_that("list_gene_sets returns a non-empty tibble for hallmark", {
  skip_if_not_installed("msigdbr")
  skip_if_not_installed("clusterProfiler")
  out <- list_gene_sets("hallmark")
  expect_s3_class(out, "tbl_df")
  expect_true(all(c("database", "pathway_id", "pathway_name", "gene_symbol") %in%
                    colnames(out)))
  expect_gt(nrow(out), 100L)
  expect_equal(unique(out$database), "hallmark")
})

test_that("list_gene_sets accepts organism aliases", {
  skip_if_not_installed("msigdbr")
  skip_if_not_installed("clusterProfiler")
  a <- list_gene_sets("hallmark", organism = "Hs")
  b <- list_gene_sets("hallmark", organism = "Homo sapiens")
  expect_equal(nrow(a), nrow(b))
})

test_that("list_gene_sets rejects unknown databases / organisms", {
  expect_error(list_gene_sets("not_a_db"), "Unsupported enrichment database")
  expect_error(list_gene_sets("hallmark", organism = "Xx"), "Unsupported organism")
})

# ---- run_enrichment: ORA ----------------------------------------------

test_that("run_enrichment ORA returns a standardized bundle", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("msigdbr")
  b <- build_test_diff_bundle()
  enr <- run_enrichment(b, type = "ora", database = "hallmark",
                        direction = "both", p_cutoff = 0.1)
  expect_true(is_analysis_bundle(enr))
  expect_equal(enr$analysis_name, "run_enrichment")
  expect_named(enr$results, c("enrich_result_df", "enrich_object"))
  expect_s3_class(enr$results$enrich_result_df, "data.frame")
  expect_true(all(ENRICH_RESULT_REQUIRED_COLS %in% colnames(enr$results$enrich_result_df)))
  res_types <- unique(enr$results$enrich_result_df$result_type)
  expect_true(length(res_types) == 0L || identical(res_types, "ora"))
  expect_equal(enr$params$type, "ora")
})

test_that("run_enrichment ORA respects direction = 'up'", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("msigdbr")
  b <- build_test_diff_bundle()
  enr <- run_enrichment(b, type = "ora", database = "hallmark",
                        direction = "up", p_cutoff = 0.1)
  res <- enr$results$enrich_result_df
  if (nrow(res) > 0L) {
    expect_true(all(res$direction == "up"))
  }
  expect_true(all(names(enr$results$enrich_object) %in%
                    c("up__hallmark")))
})

# ---- run_enrichment: GSEA ---------------------------------------------

test_that("run_enrichment GSEA returns a standardized bundle", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("msigdbr")
  b <- build_test_diff_bundle()
  enr <- run_enrichment(b, type = "gsea", database = "hallmark",
                        p_cutoff = 1, min_size = 5L)
  expect_true(is_analysis_bundle(enr))
  res <- enr$results$enrich_result_df
  expect_true(all(ENRICH_RESULT_REQUIRED_COLS %in% colnames(res)))
  if (nrow(res) > 0L) {
    expect_equal(unique(res$result_type), "gsea")
    # NA NES (zero-variance pathway) maps to NA direction; up/down are
    # the only other allowed values.
    expect_true(all(res$direction %in% c("up", "down") |
                    is.na(res$direction)))
  }
})

# ---- filter_enrich_results --------------------------------------------

test_that("filter_enrich_results respects p_cutoff and direction", {
  df <- data.frame(
    database = "hallmark",
    result_type = "gsea",
    comparison = "x",
    pathway_id = paste0("P", 1:6),
    pathway_name = paste0("P", 1:6),
    effect = c(1, 2, -1, -2, 0.5, -0.5),
    effect_type = "nes",
    direction = c("up", "up", "down", "down", "up", "down"),
    p_value = c(0.01, 0.5, 0.001, 0.2, 0.04, 0.07),
    adj_p_value = c(0.02, 0.6, 0.002, 0.3, 0.08, 0.1),
    q_value = NA_real_,
    gene_set_size = 30,
    overlap_size = 10,
    overlap_features = NA_character_,
    leading_features = NA_character_,
    source_label = "gsea_hallmark",
    stringsAsFactors = FALSE
  )
  up <- filter_enrich_results(df, p_cutoff = 0.05, p_preference = "raw",
                              direction = "up")
  expect_true(all(up$direction == "up"))
  expect_true(all(up$p_value < 0.05))
  expect_lte(nrow(up), 2L)

  adj <- filter_enrich_results(df, p_cutoff = 0.05, p_preference = "adjusted")
  expect_true(all(adj$adj_p_value < 0.05))
})

test_that("filter_enrich_results errors on missing p column", {
  df <- data.frame(foo = 1, bar = 2)
  expect_error(filter_enrich_results(df, p_preference = "adjusted"),
               "No column found")
})

# ---- run_gsva ----------------------------------------------------------

test_that("run_gsva returns a pathways-by-samples matrix", {
  skip_if_not_installed("GSVA")
  skip_if_not_installed("msigdbr")
  x <- make_proteo_enrich()
  g <- run_gsva(x, database = "hallmark", method = "gsva",
                min_size = 5L, max_size = 500L)
  expect_true(is_analysis_bundle(g))
  expect_equal(g$analysis_name, "run_gsva")
  mat <- g$results$gsva_matrix
  expect_true(is.matrix(mat))
  expect_equal(ncol(mat), ncol(x$expr_mat))
  expect_gt(nrow(mat), 0L)
})

test_that("run_gsva accepts a user-supplied gene_sets list", {
  skip_if_not_installed("GSVA")
  x <- make_proteo_enrich()
  gs <- list(
    set_a = c("TP53", "EGFR", "MYC", "AKT1", "BRCA1"),
    set_b = c("VEGFA", "TNF", "IL6", "JUN", "FOS")
  )
  g <- run_gsva(x, gene_sets = gs, method = "gsva",
                min_size = 2L, max_size = 500L)
  expect_equal(rownames(g$results$gsva_matrix), names(gs))
})

# ---- plots -------------------------------------------------------------

test_that("plot_enrichment returns a ggplot from an ORA bundle", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("msigdbr")
  b <- build_test_diff_bundle()
  enr <- run_enrichment(b, type = "ora", database = "hallmark",
                        direction = "both", p_cutoff = 0.5)
  p <- plot_enrichment(enr, top_n = 5L, view = "dot")
  expect_s3_class(p, "ggplot")
  p2 <- plot_enrichment(enr, top_n = 5L, view = "bar")
  expect_s3_class(p2, "ggplot")
})

test_that("plot_gsva_heatmap returns an object for a GSVA bundle", {
  skip_if_not_installed("GSVA")
  skip_if_not_installed("msigdbr")
  x <- make_proteo_enrich()
  g <- run_gsva(x, database = "hallmark", method = "gsva",
                min_size = 5L, max_size = 500L)
  out <- plot_gsva_heatmap(g, top_n = 5L)
  # ComplexHeatmap path returns a Heatmap; ggplot fallback returns ggplot.
  expect_true(methods::is(out, "Heatmap") || inherits(out, "ggplot"))
})
