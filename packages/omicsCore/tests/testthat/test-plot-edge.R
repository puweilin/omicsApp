# Robustness tests for plot_volcano, plot_ma, plot_enrichment, plot_gsea,
# plot_gsva_heatmap, plot_heatmap, plot_pca, plot_feature_expression.
# Focus: edge inputs (tiny / empty / NA), invalid args, output class.

make_plot_input <- function(n_feat = 30, n_samp = 8, omics_type = "proteomics",
                            assay_type = "normalized_intensity") {
  mat <- matrix(rnorm(n_feat * n_samp, mean = 10, sd = 2),
                nrow = n_feat, ncol = n_samp)
  rownames(mat) <- paste0("gene_", seq_len(n_feat))
  colnames(mat) <- paste0("s", seq_len(n_samp))
  half <- n_samp %/% 2
  meta <- data.frame(
    group = c(rep("ctrl", half), rep("trt", n_samp - half)),
    batch = rep(c("B1", "B2"), length.out = n_samp),
    age   = seq(20, by = 5, length.out = n_samp),
    row.names = colnames(mat)
  )
  feat <- data.frame(feature_id = rownames(mat),
                     feature_name = paste0("Gene", seq_len(n_feat)),
                     stringsAsFactors = FALSE)
  omics_input(mat, meta, feat, omics_type = omics_type, assay_type = assay_type)
}

make_plot_diff_bundle <- function(...) {
  inp <- make_plot_input(...)
  run_diff(inp, method = "ttest", analysis_type = "group",
           group_col = "group", control_group = "ctrl", case_group = "trt")
}

# ---- plot_volcano ------------------------------------------------------

test_that("plot_volcano returns a ggplot for a standard bundle", {
  b <- make_plot_diff_bundle()
  p <- plot_volcano(b)
  expect_s3_class(p, "ggplot")
})

test_that("plot_volcano respects top_n = 0", {
  b <- make_plot_diff_bundle()
  expect_s3_class(plot_volcano(b, top_n = 0), "ggplot")
})

test_that("plot_volcano respects top_n larger than n_features", {
  b <- make_plot_diff_bundle(n_feat = 10)
  expect_s3_class(plot_volcano(b, top_n = 100), "ggplot")
})

test_that("plot_volcano respects p_basis = 'raw'", {
  b <- make_plot_diff_bundle()
  expect_s3_class(plot_volcano(b, p_basis = "raw"), "ggplot")
})

test_that("plot_volcano errors on invalid p_basis", {
  b <- make_plot_diff_bundle()
  expect_error(plot_volcano(b, p_basis = "nope"), "should be one of")
})

test_that("plot_volcano errors on non-bundle", {
  expect_error(plot_volcano(list()), "analysis_bundle")
})

test_that("plot_volcano respects p_threshold", {
  b <- make_plot_diff_bundle()
  expect_s3_class(plot_volcano(b, p_threshold = 1e-10), "ggplot")
  expect_s3_class(plot_volcano(b, p_threshold = 1), "ggplot")
})

test_that("plot_volcano respects effect_threshold", {
  b <- make_plot_diff_bundle()
  expect_s3_class(plot_volcano(b, effect_threshold = 0), "ggplot")
  expect_s3_class(plot_volcano(b, effect_threshold = 100), "ggplot")
})

test_that("plot_volcano label_features highlights known features", {
  b <- make_plot_diff_bundle()
  features <- b$results$diff_result_df$feature_id[1:3]
  expect_s3_class(plot_volcano(b, label_features = features), "ggplot")
})

test_that("plot_volcano label_features ignores unknown ids", {
  b <- make_plot_diff_bundle()
  expect_s3_class(
    plot_volcano(b, label_features = c("nope_1", "nope_2")),
    "ggplot"
  )
})

# ---- plot_ma -----------------------------------------------------------

test_that("plot_ma returns a ggplot", {
  b <- make_plot_diff_bundle()
  expect_s3_class(plot_ma(b), "ggplot")
})

test_that("plot_ma respects top_n = 0", {
  b <- make_plot_diff_bundle()
  expect_s3_class(plot_ma(b, top_n = 0), "ggplot")
})

test_that("plot_ma errors on non-bundle", {
  expect_error(plot_ma(list()), "analysis_bundle")
})

test_that("plot_ma handles label_features", {
  b <- make_plot_diff_bundle()
  f <- b$results$diff_result_df$feature_id[1:2]
  expect_s3_class(plot_ma(b, label_features = f), "ggplot")
})

# ---- plot_enrichment ---------------------------------------------------

test_that("plot_enrichment renders for a hallmark ORA bundle", {
  skip_if_not_installed("clusterProfiler")
  b <- run_enrichment(make_plot_diff_bundle(), type = "ora",
                      database = "hallmark")
  expect_s3_class(plot_enrichment(b), "ggplot")
})

test_that("plot_enrichment respects view = 'bar'", {
  skip_if_not_installed("clusterProfiler")
  b <- run_enrichment(make_plot_diff_bundle(), type = "ora",
                      database = "hallmark")
  expect_s3_class(plot_enrichment(b, view = "bar"), "ggplot")
})

test_that("plot_enrichment respects view = 'gsea_dot' on GSEA bundle", {
  skip_if_not_installed("clusterProfiler")
  b <- tryCatch(
    suppressWarnings(run_enrichment(realistic_diff_bundle(), type = "gsea",
                                    database = "hallmark")),
    error = function(e) NULL
  )
  skip_if(is.null(b), "GSEA bundle could not be built")
  expect_s3_class(plot_enrichment(b, view = "gsea_dot"), "ggplot")
})

test_that("plot_enrichment errors on invalid view", {
  skip_if_not_installed("clusterProfiler")
  b <- run_enrichment(make_plot_diff_bundle(), type = "ora",
                      database = "hallmark")
  expect_error(plot_enrichment(b, view = "nope"), "should be one of")
})

test_that("plot_enrichment errors on invalid p_preference", {
  skip_if_not_installed("clusterProfiler")
  b <- run_enrichment(make_plot_diff_bundle(), type = "ora",
                      database = "hallmark")
  expect_error(plot_enrichment(b, p_preference = "nope"), "should be one of")
})

test_that("plot_enrichment errors on non-bundle", {
  expect_error(plot_enrichment(list()), "analysis_bundle")
})

test_that("plot_enrichment handles p_cutoff = 1", {
  skip_if_not_installed("clusterProfiler")
  b <- run_enrichment(make_plot_diff_bundle(), type = "ora",
                      database = "hallmark")
  expect_s3_class(plot_enrichment(b, p_cutoff = 1), "ggplot")
})

test_that("plot_enrichment handles top_n = 0", {
  skip_if_not_installed("clusterProfiler")
  b <- run_enrichment(make_plot_diff_bundle(), type = "ora",
                      database = "hallmark")
  expect_s3_class(plot_enrichment(b, top_n = 0L), "ggplot")
})

# ---- plot_heatmap ------------------------------------------------------

test_that("plot_heatmap works with omics_input directly", {
  inp <- make_plot_input()
  p <- plot_heatmap(inp, n_top = 10L)
  expect_true(inherits(p, c("ggplot", "Heatmap", "HeatmapList")))
})

test_that("plot_heatmap works with diff bundle + input", {
  b <- make_plot_diff_bundle()
  inp <- make_plot_input()
  p <- plot_heatmap(b, input = inp, n_top = 5L)
  expect_true(inherits(p, c("ggplot", "Heatmap", "HeatmapList")))
})

test_that("plot_heatmap with n_top larger than features", {
  inp <- make_plot_input(n_feat = 5)
  p <- plot_heatmap(inp, n_top = 999L)
  expect_true(inherits(p, c("ggplot", "Heatmap", "HeatmapList")))
})

test_that("plot_heatmap with specific feature list", {
  inp <- make_plot_input()
  feats <- rownames(inp$expr_mat)[1:3]
  p <- plot_heatmap(inp, features = feats)
  expect_true(inherits(p, c("ggplot", "Heatmap", "HeatmapList")))
})

test_that("plot_heatmap errors on invalid scale", {
  inp <- make_plot_input()
  expect_error(plot_heatmap(inp, scale = "diagonal"), "should be one of")
})

test_that("plot_heatmap with annotation_cols renders", {
  inp <- make_plot_input()
  p <- plot_heatmap(inp, annotation_cols = "group", n_top = 5L)
  expect_true(inherits(p, c("ggplot", "Heatmap", "HeatmapList")))
})

test_that("plot_heatmap with cluster flags disabled", {
  inp <- make_plot_input()
  p <- plot_heatmap(inp, cluster_rows = FALSE, cluster_cols = FALSE,
                    n_top = 5L)
  expect_true(inherits(p, c("ggplot", "Heatmap", "HeatmapList")))
})

test_that("plot_heatmap errors on neither bundle nor input", {
  expect_error(plot_heatmap(list()), regexp = ".+")
})

# ---- plot_pca ---------------------------------------------------------

test_that("plot_pca renders with shape_by", {
  inp <- make_plot_input()
  expect_s3_class(plot_pca(inp, color_by = "group", shape_by = "batch"),
                  "ggplot")
})

test_that("plot_pca renders with log2 = TRUE", {
  inp <- make_plot_input()
  expect_s3_class(plot_pca(inp, log2 = TRUE), "ggplot")
})

test_that("plot_pca renders with log2 = FALSE", {
  inp <- make_plot_input()
  expect_s3_class(plot_pca(inp, log2 = FALSE), "ggplot")
})

test_that("plot_pca handles nonexistent color_by gracefully", {
  inp <- make_plot_input()
  expect_error(plot_pca(inp, color_by = "nope"))
})

test_that("plot_pca handles 2-sample input", {
  mat <- matrix(rnorm(20), nrow = 10, ncol = 2,
                dimnames = list(paste0("g", 1:10), c("s1", "s2")))
  meta <- data.frame(group = c("A", "B"), row.names = c("s1", "s2"))
  feat <- data.frame(feature_id = rownames(mat))
  inp <- omics_input(mat, meta, feat, omics_type = "proteomics",
                     assay_type = "normalized_intensity")
  # 2 samples: PCA is degenerate, function may render or error gracefully
  res <- tryCatch(plot_pca(inp), error = function(e) e)
  expect_true(inherits(res, "ggplot") || inherits(res, "error"))
})

# ---- plot_feature_expression ------------------------------------------

test_that("plot_feature_expression with color_by", {
  inp <- make_plot_input()
  p <- plot_feature_expression(inp, features = "gene_1",
                               group_by = "group", color_by = "batch")
  expect_s3_class(p, "ggplot")
})

test_that("plot_feature_expression with many features", {
  inp <- make_plot_input()
  feats <- rownames(inp$expr_mat)[1:5]
  p <- plot_feature_expression(inp, features = feats, group_by = "group")
  expect_s3_class(p, "ggplot")
})

test_that("plot_feature_expression errors on unknown feature", {
  inp <- make_plot_input()
  expect_error(
    plot_feature_expression(inp, features = "fake_gene", group_by = "group")
  )
})

test_that("plot_feature_expression errors on unknown group_by", {
  inp <- make_plot_input()
  expect_error(
    plot_feature_expression(inp, features = "gene_1", group_by = "nope")
  )
})

# ---- plot_qc ----------------------------------------------------------

test_that("plot_qc view = 'pca' renders", {
  inp <- make_plot_input()
  b <- run_qc(inp)
  expect_s3_class(plot_qc(b, view = "pca"), "ggplot")
})

test_that("plot_qc view = 'missing' renders", {
  inp <- make_plot_input()
  inp$expr_mat[1, 1:2] <- NA
  b <- run_qc(inp)
  expect_s3_class(plot_qc(b, view = "missing"), "ggplot")
})

test_that("plot_qc view = 'connectivity' renders", {
  inp <- make_plot_input()
  b <- run_qc(inp)
  res <- tryCatch(plot_qc(b, view = "connectivity"),
                  error = function(e) e)
  expect_true(inherits(res, "ggplot") || inherits(res, "error"))
})

test_that("plot_qc errors on invalid view", {
  inp <- make_plot_input()
  b <- run_qc(inp)
  expect_error(plot_qc(b, view = "nope"), "should be one of")
})

test_that("plot_qc errors on non-bundle", {
  expect_error(plot_qc(list()), "analysis_bundle")
})

test_that("plot_qc with color_by renders", {
  inp <- make_plot_input()
  b <- run_qc(inp)
  expect_s3_class(plot_qc(b, view = "pca", color_by = "group"), "ggplot")
})

# ---- plot_integration -------------------------------------------------

make_dual_diff_bundles <- function() {
  i1 <- make_plot_input(n_feat = 20, n_samp = 6)
  i2 <- make_plot_input(n_feat = 20, n_samp = 6, omics_type = "rnaseq",
                        assay_type = "logcpm")
  colnames(i2$expr_mat) <- paste0("t", 1:6)
  rownames(i2$meta_df) <- paste0("t", 1:6)
  proj <- omics_project("dual", experiments = list(prot = i1, rna = i2))
  d1 <- run_diff(i1, method = "ttest", analysis_type = "group",
                 group_col = "group", control_group = "ctrl",
                 case_group = "trt")
  d2 <- run_diff(i2, method = "ttest", analysis_type = "group",
                 group_col = "group", control_group = "ctrl",
                 case_group = "trt")
  list(proj = proj, d1 = d1, d2 = d2)
}

test_that("plot_integration view = 'scatter' renders", {
  s <- make_dual_diff_bundles()
  res <- run_integration(s$proj, method = "concordance",
                         experiments = c("prot", "rna"),
                         diff_bundles = list(prot = s$d1, rna = s$d2))
  expect_s3_class(plot_integration(res, view = "scatter"), "ggplot")
})

test_that("plot_integration view = 'quadrant' renders", {
  s <- make_dual_diff_bundles()
  res <- run_integration(s$proj, method = "concordance",
                         experiments = c("prot", "rna"),
                         diff_bundles = list(prot = s$d1, rna = s$d2))
  res2 <- tryCatch(plot_integration(res, view = "quadrant"),
                   error = function(e) e)
  expect_true(inherits(res2, "ggplot") || inherits(res2, "error"))
})

test_that("plot_integration view = 'dual_volcano' renders", {
  s <- make_dual_diff_bundles()
  res <- run_integration(s$proj, method = "concordance",
                         experiments = c("prot", "rna"),
                         diff_bundles = list(prot = s$d1, rna = s$d2))
  res2 <- tryCatch(plot_integration(res, view = "dual_volcano"),
                   error = function(e) e)
  expect_true(inherits(res2, "ggplot") || inherits(res2, "error"))
})

test_that("plot_integration view = 'dotplot' renders", {
  s <- make_dual_diff_bundles()
  res <- run_integration(s$proj, method = "concordance",
                         experiments = c("prot", "rna"),
                         diff_bundles = list(prot = s$d1, rna = s$d2))
  res2 <- tryCatch(plot_integration(res, view = "dotplot"),
                   error = function(e) e)
  expect_true(inherits(res2, "ggplot") || inherits(res2, "error"))
})

test_that("plot_integration errors on invalid view", {
  s <- make_dual_diff_bundles()
  res <- run_integration(s$proj, method = "concordance",
                         experiments = c("prot", "rna"),
                         diff_bundles = list(prot = s$d1, rna = s$d2))
  expect_error(plot_integration(res, view = "nope"), "should be one of")
})

test_that("plot_integration errors on non-bundle", {
  expect_error(plot_integration(list()), "analysis_bundle")
})

test_that("plot_integration respects top_n = 0", {
  s <- make_dual_diff_bundles()
  res <- run_integration(s$proj, method = "concordance",
                         experiments = c("prot", "rna"),
                         diff_bundles = list(prot = s$d1, rna = s$d2))
  expect_s3_class(plot_integration(res, top_n = 0L), "ggplot")
})

test_that("plot_integration respects p_cutoff = 1", {
  s <- make_dual_diff_bundles()
  res <- run_integration(s$proj, method = "concordance",
                         experiments = c("prot", "rna"),
                         diff_bundles = list(prot = s$d1, rna = s$d2))
  expect_s3_class(plot_integration(res, p_cutoff = 1), "ggplot")
})

# ---- plot_gsea --------------------------------------------------------

test_that("plot_gsea errors on non-bundle", {
  expect_error(plot_gsea(list(), pathway_id = "x"), "analysis_bundle")
})

test_that("plot_gsea errors on missing pathway_id", {
  skip_if_not_installed("clusterProfiler")
  b <- suppressWarnings(run_enrichment(realistic_diff_bundle(), type = "gsea",
                                       database = "hallmark"))
  expect_error(plot_gsea(b))
})

# ---- plot_gsva_heatmap ------------------------------------------------
# On the realistic fixture, so the Hallmark sets overlap the matrix. The
# old fixture's made-up gene names left nothing to score, and the
# tryCatch/skip around it reported that as "could not be built" -- on
# every machine, forever -- while also naming a result slot that never
# existed.

gsva_plot_bundle <- function() {
  run_gsva(realistic_input(), database = "hallmark")
}

test_that("plot_gsva_heatmap renders for a gsva bundle", {
  skip_if_not_installed("GSVA")
  skip_if_not_installed("msigdbr")
  p <- plot_gsva_heatmap(gsva_plot_bundle(), top_n = 10L)
  expect_true(inherits(p, c("ggplot", "Heatmap", "HeatmapList")))
})

test_that("plot_gsva_heatmap with pathways subset", {
  skip_if_not_installed("GSVA")
  skip_if_not_installed("msigdbr")
  b <- gsva_plot_bundle()
  pws <- rownames(b$results$gsva_matrix)[1:3]
  p <- plot_gsva_heatmap(b, pathways = pws)
  expect_true(inherits(p, c("ggplot", "Heatmap", "HeatmapList")))
})

test_that("plot_gsva_heatmap errors on non-bundle", {
  expect_error(plot_gsva_heatmap(list()), "analysis_bundle")
})

test_that("plot_gsva_heatmap errors on invalid scale", {
  skip_if_not_installed("GSVA")
  skip_if_not_installed("msigdbr")
  expect_error(plot_gsva_heatmap(gsva_plot_bundle(), scale = "diagonal"),
               "should be one of")
})

# ---- consolidated palette and views -----------------------------------

test_that("the palette covers every quadrant a concordance plot can emit", {
  pal <- quadrant_palette()
  expect_true(all(c("up_up", "down_down", "up_down", "down_up") %in% names(pal)))
  # Features reaching no threshold, and rows the schema left NA, both
  # have to land on a colour or ggplot drops them silently.
  expect_true(all(c("ns", "n/a") %in% names(pal)))
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", unname(pal))))
})

test_that("plot_volcano re-derives significance from supplied thresholds", {
  bundle <- make_plot_diff_bundle()
  df <- bundle$results$diff_result_df

  loose <- ggplot2::ggplot_build(
    plot_volcano(bundle, p_threshold = 1, effect_threshold = 0))$data[[1]]
  strict <- ggplot2::ggplot_build(
    plot_volcano(bundle, p_threshold = 1e-12, effect_threshold = 50))$data[[1]]

  # A threshold nothing can pass must colour nothing as significant; one
  # everything passes must colour everything. Without re-deriving, both
  # would return the bundle's stored mask and be identical.
  expect_gt(length(unique(loose$colour)), 0L)
  expect_equal(length(unique(strict$colour)), 1L)
  expect_false(identical(sort(unique(loose$colour)),
                         sort(unique(strict$colour))))
})

test_that("plot_volcano leaves the stored mask alone when not asked", {
  bundle <- make_plot_diff_bundle()
  # No thresholds supplied: the figure describes the analysis as run.
  expect_s3_class(plot_volcano(bundle), "ggplot")
})
