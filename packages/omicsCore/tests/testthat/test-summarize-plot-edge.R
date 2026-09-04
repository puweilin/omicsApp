# Tests for summarize_omics, plot_heatmap, and related functions.

make_test_input <- function() {
  mat <- matrix(rnorm(50, mean = 10, sd = 2), nrow = 10, ncol = 5)
  rownames(mat) <- paste0("gene_", 1:10)
  colnames(mat) <- paste0("sample_", 1:5)
  meta <- data.frame(
    group = c("A", "A", "B", "B", "B"),
    row.names = colnames(mat)
  )
  feat <- data.frame(
    feature_id = rownames(mat),
    feature_name = paste0("Gene", 1:10),
    stringsAsFactors = FALSE
  )
  omics_input(mat, meta, feat, omics_type = "proteomics", assay_type = "normalized_intensity")
}

# ---- summarize_omics ---------------------------------------------------

test_that("summarize_omics returns expected structure for omics_input", {
  inp <- make_test_input()
  res <- summarize_omics(inp)
  expect_type(res, "list")
  expect_true("feature_count" %in% names(res) || "n_features" %in% names(res) ||
              "nrow" %in% names(res))
})

test_that("summarize_omics returns expected structure for omics_project", {
  inp <- make_test_input()
  proj <- omics_project("test", experiments = list(prot = inp))
  res <- summarize_omics(proj)
  expect_type(res, "list")
})

test_that("summarize_omics errors on invalid input", {
  expect_error(summarize_omics(list()), "omics_input.*omics_project")
})

test_that("summarize_omics with omics_input reports sample count", {
  inp <- make_test_input()
  res <- summarize_omics(inp)
  # should mention the number of samples somewhere
  flat <- paste(unlist(res), collapse = " ")
  expect_true(grepl("5", flat))
})

# ---- plot_heatmap ------------------------------------------------------

test_that("plot_heatmap renders a heatmap object", {
  inp <- make_test_input()
  p <- plot_heatmap(inp)
  expect_s4_class(p, "Heatmap")
})

test_that("plot_heatmap errors on invalid input", {
  expect_error(plot_heatmap(list()), "omics_input")
})

# ---- plot_pca ----------------------------------------------------------

test_that("plot_pca renders a ggplot", {
  inp <- make_test_input()
  p <- plot_pca(inp)
  expect_s3_class(p, "ggplot")
})

test_that("plot_pca uses color_by when provided", {
  inp <- make_test_input()
  p <- plot_pca(inp, color_by = "group")
  expect_s3_class(p, "ggplot")
})

test_that("plot_pca errors on invalid input", {
  expect_error(plot_pca(list()), "omics_input")
})

# ---- plot_qc -----------------------------------------------------------

test_that("plot_qc renders a ggplot from a qc bundle", {
  inp <- make_test_input()
  bundle <- run_qc(inp)
  p <- plot_qc(bundle)
  expect_s3_class(p, "ggplot")
})

test_that("plot_qc errors on invalid input", {
  expect_error(plot_qc(list()), "analysis_bundle")
})

# ---- plot_feature_expression -------------------------------------------

test_that("plot_feature_expression renders for a valid feature", {
  inp <- make_test_input()
  p <- plot_feature_expression(inp, features = "gene_1", group_by = "group")
  expect_s3_class(p, "ggplot")
})

test_that("plot_feature_expression works with multiple features", {
  inp <- make_test_input()
  p <- plot_feature_expression(inp, features = c("gene_1", "gene_2"),
                               group_by = "group")
  expect_s3_class(p, "ggplot")
})

test_that("plot_feature_expression errors on invalid input", {
  expect_error(
    plot_feature_expression(list(), features = "x", group_by = "group"),
    "omics_input"
  )
})

# ---- theme_omicsCore ---------------------------------------------------

test_that("theme_omicsCore returns a ggplot2 theme", {
  th <- theme_omicsCore()
  expect_s3_class(th, "theme")
})
