# Tests for plot_heatmap() + summarize_omics() (slice 1H).

make_heatmap_input <- function(n_features = 30L, n_per_group = 4L) {
  set.seed(2029)
  n_samples <- 2L * n_per_group
  expr <- matrix(
    rnorm(n_features * n_samples, mean = 5, sd = 1.0),
    nrow = n_features,
    dimnames = list(paste0("f", seq_len(n_features)),
                    paste0("s", seq_len(n_samples)))
  )
  # Make the top 5 features highly variable so they win the variance pick.
  expr[1:5, ] <- expr[1:5, ] * 5
  meta <- data.frame(group = rep(c("ctrl", "case"), each = n_per_group),
                     row.names = paste0("s", seq_len(n_samples)),
                     stringsAsFactors = FALSE)
  feat <- data.frame(feature_id = paste0("f", seq_len(n_features)),
                     feature_symbol = paste0("GENE", seq_len(n_features)),
                     row.names = paste0("f", seq_len(n_features)),
                     stringsAsFactors = FALSE)
  omics_input(expr, meta, feat, omics_type = "proteomics",
              assay_type = "normalized_intensity")
}

# ---- plot_heatmap ----------------------------------------------------

test_that("plot_heatmap accepts an omics_input directly", {
  x <- make_heatmap_input()
  out <- plot_heatmap(x, n_top = 10)
  expect_true(inherits(out, c("Heatmap", "HeatmapList", "ggplot")))
})

test_that("plot_heatmap dispatches on a diff bundle", {
  x <- make_heatmap_input()
  b <- run_diff(x, method = "ttest", analysis_type = "group",
                group_col = "group",
                control_group = "ctrl", case_group = "case")
  out <- plot_heatmap(b, input = x, n_top = 8)
  expect_true(inherits(out, c("Heatmap", "HeatmapList", "ggplot")))
})

test_that("plot_heatmap rejects a diff bundle without an input", {
  x <- make_heatmap_input()
  b <- run_diff(x, method = "ttest", analysis_type = "group",
                group_col = "group",
                control_group = "ctrl", case_group = "case")
  expect_error(plot_heatmap(b), "omics_input")
})

test_that("plot_heatmap rejects other objects", {
  expect_error(plot_heatmap(list()), "omics_input.*analysis_bundle")
})

test_that("plot_heatmap honors a forced feature list", {
  x <- make_heatmap_input()
  out <- plot_heatmap(x, features = c("f1", "f2", "f3"))
  # Hard to inspect ComplexHeatmap rows portably, so just check no error.
  expect_true(inherits(out, c("Heatmap", "HeatmapList", "ggplot")))
})

# ---- summarize_omics --------------------------------------------------

test_that("summarize_omics returns one row per omics_input", {
  x <- make_heatmap_input()
  s <- summarize_omics(x)
  expect_s3_class(s, "tbl_df")
  expect_equal(nrow(s), 1L)
  expect_identical(s$omics_type, "proteomics")
  expect_identical(s$assay_type, "normalized_intensity")
  expect_equal(s$n_samples, 8L)
  expect_equal(s$n_features, 30L)
  expect_equal(s$n_missing, 0L)
})

test_that("summarize_omics counts NAs correctly", {
  x <- make_heatmap_input()
  x$expr_mat[1, 1] <- NA_real_
  x$expr_mat[2, 2] <- NA_real_
  s <- summarize_omics(x)
  expect_equal(s$n_missing, 2L)
  expect_equal(round(s$missing_pct, 3),
               round(100 * 2 / (30 * 8), 3))
})

test_that("summarize_omics handles a multi-experiment project", {
  a <- make_heatmap_input()
  b <- make_heatmap_input(n_features = 20L)
  rownames(b$expr_mat) <- paste0("g", seq_len(nrow(b$expr_mat)))
  rownames(b$feature_df) <- paste0("g", seq_len(nrow(b$feature_df)))
  b$feature_df$feature_id <- paste0("g", seq_len(nrow(b$feature_df)))
  b$omics_type <- "rnaseq"
  p <- omics_project("demo",
                     experiments = list(proteo = a, rna = b))
  s <- summarize_omics(p)
  expect_equal(nrow(s), 2L)
  expect_setequal(s$tag, c("proteo", "rna"))
  expect_setequal(s$omics_type, c("proteomics", "rnaseq"))
})

test_that("summarize_omics returns empty rows for an empty project", {
  p <- omics_project("empty")
  s <- summarize_omics(p)
  expect_s3_class(s, "tbl_df")
  expect_equal(nrow(s), 0L)
  expect_true("tag" %in% colnames(s))
})

test_that("summarize_omics rejects other inputs", {
  expect_error(summarize_omics(list()), "omics_input.*omics_project")
})
