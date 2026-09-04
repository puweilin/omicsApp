# Test fixtures: build small but realistic omics_input objects.
# `make_proteo()` injects log2-scale signal on the first 5 features of the case
# group; `make_rna()` triples raw counts on the first 5 features.

make_proteo <- function(n_features = 30, n_per_group = 5,
                        omics_type = "proteomics") {
  set.seed(2025)
  n_samples <- 2 * n_per_group
  feat_ids <- paste0("g", seq_len(n_features))
  samp_ids <- paste0("s", seq_len(n_samples))
  groups <- rep(c("ctrl", "case"), each = n_per_group)
  expr <- matrix(
    rnorm(n_features * n_samples, mean = 5, sd = 1.0),
    nrow = n_features,
    dimnames = list(feat_ids, samp_ids)
  )
  expr[1:5, groups == "case"] <- expr[1:5, groups == "case"] + 3
  meta <- data.frame(
    group = groups,
    age = c(20, 25, 30, 35, 40, 22, 27, 32, 37, 42),
    sex = rep(c("M", "F"), length.out = n_samples),
    row.names = samp_ids,
    stringsAsFactors = FALSE
  )
  feat <- data.frame(
    feature_id = feat_ids,
    feature_symbol = paste0("S", seq_len(n_features)),
    row.names = feat_ids,
    stringsAsFactors = FALSE
  )
  omics_input(expr, meta, feat, omics_type = omics_type,
              assay_type = "normalized_intensity")
}

make_rna <- function(n_features = 30, n_per_group = 5) {
  set.seed(2025)
  n_samples <- 2 * n_per_group
  feat_ids <- paste0("g", seq_len(n_features))
  samp_ids <- paste0("s", seq_len(n_samples))
  groups <- rep(c("ctrl", "case"), each = n_per_group)
  expr <- matrix(
    rpois(n_features * n_samples, lambda = 100),
    nrow = n_features,
    dimnames = list(feat_ids, samp_ids)
  )
  expr[1:5, groups == "case"] <- expr[1:5, groups == "case"] * 3L
  meta <- data.frame(
    group = groups,
    age = c(20, 25, 30, 35, 40, 22, 27, 32, 37, 42),
    row.names = samp_ids,
    stringsAsFactors = FALSE
  )
  feat <- data.frame(
    feature_id = feat_ids,
    feature_symbol = paste0("ENS", seq_len(n_features)),
    row.names = feat_ids,
    stringsAsFactors = FALSE
  )
  omics_input(expr, meta, feat, omics_type = "rnaseq", assay_type = "raw_count")
}

# ---- run_diff: t-test backend -------------------------------------------

test_that("run_diff t-test returns a well-formed analysis_bundle", {
  x <- make_proteo()
  b <- run_diff(x, method = "ttest", analysis_type = "group",
                group_col = "group", control_group = "ctrl", case_group = "case")
  expect_true(is_analysis_bundle(b))
  expect_equal(b$analysis_name, "run_diff")
  expect_named(b$results, c("diff_result_df", "diff_raw_df", "diff_object"))
  expect_s3_class(b$results$diff_result_df, "data.frame")
  expect_true(all(DIFF_RESULT_REQUIRED_COLS %in% colnames(b$results$diff_result_df)))
  expect_equal(b$params$method, "ttest")
  expect_equal(b$params$comparison, "case_vs_ctrl")
  expect_null(b$results$diff_object)
})

test_that("run_diff t-test detects the injected signal", {
  x <- make_proteo()
  b <- run_diff(x, method = "ttest", analysis_type = "group",
                group_col = "group", control_group = "ctrl", case_group = "case")
  res <- b$results$diff_result_df
  hits <- res$feature_id[res$adj_p_value < 0.05]
  # First 5 features carry the signal; expect most of them to be flagged.
  expect_gte(length(intersect(hits, paste0("g", 1:5))), 4L)
})

# ---- run_diff: lm backend (group + continuous) --------------------------

test_that("run_diff lm group returns standardized schema", {
  x <- make_proteo()
  b <- run_diff(x, method = "lm", analysis_type = "group",
                group_col = "group", control_group = "ctrl", case_group = "case")
  res <- b$results$diff_result_df
  expect_true(all(DIFF_RESULT_REQUIRED_COLS %in% colnames(res)))
  expect_equal(unique(res$method), "lm")
  expect_equal(unique(res$effect_type), "beta")
  expect_true(all(res$direction %in% c("up", "down", "ns")))
})

test_that("run_diff lm group accepts covariates", {
  x <- make_proteo()
  b <- run_diff(x, method = "lm", analysis_type = "group",
                group_col = "group", control_group = "ctrl", case_group = "case",
                covariates = "sex")
  expect_equal(b$params$covariates, "sex")
  expect_true(nrow(b$results$diff_result_df) == nrow(x$expr_mat))
})

test_that("run_diff lm continuous returns spearman correlation as effect", {
  x <- make_proteo()
  b <- run_diff_continuous(x, method = "lm", continuous_col = "age")
  res <- b$results$diff_result_df
  expect_equal(unique(res$effect_type), "correlation")
  expect_equal(unique(res$analysis_type), "continuous_linear")
  expect_true(all(res$direction %in% c("positive", "negative", "ns")))
})

# ---- run_diff: auto dispatch + fallback ---------------------------------

test_that("run_diff auto picks limma for proteomics when installed", {
  skip_if_not_installed("limma")
  x <- make_proteo()
  b <- run_diff(x, method = "auto", analysis_type = "group",
                group_col = "group", control_group = "ctrl", case_group = "case")
  expect_equal(b$params$method, "limma")
})

# ---- run_diff: argument validation --------------------------------------

test_that("run_diff requires group_col / control_group / case_group for group analysis", {
  x <- make_proteo()
  expect_error(
    run_diff(x, method = "ttest", analysis_type = "group"),
    "group_col"
  )
  expect_error(
    run_diff(x, method = "ttest", analysis_type = "group", group_col = "group"),
    "control_group"
  )
})

test_that("run_diff requires continuous_col for continuous analysis", {
  x <- make_proteo()
  expect_error(
    run_diff(x, method = "lm", analysis_type = "continuous"),
    "continuous_col"
  )
})

test_that("run_diff anova is limma-only", {
  x <- make_proteo()
  expect_error(
    run_diff(x, method = "ttest", analysis_type = "anova", group_col = "group"),
    "limma"
  )
})

test_that("run_diff edger backend is group-only", {
  x <- make_rna()
  expect_error(
    run_diff(x, method = "edger", analysis_type = "continuous",
             continuous_col = "age"),
    "group"
  )
})

# ---- filter_diff_results ------------------------------------------------

test_that("filter_diff_results applies cutoffs and sets is_significant", {
  x <- make_proteo()
  b <- run_diff(x, method = "ttest", analysis_type = "group",
                group_col = "group", control_group = "ctrl", case_group = "case")
  res <- b$results$diff_result_df

  filt <- filter_diff_results(res, p_cutoff = 0.05)
  expect_true(all(filt$adj_p_value < 0.05))
  expect_true(all(filt$is_significant))
  expect_lte(nrow(filt), nrow(res))

  filt_effect <- filter_diff_results(res, p_cutoff = 0.05, effect_cutoff = 1)
  expect_true(all(abs(filt_effect$effect) >= 1))
  expect_lte(nrow(filt_effect), nrow(filt))
})

test_that("filter_diff_results supports raw p-value preference", {
  x <- make_proteo()
  b <- run_diff(x, method = "ttest", analysis_type = "group",
                group_col = "group", control_group = "ctrl", case_group = "case")
  res <- b$results$diff_result_df
  filt <- filter_diff_results(res, p_cutoff = 0.05, p_preference = "raw")
  expect_true(all(filt$p_value < 0.05))
})

# ---- make_ranked_features ----------------------------------------------

test_that("make_ranked_features returns a decreasing named vector with unique names", {
  x <- make_proteo()
  b <- run_diff(x, method = "ttest", analysis_type = "group",
                group_col = "group", control_group = "ctrl", case_group = "case")
  res <- b$results$diff_result_df

  ranked <- make_ranked_features(res)
  expect_type(ranked, "double")
  expect_true(!is.null(names(ranked)))
  expect_false(anyDuplicated(names(ranked)) > 0)
  expect_equal(ranked, sort(ranked, decreasing = TRUE))
})

test_that("make_ranked_features errors on unknown columns", {
  x <- make_proteo()
  b <- run_diff(x, method = "ttest", analysis_type = "group",
                group_col = "group", control_group = "ctrl", case_group = "case")
  res <- b$results$diff_result_df
  expect_error(make_ranked_features(res, feature_col = "no_such_col"),
               "Feature column")
  expect_error(make_ranked_features(res, rank_col = "no_such_col"),
               "Rank column")
})

# ---- Plot functions: ggplot return types --------------------------------

test_that("plot_volcano returns a ggplot for a t-test bundle", {
  x <- make_proteo()
  b <- run_diff(x, method = "ttest", analysis_type = "group",
                group_col = "group", control_group = "ctrl", case_group = "case")
  p <- plot_volcano(b, top_n = 5, effect_threshold = 1)
  expect_s3_class(p, "ggplot")
})

test_that("plot_volcano honors label_features", {
  x <- make_proteo()
  b <- run_diff(x, method = "ttest", analysis_type = "group",
                group_col = "group", control_group = "ctrl", case_group = "case")
  p <- plot_volcano(b, top_n = 0, label_features = "S1")
  expect_s3_class(p, "ggplot")
})

test_that("plot_ma returns a ggplot when base_mean is populated", {
  x <- make_proteo()
  b <- run_diff(x, method = "lm", analysis_type = "group",
                group_col = "group", control_group = "ctrl", case_group = "case")
  p <- plot_ma(b, top_n = 5)
  expect_s3_class(p, "ggplot")
})

test_that("plot_ma errors when base_mean is all NA", {
  x <- make_proteo()
  b <- run_diff(x, method = "ttest", analysis_type = "group",
                group_col = "group", control_group = "ctrl", case_group = "case")
  # Force base_mean to NA to trigger the guard
  b$results$diff_result_df$base_mean <- NA_real_
  expect_error(plot_ma(b), "base_mean")
})

test_that("plot_pca returns a ggplot and supports color_by", {
  x <- make_proteo()
  p <- plot_pca(x)
  expect_s3_class(p, "ggplot")
  p2 <- plot_pca(x, color_by = "group", shape_by = "sex")
  expect_s3_class(p2, "ggplot")
})

test_that("plot_pca errors on unknown color_by column", {
  x <- make_proteo()
  expect_error(plot_pca(x, color_by = "no_such_col"), "color_by")
})

test_that("plot_feature_expression returns a ggplot for a given feature set", {
  x <- make_proteo()
  p <- plot_feature_expression(x, features = c("g1", "g2"), group_by = "group")
  expect_s3_class(p, "ggplot")
})

test_that("plot_feature_expression resolves feature_symbol as well as feature_id", {
  x <- make_proteo()
  p <- plot_feature_expression(x, features = c("S1", "S2"), group_by = "group")
  expect_s3_class(p, "ggplot")
})

test_that("plot_feature_expression errors when no features match", {
  x <- make_proteo()
  expect_error(
    plot_feature_expression(x, features = "no_such_feature", group_by = "group"),
    "None of the requested features"
  )
})

# ---- Bioconductor backends (Suggests, optional) -------------------------

test_that("run_diff limma group returns a well-formed bundle", {
  skip_if_not_installed("limma")
  x <- make_proteo()
  b <- run_diff(x, method = "limma", analysis_type = "group",
                group_col = "group", control_group = "ctrl", case_group = "case")
  res <- b$results$diff_result_df
  expect_true(all(DIFF_RESULT_REQUIRED_COLS %in% colnames(res)))
  expect_equal(unique(res$method), "limma")
  expect_equal(unique(res$effect_type), "log2FC")
})

test_that("run_diff limma continuous runs and populates correlation effect", {
  skip_if_not_installed("limma")
  x <- make_proteo()
  b <- run_diff_continuous(x, method = "limma", continuous_col = "age")
  res <- b$results$diff_result_df
  expect_equal(unique(res$effect_type), "correlation")
  expect_equal(unique(res$analysis_type), "continuous_linear")
})

test_that("run_diff limma anova returns F-statistic effect", {
  skip_if_not_installed("limma")
  x <- make_proteo()
  # Need at least 3 levels for anova; replace `group`.
  x$meta_df$group3 <- rep(c("A", "B", "C"), length.out = nrow(x$meta_df))
  b <- run_diff(x, method = "limma", analysis_type = "anova",
                group_col = "group3")
  res <- b$results$diff_result_df
  expect_equal(unique(res$effect_type), "F_statistic")
})

test_that("run_diff deseq2 group returns log2FC effect on raw counts", {
  skip_if_not_installed("DESeq2")
  x <- make_rna()
  b <- run_diff(x, method = "deseq2", analysis_type = "group",
                group_col = "group", control_group = "ctrl", case_group = "case")
  res <- b$results$diff_result_df
  expect_equal(unique(res$method), "deseq2")
  expect_equal(unique(res$effect_type), "log2FC")
  hits <- res$feature_id[!is.na(res$adj_p_value) & res$adj_p_value < 0.05]
  expect_gte(length(intersect(hits, paste0("g", 1:5))), 3L)
})

test_that("run_diff edger group returns log2FC effect on raw counts", {
  skip_if_not_installed("edgeR")
  x <- make_rna()
  b <- run_diff(x, method = "edger", analysis_type = "group",
                group_col = "group", control_group = "ctrl", case_group = "case")
  res <- b$results$diff_result_df
  expect_equal(unique(res$method), "edger")
  expect_equal(unique(res$effect_type), "log2FC")
})
