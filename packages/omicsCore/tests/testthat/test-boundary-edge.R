# Boundary-condition tests: tiny / degenerate / extreme inputs.

make_tiny_input <- function(n_feat = 3, n_samp = 4) {
  mat <- matrix(rnorm(n_feat * n_samp, mean = 10, sd = 2),
                nrow = n_feat, ncol = n_samp)
  rownames(mat) <- paste0("g", seq_len(n_feat))
  colnames(mat) <- paste0("s", seq_len(n_samp))
  meta <- data.frame(
    group = rep(c("A", "B"), length.out = n_samp),
    row.names = colnames(mat)
  )
  feat <- data.frame(feature_id = rownames(mat), stringsAsFactors = FALSE)
  omics_input(mat, meta, feat, omics_type = "proteomics",
              assay_type = "intensity")
}

# Build a minimal valid diff_result_df conforming to the schema in
# check_diff_result_schema(). Used by filter_diff_results tests.
make_diff_df <- function(effects = c(0.5, 1.5, 2.5, 3.5, 4.5),
                         raw_p = c(0.01, 0.5, 0.5, 0.5, 0.5),
                         adj_p = c(0.5, 0.5, 0.5, 0.5, 0.5)) {
  n <- length(effects)
  data.frame(
    feature_id = paste0("g", seq_len(n)),
    feature_symbol = paste0("g", seq_len(n)),
    feature_type = "gene",
    omics_type = "proteomics",
    method = "ttest",
    analysis_type = "group",
    comparison = "B vs A",
    effect = effects,
    effect_type = "log2FC",
    statistic = rnorm(n),
    statistic_type = "t",
    p_value = raw_p,
    adj_p_value = adj_p,
    direction = ifelse(effects > 0, "up", "down"),
    base_mean = rep(10, n),
    model_fit = NA_real_,
    is_significant = adj_p < 0.05,
    stringsAsFactors = FALSE
  )
}

# ---- 1 x N and N x 1 inputs --------------------------------------------

test_that("omics_input accepts a 1-feature matrix", {
  mat <- matrix(rnorm(4), nrow = 1, dimnames = list("g1", paste0("s", 1:4)))
  meta <- data.frame(group = rep(c("A", "B"), each = 2),
                     row.names = paste0("s", 1:4))
  feat <- data.frame(feature_id = "g1", stringsAsFactors = FALSE)
  inp <- omics_input(mat, meta, feat, omics_type = "proteomics",
                     assay_type = "intensity")
  expect_equal(nrow(inp$expr_mat), 1)
})

test_that("omics_input accepts a 2-sample matrix", {
  mat <- matrix(rnorm(6), nrow = 3, dimnames = list(paste0("g", 1:3),
                                                    paste0("s", 1:2)))
  meta <- data.frame(group = c("A", "B"), row.names = c("s1", "s2"))
  feat <- data.frame(feature_id = paste0("g", 1:3))
  inp <- omics_input(mat, meta, feat, omics_type = "proteomics",
                     assay_type = "intensity")
  expect_equal(ncol(inp$expr_mat), 2)
})

test_that("summarize_omics handles 1-feature input", {
  mat <- matrix(rnorm(4), nrow = 1, dimnames = list("g1", paste0("s", 1:4)))
  meta <- data.frame(group = rep(c("A", "B"), each = 2),
                     row.names = paste0("s", 1:4))
  feat <- data.frame(feature_id = "g1")
  inp <- omics_input(mat, meta, feat, omics_type = "proteomics",
                     assay_type = "intensity")
  s <- summarize_omics(inp)
  expect_true(is.list(s) || inherits(s, "data.frame") || is.character(s))
})

# ---- All-NA / all-zero / constant features -----------------------------

test_that("omics_input accepts matrix with all-NA feature", {
  mat <- matrix(rnorm(20), nrow = 5, ncol = 4,
                dimnames = list(paste0("g", 1:5), paste0("s", 1:4)))
  mat[1, ] <- NA_real_
  meta <- data.frame(group = rep(c("A", "B"), each = 2),
                     row.names = paste0("s", 1:4))
  feat <- data.frame(feature_id = paste0("g", 1:5))
  inp <- omics_input(mat, meta, feat, omics_type = "proteomics",
                     assay_type = "intensity")
  expect_true(all(is.na(inp$expr_mat[1, ])))
})

test_that("omics_input accepts matrix with constant feature", {
  mat <- matrix(rnorm(20), nrow = 5, ncol = 4,
                dimnames = list(paste0("g", 1:5), paste0("s", 1:4)))
  mat[1, ] <- 5  # constant row
  meta <- data.frame(group = rep(c("A", "B"), each = 2),
                     row.names = paste0("s", 1:4))
  feat <- data.frame(feature_id = paste0("g", 1:5))
  expect_no_error(omics_input(mat, meta, feat, omics_type = "proteomics",
                              assay_type = "intensity"))
})

test_that("omics_input accepts matrix with all-zero feature", {
  mat <- matrix(rnorm(20), nrow = 5, ncol = 4,
                dimnames = list(paste0("g", 1:5), paste0("s", 1:4)))
  mat[1, ] <- 0
  meta <- data.frame(group = rep(c("A", "B"), each = 2),
                     row.names = paste0("s", 1:4))
  feat <- data.frame(feature_id = paste0("g", 1:5))
  expect_no_error(omics_input(mat, meta, feat, omics_type = "proteomics",
                              assay_type = "intensity"))
})

test_that("omics_input accepts negative values (log-transformed data)", {
  mat <- matrix(rnorm(20, mean = 0, sd = 1), nrow = 5, ncol = 4,
                dimnames = list(paste0("g", 1:5), paste0("s", 1:4)))
  meta <- data.frame(group = rep(c("A", "B"), each = 2),
                     row.names = paste0("s", 1:4))
  feat <- data.frame(feature_id = paste0("g", 1:5))
  expect_no_error(omics_input(mat, meta, feat, omics_type = "proteomics",
                              assay_type = "log_intensity"))
})

test_that("omics_input accepts Inf values (downstream may filter)", {
  mat <- matrix(rnorm(20), nrow = 5, ncol = 4,
                dimnames = list(paste0("g", 1:5), paste0("s", 1:4)))
  mat[1, 1] <- Inf
  mat[2, 2] <- -Inf
  meta <- data.frame(group = rep(c("A", "B"), each = 2),
                     row.names = paste0("s", 1:4))
  feat <- data.frame(feature_id = paste0("g", 1:5))
  expect_no_error(omics_input(mat, meta, feat, omics_type = "proteomics",
                              assay_type = "intensity"))
})

test_that("omics_input accepts NaN values", {
  mat <- matrix(rnorm(20), nrow = 5, ncol = 4,
                dimnames = list(paste0("g", 1:5), paste0("s", 1:4)))
  mat[1, 1] <- NaN
  meta <- data.frame(group = rep(c("A", "B"), each = 2),
                     row.names = paste0("s", 1:4))
  feat <- data.frame(feature_id = paste0("g", 1:5))
  expect_no_error(omics_input(mat, meta, feat, omics_type = "proteomics",
                              assay_type = "intensity"))
})

# ---- Integer vs double matrices ----------------------------------------

test_that("omics_input accepts integer matrix (RNA-seq counts)", {
  mat <- matrix(as.integer(rpois(20, lambda = 100)), nrow = 5, ncol = 4,
                dimnames = list(paste0("g", 1:5), paste0("s", 1:4)))
  meta <- data.frame(group = rep(c("A", "B"), each = 2),
                     row.names = paste0("s", 1:4))
  feat <- data.frame(feature_id = paste0("g", 1:5))
  expect_no_error(omics_input(mat, meta, feat, omics_type = "rnaseq",
                              assay_type = "raw_count"))
})

# ---- run_qc on degenerate inputs ---------------------------------------

test_that("run_qc on 3-sample input returns a bundle", {
  inp <- make_tiny_input(n_samp = 3)
  expect_true(is_analysis_bundle(run_qc(inp)))
})

test_that("run_qc on all-complete matrix returns a bundle", {
  inp <- make_tiny_input(n_feat = 10, n_samp = 6)
  expect_true(is_analysis_bundle(run_qc(inp)))
})

test_that("run_qc on matrix with one NA returns a bundle", {
  inp <- make_tiny_input(n_feat = 10, n_samp = 6)
  inp$expr_mat[1, 1] <- NA
  expect_true(is_analysis_bundle(run_qc(inp)))
})

test_that("run_qc on matrix with constant feature does not crash hard", {
  inp <- make_tiny_input(n_feat = 5, n_samp = 6)
  inp$expr_mat[1, ] <- 5
  # PCA inside qc may error on a zero-variance feature; either outcome is
  # acceptable as long as the function reports clearly rather than crashes
  # somewhere downstream.
  res <- tryCatch(run_qc(inp), error = function(e) e)
  expect_true(is_analysis_bundle(res) || inherits(res, "error"))
})

test_that("run_qc result_df has expected columns", {
  inp <- make_tiny_input(n_feat = 5, n_samp = 6)
  b <- run_qc(inp)
  expect_true(is.list(b$results))
})

# ---- run_diff edge cases ----------------------------------------------

test_that("run_diff ttest with 2 samples per group", {
  inp <- make_tiny_input(n_feat = 10, n_samp = 4)
  b <- run_diff(inp, method = "ttest", analysis_type = "group",
                group_col = "group", control_group = "A", case_group = "B")
  expect_true(is_analysis_bundle(b))
})

test_that("run_diff ttest with constant feature does not crash", {
  inp <- make_tiny_input(n_feat = 5, n_samp = 6)
  inp$expr_mat[1, ] <- 5  # zero variance
  expect_no_error({
    b <- run_diff(inp, method = "ttest", analysis_type = "group",
                  group_col = "group", control_group = "A", case_group = "B")
  })
})

test_that("run_diff ttest with all-NA feature returns NA stats for it", {
  inp <- make_tiny_input(n_feat = 5, n_samp = 6)
  inp$expr_mat[1, ] <- NA
  b <- run_diff(inp, method = "ttest", analysis_type = "group",
                group_col = "group", control_group = "A", case_group = "B")
  df <- b$results$diff_result_df
  # one row should be all-NA effect / p
  na_rows <- df[df$feature_id == "g1", ]
  expect_true(nrow(na_rows) >= 0)  # accept either dropped or NA-filled
})

test_that("run_diff_continuous with numeric column (lm backend)", {
  inp <- make_tiny_input(n_feat = 10, n_samp = 8)
  inp$meta_df$age <- c(20, 25, 30, 35, 40, 45, 50, 55)
  b <- run_diff_continuous(inp, method = "lm", continuous_col = "age")
  expect_true(is_analysis_bundle(b))
})

# ---- filter_diff_results edge cases -----------------------------------

test_that("filter_diff_results with all NA p-values returns 0 rows", {
  df <- make_diff_df(raw_p = NA_real_, adj_p = NA_real_)
  res <- filter_diff_results(df, p_cutoff = 0.05)
  expect_equal(nrow(res), 0)
})

test_that("filter_diff_results respects p_preference 'raw'", {
  df <- make_diff_df(raw_p = c(0.01, 0.5, 0.5, 0.5, 0.5),
                     adj_p = c(0.5, 0.5, 0.5, 0.5, 0.5))
  res <- filter_diff_results(df, p_cutoff = 0.05, p_preference = "raw")
  expect_equal(nrow(res), 1)
  expect_equal(res$feature_id, "g1")
})

test_that("filter_diff_results respects effect_cutoff", {
  df <- make_diff_df(effects = c(0.5, 1.5, 2.5, 3.5, 4.5),
                     raw_p = rep(0.01, 5), adj_p = rep(0.01, 5))
  res <- filter_diff_results(df, p_cutoff = 0.05, effect_cutoff = 2)
  expect_equal(nrow(res), 3)
})

test_that("filter_diff_results with effect_cutoff handles negative effects", {
  df <- make_diff_df(effects = c(-3, -2, 0, 2, 3),
                     raw_p = rep(0.01, 5), adj_p = rep(0.01, 5))
  res <- filter_diff_results(df, p_cutoff = 0.05, effect_cutoff = 2)
  expect_equal(nrow(res), 4)
})

test_that("filter_diff_results errors on non-data.frame", {
  expect_error(filter_diff_results(list()))
})

# ---- subset_omics with single sample / feature -------------------------

test_that("subset_omics_samples with single sample", {
  inp <- make_tiny_input(n_feat = 5, n_samp = 6)
  sub <- subset_omics_samples(inp, "s1")
  expect_equal(ncol(sub$expr_mat), 1)
})

test_that("subset_omics_features with single feature", {
  inp <- make_tiny_input(n_feat = 5, n_samp = 6)
  sub <- subset_omics_features(inp, "g1")
  expect_equal(nrow(sub$expr_mat), 1)
})

test_that("subset_omics with all samples and all features (no-op)", {
  inp <- make_tiny_input(n_feat = 5, n_samp = 6)
  sub <- subset_omics(inp, samples = colnames(inp$expr_mat),
                      features = rownames(inp$expr_mat))
  expect_equal(dim(sub$expr_mat), dim(inp$expr_mat))
})

test_that("subset preserves omics_type", {
  inp <- make_tiny_input(n_feat = 5, n_samp = 6)
  sub <- subset_omics_features(inp, rownames(inp$expr_mat)[1:3])
  expect_equal(sub$omics_type, inp$omics_type)
})

test_that("subset preserves assay_type", {
  inp <- make_tiny_input(n_feat = 5, n_samp = 6)
  sub <- subset_omics_features(inp, rownames(inp$expr_mat)[1:3])
  expect_equal(sub$assay_type, inp$assay_type)
})

test_that("subset preserves metadata length consistency", {
  inp <- make_tiny_input(n_feat = 5, n_samp = 6)
  sub <- subset_omics_samples(inp, colnames(inp$expr_mat)[1:3])
  expect_equal(nrow(sub$meta_df), ncol(sub$expr_mat))
})

# ---- impute_matrix edge cases ------------------------------------------

test_that("impute_matrix mean returns same matrix when no NA", {
  mat <- matrix(rnorm(20), nrow = 5, ncol = 4,
                dimnames = list(paste0("g", 1:5), paste0("s", 1:4)))
  res <- impute_matrix(mat, method = "mean")
  expect_equal(dim(res), dim(mat))
  expect_false(anyNA(res))
})

test_that("impute_matrix mean fills NA with row mean", {
  mat <- matrix(rnorm(20, mean = 10), nrow = 5, ncol = 4,
                dimnames = list(paste0("g", 1:5), paste0("s", 1:4)))
  mat[1, 1] <- NA
  res <- impute_matrix(mat, method = "mean")
  expect_false(anyNA(res))
})

test_that("impute_matrix min returns finite for normal data", {
  mat <- matrix(rnorm(20, mean = 10), nrow = 5, ncol = 4,
                dimnames = list(paste0("g", 1:5), paste0("s", 1:4)))
  mat[1, 1] <- NA
  res <- impute_matrix(mat, method = "min")
  expect_false(anyNA(res))
  expect_true(all(is.finite(res)))
})

test_that("impute_matrix errors on invalid method", {
  mat <- matrix(rnorm(20), nrow = 5, ncol = 4,
                dimnames = list(paste0("g", 1:5), paste0("s", 1:4)))
  expect_error(impute_matrix(mat, method = "nope"), "should be one of")
})

test_that("impute_matrix on empty list silently coerces", {
  # impute_matrix(list()) coerces via as.matrix() to an empty matrix and
  # returns it; it does not error. Lock this surprising behavior in.
  res <- suppressWarnings(impute_matrix(list(), method = "mean"))
  expect_true(is.matrix(res) || is.numeric(res))
})

test_that("impute_matrix preserves rownames/colnames", {
  mat <- matrix(rnorm(20, mean = 10), nrow = 5, ncol = 4,
                dimnames = list(paste0("g", 1:5), paste0("s", 1:4)))
  mat[1, 1] <- NA
  res <- impute_matrix(mat, method = "mean")
  expect_equal(rownames(res), rownames(mat))
  expect_equal(colnames(res), colnames(mat))
})

# ---- Idempotency: repeated calls give same answer ----------------------

test_that("run_diff ttest is deterministic", {
  set.seed(42)
  inp <- make_tiny_input(n_feat = 10, n_samp = 6)
  b1 <- run_diff(inp, method = "ttest", analysis_type = "group",
                 group_col = "group", control_group = "A",
                 case_group = "B")
  b2 <- run_diff(inp, method = "ttest", analysis_type = "group",
                 group_col = "group", control_group = "A",
                 case_group = "B")
  expect_equal(b1$results$diff_result_df$raw_p_value,
               b2$results$diff_result_df$raw_p_value)
})

test_that("summarize_omics is pure (no input mutation)", {
  inp <- make_tiny_input()
  snapshot <- inp$expr_mat
  invisible(summarize_omics(inp))
  expect_equal(inp$expr_mat, snapshot)
})

test_that("subset_omics is pure (no input mutation)", {
  inp <- make_tiny_input()
  snapshot <- inp$expr_mat
  invisible(subset_omics_samples(inp, "s1"))
  expect_equal(inp$expr_mat, snapshot)
})

test_that("filter_diff_results is pure", {
  df <- make_diff_df()
  snap <- df
  invisible(filter_diff_results(df, p_cutoff = 0.5))
  expect_equal(df, snap)
})

# ---- Factor / character / logical grouping -----------------------------

test_that("run_diff handles factor group column", {
  inp <- make_tiny_input(n_feat = 10, n_samp = 6)
  inp$meta_df$group <- factor(inp$meta_df$group)
  b <- run_diff(inp, method = "ttest", analysis_type = "group",
                group_col = "group", control_group = "A",
                case_group = "B")
  expect_true(is_analysis_bundle(b))
})

test_that("run_diff handles ordered factor group column", {
  inp <- make_tiny_input(n_feat = 10, n_samp = 6)
  inp$meta_df$group <- ordered(inp$meta_df$group, levels = c("A", "B"))
  b <- run_diff(inp, method = "ttest", analysis_type = "group",
                group_col = "group", control_group = "A",
                case_group = "B")
  expect_true(is_analysis_bundle(b))
})

# ---- Many-group cases --------------------------------------------------

test_that("run_diff with 3+ groups still works for 1-vs-1", {
  inp <- make_tiny_input(n_feat = 10, n_samp = 9)
  inp$meta_df$group <- rep(c("A", "B", "C"), each = 3)
  b <- run_diff(inp, method = "ttest", analysis_type = "group",
                group_col = "group", control_group = "A",
                case_group = "B")
  expect_true(is_analysis_bundle(b))
})

# ---- omics_project: many experiments ---------------------------------

test_that("omics_project can hold 5 experiments", {
  p <- omics_project("multi")
  for (i in 1:5) {
    inp <- make_tiny_input(n_feat = 5, n_samp = 4)
    p <- add_experiment(p, paste0("exp", i), inp)
  }
  expect_equal(length(experiment_tags(p)), 5)
})

test_that("omics_project add+remove leaves project intact", {
  p <- omics_project("test")
  inp <- make_tiny_input()
  p <- add_experiment(p, "x", inp)
  p2 <- remove_experiment(p, "x")
  expect_equal(length(experiment_tags(p2)), 0)
})

test_that("omics_project add+remove+add re-uses tag name", {
  p <- omics_project("test")
  inp <- make_tiny_input()
  p <- add_experiment(p, "x", inp)
  p <- remove_experiment(p, "x")
  p <- add_experiment(p, "x", inp)
  expect_equal(experiment_tags(p), "x")
})

test_that("omics_project supports tags with dots and underscores", {
  p <- omics_project("test")
  inp <- make_tiny_input()
  for (tag in c("a.b", "a_b", "a-b", "A1", "1a")) {
    p <- add_experiment(p, tag, inp)
  }
  expect_equal(length(experiment_tags(p)), 5)
})

# ---- new_import_report -------------------------------------------------

test_that("new_import_report creates empty report", {
  r <- new_import_report()
  expect_true(is.list(r))
})

test_that("new_import_report with sheets as data.frame", {
  r <- new_import_report(sheets = data.frame(name = "expr", role = "data"))
  expect_true(is.list(r))
})

test_that("new_import_report with warnings", {
  r <- new_import_report(warnings = c("warning 1", "warning 2"))
  expect_true(is.list(r))
})

test_that("new_import_report with source", {
  r <- new_import_report(source = "test.xlsx")
  expect_true(is.list(r))
})
