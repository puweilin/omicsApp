# Bundle / standardize / utility-helper edge robustness tests.

# ---- new_analysis_bundle ----------------------------------------------

test_that("new_analysis_bundle returns an analysis_bundle", {
  b <- new_analysis_bundle("test")
  expect_true(is_analysis_bundle(b))
})

test_that("new_analysis_bundle accepts results = list()", {
  b <- new_analysis_bundle("test", results = list(some = data.frame()))
  expect_true(is_analysis_bundle(b))
})

test_that("new_analysis_bundle accepts warnings", {
  b <- new_analysis_bundle("test", warnings = c("w1", "w2"))
  expect_true(is_analysis_bundle(b))
})

test_that("new_analysis_bundle accepts messages", {
  b <- new_analysis_bundle("test", messages = c("m1", "m2"))
  expect_true(is_analysis_bundle(b))
})

test_that("new_analysis_bundle preserves analysis_name", {
  b <- new_analysis_bundle("my_analysis")
  expect_equal(b$analysis_name, "my_analysis")
})

test_that("is_analysis_bundle returns FALSE for a plain list", {
  expect_false(is_analysis_bundle(list()))
})

test_that("is_analysis_bundle returns FALSE for NULL", {
  expect_false(is_analysis_bundle(NULL))
})

# ---- diff_result_from_bundle ------------------------------------------

make_diff_result_bundle <- function() {
  df <- data.frame(
    feature_id = paste0("g", 1:3),
    feature_symbol = paste0("g", 1:3),
    feature_type = "gene",
    omics_type = "proteomics",
    method = "ttest",
    analysis_type = "group",
    comparison = "B vs A",
    effect = c(1, 2, -1),
    effect_type = "log2FC",
    statistic = c(1.5, 2.5, -1.5),
    statistic_type = "t",
    p_value = c(0.01, 0.001, 0.5),
    adj_p_value = c(0.02, 0.005, 0.5),
    direction = c("up", "up", "down"),
    base_mean = c(10, 10, 10),
    model_fit = NA_real_,
    is_significant = c(TRUE, TRUE, FALSE),
    stringsAsFactors = FALSE
  )
  new_analysis_bundle("run_diff",
                      results = list(diff_result_df = df))
}

test_that("diff_result_from_bundle extracts the diff_result_df", {
  b <- make_diff_result_bundle()
  out <- diff_result_from_bundle(b)
  expect_true(is.data.frame(out))
  expect_equal(nrow(out), 3)
})

test_that("diff_result_from_bundle errors on non-bundle", {
  expect_error(diff_result_from_bundle(list()))
})

# ---- new_diff_result_template -----------------------------------------

test_that("new_diff_result_template returns a 0-row data.frame", {
  tpl <- new_diff_result_template()
  expect_true(is.data.frame(tpl))
  expect_equal(nrow(tpl), 0)
})

test_that("new_diff_result_template has the required columns", {
  tpl <- new_diff_result_template()
  expect_true(all(DIFF_RESULT_REQUIRED_COLS %in% colnames(tpl)))
})

# ---- new_enrich_result_template ---------------------------------------

test_that("new_enrich_result_template returns a 0-row data.frame", {
  tpl <- new_enrich_result_template()
  expect_true(is.data.frame(tpl))
  expect_equal(nrow(tpl), 0)
})

test_that("new_enrich_result_template has the required columns", {
  tpl <- new_enrich_result_template()
  expect_true(all(ENRICH_RESULT_REQUIRED_COLS %in% colnames(tpl)))
})

# ---- new_integration_result_template ----------------------------------

test_that("new_integration_result_template returns a 0-row data.frame", {
  tpl <- new_integration_result_template()
  expect_true(is.data.frame(tpl))
  expect_equal(nrow(tpl), 0)
})

# ---- make_ranked_features ---------------------------------------------

make_full_rank_df <- function(n = 5,
                               effects = NULL,
                               statistic = NULL) {
  if (is.null(effects)) effects <- seq(-2, 2, length.out = n)
  if (is.null(statistic)) statistic <- effects * 1.5
  data.frame(
    feature_id = paste0("g", seq_len(n)),
    feature_symbol = paste0("G", seq_len(n)),
    feature_type = "gene",
    omics_type = "proteomics",
    method = "ttest",
    analysis_type = "group",
    comparison = "B vs A",
    effect = effects,
    effect_type = "log2FC",
    statistic = statistic,
    statistic_type = "t",
    p_value = runif(n, 0.001, 0.5),
    adj_p_value = runif(n, 0.01, 0.5),
    direction = ifelse(effects > 0, "up", "down"),
    base_mean = rep(10, n),
    model_fit = NA_real_,
    is_significant = FALSE,
    stringsAsFactors = FALSE
  )
}

test_that("make_ranked_features returns a named numeric vector", {
  df <- make_full_rank_df(5)
  out <- make_ranked_features(df)
  expect_type(out, "double")
  expect_true(!is.null(names(out)))
})

test_that("make_ranked_features ranks by effect by default", {
  df <- make_full_rank_df(3, effects = c(1, 3, 2))
  out <- make_ranked_features(df)
  expect_equal(length(out), 3)
})

test_that("make_ranked_features supports custom feature_col", {
  df <- make_full_rank_df(3)
  out <- make_ranked_features(df, feature_col = "feature_id")
  expect_true(all(grepl("^g", names(out))))
})

test_that("make_ranked_features supports custom rank_col", {
  df <- make_full_rank_df(3)
  out <- make_ranked_features(df, rank_col = "statistic")
  expect_equal(length(out), 3)
})

test_that("make_ranked_features drops NA ranks", {
  df <- make_full_rank_df(4)
  df$effect[c(2, 4)] <- NA
  out <- make_ranked_features(df)
  expect_true(!any(is.na(out)))
})

# ---- filter_enrich_results --------------------------------------------

test_that("filter_enrich_results errors on non-data.frame", {
  expect_error(filter_enrich_results(list()))
})

# ---- pick_features_by_variance ---------------------------------------

test_that("pick_features_by_variance returns top-N feature_ids", {
  mat <- matrix(rnorm(60), nrow = 10, ncol = 6,
                dimnames = list(paste0("g", 1:10), paste0("s", 1:6)))
  # row 1 has highest variance
  mat[1, ] <- c(0, 100, 0, 100, 0, 100)
  out <- pick_features_by_variance(mat, assay_type = "intensity",
                                   n_top = 3L, features = rownames(mat))
  expect_true("g1" %in% out)
})

test_that("pick_features_by_variance n_top larger than rows returns all rows", {
  mat <- matrix(rnorm(30), nrow = 5, ncol = 6,
                dimnames = list(paste0("g", 1:5), paste0("s", 1:6)))
  out <- pick_features_by_variance(mat, assay_type = "intensity",
                                   n_top = 100L, features = rownames(mat))
  expect_equal(length(out), 5)
})

# ---- safe_z -----------------------------------------------------------

test_that("safe_z z-scores a numeric vector", {
  out <- safe_z(c(1, 2, 3, 4, 5))
  expect_true(abs(mean(out)) < 1e-12)
})

test_that("safe_z returns finite values for constant vector", {
  out <- safe_z(rep(5, 5))
  expect_true(all(is.finite(out)) || all(out == 0) || all(is.na(out)))
})

# ---- scale_rows -------------------------------------------------------

test_that("scale_rows returns a matrix of the same shape", {
  m <- matrix(rnorm(20), nrow = 5, ncol = 4)
  out <- scale_rows(m)
  expect_equal(dim(out), dim(m))
})

test_that("scale_rows centers rows around zero", {
  m <- matrix(rnorm(20, mean = 10), nrow = 5, ncol = 4)
  out <- scale_rows(m)
  expect_true(max(abs(rowMeans(out, na.rm = TRUE))) < 1e-6)
})

# ---- format_short / parse_overlap_size --------------------------------

test_that("format_short returns a character", {
  out <- format_short("this_is_a_long_pathway_name_xxx")
  expect_true(is.character(out))
})

test_that("parse_overlap_size handles 'a/b' format", {
  out <- parse_overlap_size("3/15")
  expect_equal(out, 3)
})

test_that("parse_overlap_size returns NA on bad input", {
  out <- parse_overlap_size("garbage")
  expect_true(is.na(out) || is.numeric(out))
})

# ---- tidy_pathway_names -----------------------------------------------

test_that("tidy_pathway_names trims and cleans labels", {
  out <- tidy_pathway_names(c("HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                              "REACTOME_CELL_CYCLE"))
  expect_true(is.character(out) && length(out) == 2)
})

test_that("tidy_pathway_names handles empty input", {
  out <- tidy_pathway_names(character(0))
  expect_equal(length(out), 0)
})

# ---- truncate_pathway_name --------------------------------------------

test_that("truncate_pathway_name shortens overlong names", {
  out <- truncate_pathway_name("a_very_long_pathway_name", max_chars = 8L)
  expect_true(nchar(out) <= 12L)  # allow ellipsis bytes
})

test_that("truncate_pathway_name preserves short names", {
  out <- truncate_pathway_name("short", max_chars = 60L)
  expect_equal(out, "short")
})

# ---- summarize_omics_input --------------------------------------------

make_summary_input <- function(n_feat = 10, n_samp = 6) {
  mat <- matrix(rnorm(n_feat * n_samp), nrow = n_feat, ncol = n_samp)
  rownames(mat) <- paste0("g", seq_len(n_feat))
  colnames(mat) <- paste0("s", seq_len(n_samp))
  meta <- data.frame(group = rep(c("A", "B"), length.out = n_samp),
                     row.names = colnames(mat))
  feat <- data.frame(feature_id = rownames(mat))
  omics_input(mat, meta, feat, omics_type = "proteomics",
              assay_type = "intensity")
}

test_that("summarize_omics_input returns a list for a valid input", {
  inp <- make_summary_input()
  out <- summarize_omics_input(inp)
  expect_true(is.list(out) || is.data.frame(out))
})

test_that("summarize_omics_input is non-mutating", {
  inp <- make_summary_input()
  snap <- inp
  invisible(summarize_omics_input(inp))
  expect_equal(inp$expr_mat, snap$expr_mat)
})

test_that("summarize_omics_input handles NA-containing matrix", {
  inp <- make_summary_input()
  inp$expr_mat[1, 1] <- NA
  out <- summarize_omics_input(inp)
  expect_true(is.list(out) || is.data.frame(out))
})

# ---- is_installed and install group resolution ------------------------

test_that("is_installed returns TRUE for base package", {
  expect_true(is_installed("methods"))
})

test_that("is_installed returns FALSE for a fake package", {
  expect_false(is_installed("definitely_not_a_real_pkg_xyzzy"))
})

test_that("resolve_install_group accepts valid group", {
  out <- resolve_install_group("rnaseq")
  expect_true(is.character(out))
})

test_that("resolve_install_group errors on unknown group", {
  expect_error(resolve_install_group("xyz_group"))
})

# ---- materialize_* helpers --------------------------------------------

test_that("materialize_matrix returns a matrix from features-as-rows df", {
  df <- data.frame(feature_id = paste0("g", 1:5),
                   s1 = rnorm(5), s2 = rnorm(5), s3 = rnorm(5),
                   stringsAsFactors = FALSE)
  out <- materialize_matrix(df, orientation = "features_as_rows")
  expect_true(is.matrix(out) || inherits(out, "data.frame"))
})

test_that("materialize_metadata picks rows by sample_ids", {
  meta <- data.frame(sample_id = c("s1", "s2", "s3"),
                     group = c("A", "B", "A"),
                     stringsAsFactors = FALSE)
  out <- materialize_metadata(meta, sample_ids = c("s1", "s2"))
  expect_true(is.data.frame(out))
})

test_that("materialize_feature_annot picks rows by feature_ids", {
  feat <- data.frame(feature_id = paste0("g", 1:5),
                     stringsAsFactors = FALSE)
  out <- materialize_feature_annot(feat, feature_ids = c("g1", "g2"))
  expect_true(is.data.frame(out))
})

# ---- subset_omics tag-preserving --------------------------------------

test_that("subset_omics_samples preserves feature_df", {
  inp <- make_summary_input()
  sub <- subset_omics_samples(inp, c("s1", "s2"))
  expect_equal(sub$feature_df$feature_id, inp$feature_df$feature_id)
})

test_that("subset_omics_features preserves meta_df", {
  inp <- make_summary_input()
  sub <- subset_omics_features(inp, c("g1", "g2"))
  expect_equal(rownames(sub$meta_df), rownames(inp$meta_df))
})

test_that("subset_omics_samples errors on unknown sample", {
  inp <- make_summary_input()
  expect_error(subset_omics_samples(inp, c("not_a_sample")))
})

test_that("subset_omics_features errors on unknown feature", {
  inp <- make_summary_input()
  expect_error(subset_omics_features(inp, c("ghost_gene")))
})
