# Validation & schema-check robustness tests.

make_valid_input <- function(n_feat = 5, n_samp = 6) {
  mat <- matrix(rnorm(n_feat * n_samp), nrow = n_feat, ncol = n_samp)
  rownames(mat) <- paste0("g", seq_len(n_feat))
  colnames(mat) <- paste0("s", seq_len(n_samp))
  half <- n_samp %/% 2
  meta <- data.frame(
    group = c(rep("A", half), rep("B", n_samp - half)),
    paired = c(seq_len(half), seq_len(n_samp - half)),
    row.names = colnames(mat)
  )
  feat <- data.frame(feature_id = rownames(mat), stringsAsFactors = FALSE)
  omics_input(mat, meta, feat, omics_type = "proteomics",
              assay_type = "intensity")
}

make_full_diff_df <- function(n = 5) {
  data.frame(
    feature_id = paste0("g", seq_len(n)),
    feature_symbol = paste0("g", seq_len(n)),
    feature_type = "gene",
    omics_type = "proteomics",
    method = "ttest",
    analysis_type = "group",
    comparison = "B vs A",
    effect = rnorm(n),
    effect_type = "log2FC",
    statistic = rnorm(n),
    statistic_type = "t",
    p_value = runif(n),
    adj_p_value = runif(n),
    direction = "up",
    base_mean = rep(10, n),
    model_fit = NA_real_,
    is_significant = FALSE,
    stringsAsFactors = FALSE
  )
}

# ---- validate_omics_input ----------------------------------------------

test_that("validate_omics_input passes for a valid input", {
  inp <- make_valid_input()
  expect_no_error(validate_omics_input(inp))
})

test_that("validate_omics_input errors on a plain list", {
  expect_error(validate_omics_input(list()))
})

test_that("validate_omics_input errors on NULL", {
  expect_error(validate_omics_input(NULL))
})

test_that("is_omics_input returns TRUE for a valid input", {
  inp <- make_valid_input()
  expect_true(is_omics_input(inp))
})

test_that("is_omics_input returns FALSE for a plain list", {
  expect_false(is_omics_input(list()))
})

test_that("is_omics_input returns FALSE for NULL", {
  expect_false(is_omics_input(NULL))
})

# ---- check_required_cols -----------------------------------------------

test_that("check_required_cols passes when all cols present", {
  df <- data.frame(a = 1, b = 2, c = 3)
  expect_no_error(check_required_cols(df, c("a", "b")))
})

test_that("check_required_cols errors on missing col", {
  df <- data.frame(a = 1, b = 2)
  expect_error(check_required_cols(df, c("a", "c")))
})

test_that("check_required_cols message mentions missing col name", {
  df <- data.frame(a = 1)
  err <- tryCatch(check_required_cols(df, c("a", "missing_x")),
                  error = function(e) e$message)
  expect_match(err, "missing_x")
})

test_that("check_required_cols accepts custom object_name in message", {
  df <- data.frame(a = 1)
  err <- tryCatch(
    check_required_cols(df, "b", object_name = "my_table"),
    error = function(e) e$message)
  expect_match(err, "my_table")
})

# ---- check_diff_result_schema -----------------------------------------

test_that("check_diff_result_schema passes for full diff df", {
  df <- make_full_diff_df()
  expect_no_error(check_diff_result_schema(df))
})

test_that("check_diff_result_schema errors on missing required col", {
  df <- make_full_diff_df()
  df$feature_id <- NULL
  expect_error(check_diff_result_schema(df))
})

test_that("check_diff_result_schema errors on non-data.frame", {
  expect_error(check_diff_result_schema(list()))
})

# ---- check_enrich_result_schema ---------------------------------------

test_that("check_enrich_result_schema errors on empty df", {
  expect_error(check_enrich_result_schema(data.frame()))
})

test_that("check_enrich_result_schema errors on non-data.frame", {
  expect_error(check_enrich_result_schema(list()))
})

# ---- check_integration_result_schema ----------------------------------

test_that("check_integration_result_schema errors on empty df", {
  expect_error(check_integration_result_schema(data.frame()))
})

# ---- check_meta_matches_expr ------------------------------------------

test_that("check_meta_matches_expr passes when meta rows = expr cols", {
  inp <- make_valid_input()
  expect_no_error(check_meta_matches_expr(inp$expr_mat, inp$meta_df))
})

test_that("check_meta_matches_expr errors when row count mismatch", {
  inp <- make_valid_input()
  bad_meta <- inp$meta_df[1:3, , drop = FALSE]
  expect_error(check_meta_matches_expr(inp$expr_mat, bad_meta))
})

test_that("check_meta_matches_expr errors on rowname mismatch", {
  inp <- make_valid_input()
  bad_meta <- inp$meta_df
  rownames(bad_meta) <- paste0("x", seq_len(nrow(bad_meta)))
  expect_error(check_meta_matches_expr(inp$expr_mat, bad_meta))
})

# ---- check_feature_matches_expr ---------------------------------------

test_that("check_feature_matches_expr passes when rows match", {
  inp <- make_valid_input()
  expect_no_error(check_feature_matches_expr(inp$expr_mat, inp$feature_df))
})

test_that("check_feature_matches_expr errors when row count mismatch", {
  inp <- make_valid_input()
  bad_feat <- inp$feature_df[1:2, , drop = FALSE]
  expect_error(check_feature_matches_expr(inp$expr_mat, bad_feat))
})

# ---- validate_two_group_pairing ---------------------------------------

test_that("validate_two_group_pairing with NULL paired_col passes regardless", {
  inp <- make_valid_input()
  # NULL paired_col → short-circuits to TRUE, even with mismatched levels.
  expect_no_error(
    validate_two_group_pairing(inp$meta_df, group_col = "group",
                               paired_col = NULL,
                               control_group = "ZZ", case_group = "B")
  )
})

test_that("validate_two_group_pairing with paired_col errors on broken pairs", {
  inp <- make_valid_input()
  # Break pairing: 3 ctrl samples paired with 3 case samples but
  # case_group label doesn't match group column.
  expect_error(
    validate_two_group_pairing(inp$meta_df, group_col = "group",
                               paired_col = "paired",
                               control_group = "A", case_group = "ZZ")
  )
})

test_that("validate_two_group_pairing with bad paired_col name errors", {
  inp <- make_valid_input()
  expect_error(
    validate_two_group_pairing(inp$meta_df, group_col = "group",
                               paired_col = "nope_col",
                               control_group = "A", case_group = "B")
  )
})

test_that("validate_two_group_pairing supports a paired column", {
  inp <- make_valid_input()
  expect_no_error(
    validate_two_group_pairing(inp$meta_df, group_col = "group",
                               paired_col = "paired",
                               control_group = "A", case_group = "B")
  )
})

# ---- validate_diff_args ------------------------------------------------

test_that("validate_diff_args passes for valid ttest/group args", {
  expect_no_error(
    validate_diff_args(analysis_type = "group", method = "ttest",
                       args = list(group_col = "group",
                                   control_group = "A",
                                   case_group = "B"))
  )
})

test_that("validate_diff_args errors when group_col missing for group analysis", {
  expect_error(
    validate_diff_args(analysis_type = "group", method = "ttest",
                       args = list(control_group = "A", case_group = "B"))
  )
})

test_that("validate_diff_args errors when continuous_col missing", {
  expect_error(
    validate_diff_args(analysis_type = "continuous", method = "lm",
                       args = list())
  )
})

test_that("validate_diff_args errors when anova uses non-limma method", {
  expect_error(
    validate_diff_args(analysis_type = "anova", method = "ttest",
                       args = list(group_col = "g"))
  )
})

test_that("validate_diff_args errors when edger used for continuous", {
  expect_error(
    validate_diff_args(analysis_type = "continuous", method = "edger",
                       args = list(continuous_col = "age"))
  )
})

# ---- validate_sample_link ---------------------------------------------

test_that("validate_sample_link errors on non-data.frame input", {
  expect_error(validate_sample_link(list(), tags = c("a")))
})

test_that("validate_sample_link errors on missing required cols", {
  sl <- data.frame(tag = "a", sample_id = "s1")
  expect_error(validate_sample_link(sl, tags = c("a")))
})

test_that("validate_sample_link passes when tag matches", {
  sl <- data.frame(tag = "exp1", sample_id = "s1", donor_id = "d1")
  expect_no_error(validate_sample_link(sl, tags = c("exp1")))
})

test_that("validate_sample_link errors on unknown tag", {
  sl <- data.frame(tag = "ghost", sample_id = "s1", donor_id = "d1")
  expect_error(validate_sample_link(sl, tags = c("exp1")))
})

# ---- coerce_continuous_col --------------------------------------------

test_that("coerce_continuous_col returns numeric vector unchanged", {
  out <- coerce_continuous_col(c(1, 2, 3), col_name = "age")
  expect_equal(out, c(1, 2, 3))
})

test_that("coerce_continuous_col coerces character-numeric", {
  out <- coerce_continuous_col(c("1", "2", "3"), col_name = "age")
  expect_equal(out, c(1, 2, 3))
})

test_that("coerce_continuous_col errors on all-non-numeric character", {
  expect_error(coerce_continuous_col(c("a", "b"), col_name = "x"))
})

# ---- coerce_to_continuous (matrix-shaped) -----------------------------

test_that("coerce_to_continuous handles raw_count matrix", {
  mat <- matrix(as.integer(rpois(20, 50)), nrow = 5)
  out <- coerce_to_continuous(mat, assay_type = "raw_count")
  expect_true(is.matrix(out))
})

test_that("coerce_to_continuous handles intensity matrix", {
  mat <- matrix(rnorm(20, mean = 10), nrow = 5)
  out <- coerce_to_continuous(mat, assay_type = "intensity")
  expect_true(is.matrix(out))
})

# ---- schema_is_supported ----------------------------------------------

test_that("schema_is_supported returns TRUE for current major version", {
  expect_true(schema_is_supported(OMP_SCHEMA_VERSION))
})

test_that("schema_is_supported returns FALSE for empty input", {
  expect_false(schema_is_supported(""))
})

test_that("schema_is_supported returns FALSE for far-future major", {
  expect_false(schema_is_supported("999.0.0"))
})

# ---- make_unique_labels -----------------------------------------------

test_that("make_unique_labels deduplicates a vector", {
  out <- make_unique_labels(c("a", "a", "b", "b", "c"))
  expect_equal(length(unique(out)), 5)
})

test_that("make_unique_labels preserves unique input", {
  out <- make_unique_labels(c("a", "b", "c"))
  expect_equal(length(unique(out)), 3)
})

test_that("make_unique_labels handles NA", {
  out <- make_unique_labels(c(NA, "a", NA))
  expect_equal(length(out), 3)
})

# ---- normalize_organism / normalize_enrich_database ------------------

test_that("normalize_organism returns standard code for Hs", {
  expect_true(is.character(normalize_organism("Hs")))
})

test_that("normalize_organism handles 'human' alias", {
  out <- normalize_organism("human")
  expect_true(is.character(out) && nchar(out) > 0)
})

test_that("normalize_enrich_database returns input for hallmark", {
  out <- normalize_enrich_database("hallmark")
  expect_true(is.character(out))
})

# ---- compute_numeric_fraction -----------------------------------------

test_that("compute_numeric_fraction is 1 for all-numeric data.frame", {
  expect_equal(compute_numeric_fraction(data.frame(a = 1:3, b = 4:6)), 1)
})

test_that("compute_numeric_fraction is 0 for all-character data.frame", {
  expect_equal(
    compute_numeric_fraction(data.frame(a = c("x", "y", "z"),
                                        stringsAsFactors = FALSE)),
    0
  )
})

test_that("compute_numeric_fraction is partial for mixed data.frame", {
  out <- compute_numeric_fraction(
    data.frame(a = 1:3, b = c("x", "y", "z"), stringsAsFactors = FALSE)
  )
  expect_true(out > 0 && out < 1)
})

test_that("compute_numeric_fraction handles zero-column df", {
  expect_equal(compute_numeric_fraction(data.frame()), 0)
})

# ---- looks_like_feature_labels ----------------------------------------

test_that("looks_like_feature_labels recognizes gene-like strings", {
  out <- looks_like_feature_labels(c("ACTB", "GAPDH", "MYC"))
  expect_true(is.logical(out) && length(out) == 1)
})

test_that("looks_like_feature_labels rejects pure numbers", {
  out <- looks_like_feature_labels(c("1", "2", "3"))
  expect_true(is.logical(out) && length(out) == 1)
})
