# Edge-case tests for validation utilities.
# Covers check_required_cols, check_meta_matches_expr,
# check_feature_matches_expr, check_paired_col, and the two
# validate_*_pairing functions.

# ---- check_required_cols -----------------------------------------------

test_that("check_required_cols passes when all columns present", {
  df <- data.frame(a = 1, b = 2, c = 3)
  expect_true(check_required_cols(df, c("a", "b")))
})

test_that("check_required_cols errors when columns missing", {
  df <- data.frame(x = 1, y = 2)
  expect_error(check_required_cols(df, c("x", "z")), "z")
})

test_that("check_required_cols includes object_name in error", {
  df <- data.frame(x = 1)
  expect_error(check_required_cols(df, "missing", object_name = "my_df"),
               "my_df")
})

test_that("check_required_cols returns invisibly", {
  df <- data.frame(a = 1)
  res <- withVisible(check_required_cols(df, "a"))
  expect_false(res$visible)
})

# ---- check_meta_matches_expr -------------------------------------------

test_that("check_meta_matches_expr passes for matching metadata", {
  expr <- matrix(1:6, nrow = 2, dimnames = list(c("g1", "g2"), c("s1", "s2", "s3")))
  meta <- data.frame(group = c("A", "A", "B"), row.names = c("s1", "s2", "s3"))
  expect_true(check_meta_matches_expr(expr, meta))
})

test_that("check_meta_matches_expr errors when expr lacks colnames", {
  expr <- matrix(1:6, nrow = 2, ncol = 3)
  meta <- data.frame(group = letters[1:3], row.names = letters[1:3])
  expect_error(check_meta_matches_expr(expr, meta), "sample column names")
})

test_that("check_meta_matches_expr errors when meta lacks rownames", {
  expr <- matrix(1:6, nrow = 2, dimnames = list(c("g1", "g2"), letters[1:3]))
  meta <- data.frame(group = letters[1:3])
  expect_error(check_meta_matches_expr(expr, meta), "rownames")
})

test_that("check_meta_matches_expr errors when samples not in meta", {
  expr <- matrix(1:6, nrow = 2, dimnames = list(c("g1", "g2"), c("s1", "s2", "s3")))
  meta <- data.frame(group = c("A", "A"), row.names = c("s1", "s2"))
  expect_error(check_meta_matches_expr(expr, meta), "Not all samples")
})

test_that("check_meta_matches_expr accepts all samples present in meta", {
  expr <- matrix(1:9, nrow = 3, dimnames = list(c("g1", "g2", "g3"), c("s1", "s2", "s3")))
  meta <- data.frame(group = c("A", "A", "B", "C"), row.names = c("s1", "s2", "s3", "s4"))
  expect_true(check_meta_matches_expr(expr, meta))
})

# ---- check_feature_matches_expr ----------------------------------------

test_that("check_feature_matches_expr passes for matching feature data", {
  expr <- matrix(1:6, nrow = 2, dimnames = list(c("f1", "f2"), c("s1", "s2", "s3")))
  feat <- data.frame(feature_id = c("f1", "f2"), name = c("a", "b"))
  expect_true(check_feature_matches_expr(expr, feat))
})

test_that("check_feature_matches_expr errors when expr lacks rownames", {
  expr <- matrix(1:6, nrow = 2, ncol = 3)
  feat <- data.frame(feature_id = c("f1", "f2"))
  expect_error(check_feature_matches_expr(expr, feat), "feature row names")
})

test_that("check_feature_matches_expr errors when feature_id_col missing", {
  expr <- matrix(1:6, nrow = 2, dimnames = list(c("f1", "f2"), c("s1", "s2", "s3")))
  feat <- data.frame(id = c("f1", "f2"))
  expect_error(check_feature_matches_expr(expr, feat, feature_id_col = "feature_id"),
               "feature_id")
})

test_that("check_feature_matches_expr errors when features not in reference", {
  expr <- matrix(1:6, nrow = 2, dimnames = list(c("f1", "f99"), c("s1", "s2", "s3")))
  feat <- data.frame(feature_id = c("f1", "f2"))
  expect_error(check_feature_matches_expr(expr, feat), "Not all features")
})

test_that("check_feature_matches_expr uses custom feature_id_col", {
  expr <- matrix(1:6, nrow = 2, dimnames = list(c("u1", "u2"), c("s1", "s2", "s3")))
  feat <- data.frame(uniprot = c("u1", "u2"), symbol = c("A", "B"))
  expect_true(check_feature_matches_expr(expr, feat, feature_id_col = "uniprot"))
})

# ---- check_paired_col --------------------------------------------------

test_that("check_paired_col passes when paired_col exists and has values", {
  meta <- data.frame(pair_id = c("P1", "P1", "P2"), row.names = c("s1", "s2", "s3"))
  expect_true(check_paired_col(meta, "pair_id"))
})

test_that("check_paired_col does nothing when paired_col is NULL", {
  meta <- data.frame(x = 1, row.names = "s1")
  expect_true(check_paired_col(meta, NULL))
})

test_that("check_paired_col errors when column missing", {
  meta <- data.frame(x = 1, row.names = "s1")
  expect_error(check_paired_col(meta, "nonexistent"), "not found")
})

test_that("check_paired_col errors when column is all NA", {
  meta <- data.frame(pair_id = c(NA, NA), row.names = c("s1", "s2"))
  expect_error(check_paired_col(meta, "pair_id"), "only missing values")
})

test_that("check_paired_col includes object_name in errors", {
  meta <- data.frame(x = 1, row.names = "s1")
  expect_error(check_paired_col(meta, "y", object_name = "custom"),
               "custom")
})

# ---- validate_two_group_pairing ----------------------------------------

test_that("validate_two_group_pairing passes for valid paired design", {
  meta <- data.frame(
    group = c("ctrl", "case", "ctrl", "case"),
    pair  = c("P1", "P1", "P2", "P2"),
    row.names = c("s1", "s2", "s3", "s4")
  )
  expect_true(validate_two_group_pairing(meta, "group", "pair", "ctrl", "case"))
})

test_that("validate_two_group_pairing no-ops when paired_col is NULL", {
  meta <- data.frame(group = c("ctrl", "case"), row.names = c("s1", "s2"))
  expect_true(validate_two_group_pairing(meta, "group", NULL, "ctrl", "case"))
})

test_that("validate_two_group_pairing errors when a pair is incomplete", {
  meta <- data.frame(
    group = c("ctrl", "case", "ctrl"),
    pair  = c("P1", "P1", "P2"),
    row.names = c("s1", "s2", "s3")
  )
  expect_error(
    validate_two_group_pairing(meta, "group", "pair", "ctrl", "case"),
    "Invalid paired design"
  )
})

test_that("validate_two_group_pairing errors when a pair has duplicate groups", {
  meta <- data.frame(
    group = c("ctrl", "ctrl", "case", "case"),
    pair  = c("P1", "P1", "P2", "P2"),
    row.names = c("s1", "s2", "s3", "s4")
  )
  # P1 has 2 ctrl, no case → not {ctrl, case}
  # P2 has 2 case, no ctrl → not {ctrl, case}
  expect_error(
    validate_two_group_pairing(meta, "group", "pair", "ctrl", "case"),
    "exactly one"
  )
})

test_that("validate_two_group_pairing handles incomplete pairs from NA filtering", {
  # Row 3 has NA group → dropped from pair checking; P2 is left with only 1
  # sample (ctrl), making it incomplete.
  meta <- data.frame(
    group = c("ctrl", "case", "case", "ctrl"),
    pair  = c("P1", "P1", "P2", "P2"),
    row.names = c("s1", "s2", "s3", "s4"),
    stringsAsFactors = FALSE
  )
  # P2 has case + ctrl → should NOT error
  expect_true(validate_two_group_pairing(meta, "group", "pair", "ctrl", "case"))
})

# ---- validate_continuous_pairing ---------------------------------------

test_that("validate_continuous_pairing passes when every pair has >=2 samples", {
  meta <- data.frame(
    pair = c("P1", "P1", "P2", "P2", "P2"),
    row.names = c("s1", "s2", "s3", "s4", "s5")
  )
  expect_true(validate_continuous_pairing(meta, "pair"))
})

test_that("validate_continuous_pairing no-ops when paired_col is NULL", {
  meta <- data.frame(x = 1, row.names = "s1")
  expect_true(validate_continuous_pairing(meta, NULL))
})

test_that("validate_continuous_pairing errors when a pair has <2 samples", {
  meta <- data.frame(
    pair = c("P1", "P1", "P2"),
    row.names = c("s1", "s2", "s3")
  )
  expect_error(
    validate_continuous_pairing(meta, "pair"),
    "at least 2"
  )
})

test_that("validate_continuous_pairing errors when paired_col missing", {
  meta <- data.frame(x = 1, row.names = "s1")
  expect_error(validate_continuous_pairing(meta, "missing_col"), "required columns")
})

# ---- validate_omics_input -----------------------------------------------

make_valid_input <- function() {
  mat <- matrix(1:12, nrow = 3, dimnames = list(paste0("g", 1:3), paste0("s", 1:4)))
  meta <- data.frame(group = 1:4, row.names = paste0("s", 1:4))
  feat <- data.frame(feature_id = paste0("g", 1:3))
  omics_input(mat, meta, feat, omics_type = "proteomics", assay_type = "normalized_intensity")
}

test_that("validate_omics_input returns TRUE for valid input", {
  inp <- make_valid_input()
  expect_true(validate_omics_input(inp))
})

test_that("validate_omics_input errors for non-omics_input objects", {
  expect_error(validate_omics_input(list()), "omics_input")
})

test_that("validate_omics_input errors when expr_mat is NULL", {
  inp <- make_valid_input()
  inp$expr_mat <- NULL
  expect_error(validate_omics_input(inp), "expr_mat")
})

test_that("validate_omics_input errors when meta_df is NULL", {
  inp <- make_valid_input()
  inp$meta_df <- NULL
  expect_error(validate_omics_input(inp), "meta_df")
})

test_that("validate_omics_input errors when feature_df is NULL", {
  inp <- make_valid_input()
  inp$feature_df <- NULL
  expect_error(validate_omics_input(inp), "feature_df")
})

test_that("validate_omics_input errors when omics_type field is missing", {
  inp <- make_valid_input()
  inp$omics_type <- NULL
  expect_error(validate_omics_input(inp), "omics_type")
})

test_that("validate_omics_input errors when assay_type field is missing", {
  inp <- make_valid_input()
  inp$assay_type <- NULL
  expect_error(validate_omics_input(inp), "assay_type")
})
