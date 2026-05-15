# Tests for package constants.
# Values must match constants.R exactly.

test_that("SUPPORTED_OMICS_TYPES contains expected types", {
  expect_type(SUPPORTED_OMICS_TYPES, "character")
  expect_true("proteomics" %in% SUPPORTED_OMICS_TYPES)
  expect_true("rnaseq" %in% SUPPORTED_OMICS_TYPES)
  expect_length(SUPPORTED_OMICS_TYPES, 2)
})

test_that("SUPPORTED_DIFF_ANALYSIS_TYPES contains expected values", {
  expect_type(SUPPORTED_DIFF_ANALYSIS_TYPES, "character")
  expect_true("group" %in% SUPPORTED_DIFF_ANALYSIS_TYPES)
  expect_true("continuous" %in% SUPPORTED_DIFF_ANALYSIS_TYPES)
  expect_true("anova" %in% SUPPORTED_DIFF_ANALYSIS_TYPES)
  expect_length(SUPPORTED_DIFF_ANALYSIS_TYPES, 3)
})

test_that("SUPPORTED_ENRICH_PREFERENCE is valid", {
  expect_type(SUPPORTED_ENRICH_PREFERENCE, "character")
  expect_length(SUPPORTED_ENRICH_PREFERENCE, 3)
  expect_true("adjusted" %in% SUPPORTED_ENRICH_PREFERENCE)
})

test_that("SUPPORTED_PREFERENCE is a character vector", {
  expect_type(SUPPORTED_PREFERENCE, "character")
  expect_length(SUPPORTED_PREFERENCE, 2)
})

test_that("DIFF_RESULT_REQUIRED_COLS is a character vector", {
  expect_type(DIFF_RESULT_REQUIRED_COLS, "character")
  expect_true("feature_id" %in% DIFF_RESULT_REQUIRED_COLS)
  expect_true("effect" %in% DIFF_RESULT_REQUIRED_COLS)
  expect_true("adj_p_value" %in% DIFF_RESULT_REQUIRED_COLS)
})

test_that("ENRICH_RESULT_REQUIRED_COLS is a character vector", {
  expect_type(ENRICH_RESULT_REQUIRED_COLS, "character")
  expect_true("pathway_name" %in% ENRICH_RESULT_REQUIRED_COLS)
  expect_true("adj_p_value" %in% ENRICH_RESULT_REQUIRED_COLS)
})

test_that("INTEGRATION_RESULT_REQUIRED_COLS is a character vector", {
  expect_type(INTEGRATION_RESULT_REQUIRED_COLS, "character")
})

test_that("IMPORT_REPORT_ROLES is a named character vector", {
  expect_type(IMPORT_REPORT_ROLES, "character")
  expect_true(length(IMPORT_REPORT_ROLES) > 0)
  expect_true(all(nzchar(names(IMPORT_REPORT_ROLES))))
})

test_that("DEFAULT_CACHE_MAX_AGE_HOURS is a positive number", {
  expect_type(DEFAULT_CACHE_MAX_AGE_HOURS, "double")
  expect_gt(DEFAULT_CACHE_MAX_AGE_HOURS, 0)
})

test_that("CANONICAL_DATABASES is a character vector", {
  expect_type(CANONICAL_DATABASES, "character")
  expect_true(length(CANONICAL_DATABASES) > 0)
})
