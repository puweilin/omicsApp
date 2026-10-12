# The standardizers were eight copies of one join-and-transmute (plus
# two global-test variants), and are now mappings handed to one shared
# function. The fixture holds engine tables taken from every backend and
# analysis type -- limma, DESeq2, edgeR, t-test and lm; two-group,
# several groups, weighted contrasts, covariates, pairing, continuous,
# spline and global tests -- trimmed to a few rows each (with untested
# rows kept), together with the standardized table the code before the
# refactor (commit cf9823b) made from each. Also: engine tables missing
# their optional columns, a feature table with no symbol or type column,
# and a table with no rows. The output must not have changed by a bit.

test_that("every standardizer gives exactly the table it gave before the refactor", {
  cases <- readRDS(test_path("fixtures", "diff-standardize-snapshot.rds"))
  expect_gt(length(cases), 25L)
  funs <- unique(vapply(cases, `[[`, character(1), "fun"))
  expect_true(all(c("standardize_limma_group_results", "standardize_limma_continuous_results",
                    "standardize_deseq2_group_results", "standardize_deseq2_continuous_results",
                    "standardize_edger_group_results", "standardize_ttest_group_results",
                    "standardize_lm_group_results", "standardize_lm_continuous_results",
                    "standardize_global_test_results") %in% funs))
  for (cs in cases) {
    got <- do.call(get(cs$fun, envir = asNamespace("omicsCore")), cs$args)
    expect_identical(got, cs$expected, label = cs$label)
  }
})

test_that("a standardizer still names the engine column that is missing", {
  raw <- data.frame(feature_id = "g1", logFC = 1, P.Value = 0.01)
  expect_error(
    standardize_limma_group_results(raw, data.frame(feature_id = "g1"), "B_vs_A"),
    "Missing required columns in raw_df: adj.P.Val", fixed = TRUE)
})
