# Which engines a caller may choose is a correctness question, not a
# convenience one: DESeq2 handed continuous intensities rounds them to
# integers and reports p-values for a model the data never fitted, and
# this package's limma backend has no voom step to make raw counts
# legitimate. Neither mistake errors, so nothing downstream can catch it.

make_input <- function(omics_type, assay_type) {
  mat <- matrix(as.numeric(1:24), nrow = 6,
                dimnames = list(paste0("g", 1:6), paste0("s", 1:4)))
  meta <- data.frame(group = c("A", "A", "B", "B"),
                     row.names = paste0("s", 1:4))
  feat <- data.frame(feature_id = paste0("g", 1:6))
  omics_input(mat, meta, feat, omics_type = omics_type,
              assay_type = assay_type)
}

test_that("raw counts get the count-native engines only", {
  m <- applicable_diff_methods(make_input("rnaseq", "raw_count"))
  expect_true(all(c("deseq2", "edger") %in% m))
  # limma without voom treats counts as continuous.
  expect_false(any(c("limma", "ttest", "lm") %in% m))
})

test_that("continuous data never offers a count model", {
  for (inp in list(make_input("proteomics", "intensity"),
                   make_input("proteomics", "normalized_intensity"),
                   make_input("rnaseq", "tpm"),
                   make_input("rnaseq", "fpkm"))) {
    m <- applicable_diff_methods(inp)
    expect_false(any(c("deseq2", "edger") %in% m),
                 info = paste(inp$omics_type, inp$assay_type))
    expect_true("limma" %in% m)
  }
})

test_that("rnaseq only counts as counts when the assay says so", {
  # omics_type alone is not enough: TPM is rnaseq and is continuous.
  expect_false("deseq2" %in% applicable_diff_methods(make_input("rnaseq", "tpm")))
  expect_true("deseq2" %in%
                applicable_diff_methods(make_input("rnaseq", "raw_count")))
})

test_that("auto is always offered", {
  expect_true("auto" %in% applicable_diff_methods(make_input("rnaseq", "raw_count")))
  expect_true("auto" %in% applicable_diff_methods(make_input("proteomics", "intensity")))
})

test_that("a continuous predictor narrows to the regression backends", {
  m <- applicable_diff_methods(make_input("proteomics", "intensity"),
                               analysis_type = "continuous")
  expect_true(all(c("limma", "lm") %in% m))
  # A two-sample t-test has no continuous predictor to regress on.
  expect_false("ttest" %in% m)
})

test_that("every offered method is one run_diff() accepts", {
  # A gate that offers a name the dispatcher rejects would trade a
  # silent wrong answer for a loud crash -- better, but still a bug.
  for (inp in list(make_input("rnaseq", "raw_count"),
                   make_input("proteomics", "intensity"))) {
    expect_true(all(applicable_diff_methods(inp) %in% SUPPORTED_DIFF_METHODS))
  }
})

test_that("auto never resolves to a method the gate would have hidden", {
  skip_if_not_installed("limma")
  inp <- make_input("proteomics", "intensity")
  chosen <- auto_select_diff_method(inp, "group")
  expect_true(chosen %in% applicable_diff_methods(inp))
})
