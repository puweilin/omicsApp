# GSVA, install_optional, and remaining utility robustness tests.
#
# The GSVA tests use the realistic fixture (helper-realistic.R): real
# gene symbols, so the Hallmark sets actually overlap the matrix. They
# used to wrap run_gsva() in tryCatch and skip on any error, which turned
# a real defect -- the symbol lookup dropped every feature -- into
# "GSVA bundle not buildable in this env", on every machine, on every
# run. They also read a result slot that has never existed. A skip
# cannot be wrong; that is exactly why it must not stand in for a
# failure.

# ---- run_gsva ----------------------------------------------------------

test_that("run_gsva returns a bundle on valid input", {
  skip_if_not_installed("GSVA")
  skip_if_not_installed("msigdbr")
  b <- run_gsva(realistic_input(), database = "hallmark")
  expect_true(is_analysis_bundle(b))
  expect_identical(b$analysis_name, "run_gsva")
  expect_identical(b$params$method, "gsva")
})

test_that("run_gsva produces a gene-set by sample score matrix", {
  skip_if_not_installed("GSVA")
  skip_if_not_installed("msigdbr")
  inp <- realistic_input()
  b <- run_gsva(inp, database = "hallmark")
  m <- b$results$gsva_matrix
  expect_true(is.matrix(m))
  expect_identical(colnames(m), colnames(inp$expr_mat))
  expect_true(all(grepl("^HALLMARK_", rownames(m))))
  expect_true(all(is.finite(m)))
  expect_true(all(abs(m) <= 1))
  expect_type(b$results$gsva_gene_sets, "list")
})

test_that("run_gsva supports ssgsea method", {
  skip_if_not_installed("GSVA")
  b <- run_gsva(realistic_input(), gene_sets = REAL_GENE_SETS,
                min_size = 5L, method = "ssgsea")
  expect_true(is_analysis_bundle(b))
  expect_identical(b$params$method, "ssgsea")
  expect_setequal(rownames(b$results$gsva_matrix), names(REAL_GENE_SETS))
})

test_that("run_gsva respects min_size and max_size", {
  skip_if_not_installed("GSVA")
  # Every supplied set has 40 members: a floor of 40 keeps them all, a
  # ceiling of 39 leaves nothing to score.
  keep <- run_gsva(realistic_input(), gene_sets = REAL_GENE_SETS,
                   min_size = 40L, max_size = 40L)
  expect_setequal(rownames(keep$results$gsva_matrix), names(REAL_GENE_SETS))
  expect_error(
    run_gsva(realistic_input(), gene_sets = REAL_GENE_SETS,
             min_size = 5L, max_size = 39L),
    "[Nn]o gene set"
  )
})

test_that("run_gsva errors on bogus database", {
  skip_if_not_installed("GSVA")
  skip_if_not_installed("msigdbr")
  expect_error(run_gsva(realistic_input(), database = "ghost_db"))
})

# ---- install_optional groups ------------------------------------------

test_that("OPTIONAL_GROUPS contains rnaseq group", {
  expect_true("rnaseq" %in% names(OPTIONAL_GROUPS))
})

test_that("OPTIONAL_GROUPS contains proteomics group", {
  expect_true("proteomics" %in% names(OPTIONAL_GROUPS))
})

test_that("OPTIONAL_GROUPS contains enrichment group", {
  expect_true("enrichment" %in% names(OPTIONAL_GROUPS))
})

test_that("OPTIONAL_GROUPS contains imputation group", {
  expect_true("imputation" %in% names(OPTIONAL_GROUPS))
})

test_that("OPTIONAL_GROUPS rnaseq lists DESeq2 and edgeR", {
  expect_true(all(c("DESeq2", "edgeR") %in% OPTIONAL_GROUPS$rnaseq))
})

test_that("OPTIONAL_GROUPS enrichment lists fgsea", {
  expect_true("fgsea" %in% OPTIONAL_GROUPS$enrichment)
})

test_that("OPTIONAL_GROUPS imputation lists missForest", {
  expect_true("missForest" %in% OPTIONAL_GROUPS$imputation)
})

test_that("check_install with no missing groups returns invisibly", {
  # base packages should always be installed
  expect_no_error(check_install("rnaseq"))
})

# ---- BIOC_PACKAGES, CANONICAL_DATABASES, SUPPORTED_* constants --------

test_that("BIOC_PACKAGES is a character vector", {
  expect_true(is.character(BIOC_PACKAGES) && length(BIOC_PACKAGES) > 0)
})

test_that("CANONICAL_DATABASES includes core databases", {
  expect_true("KEGG" %in% CANONICAL_DATABASES &&
                "Reactome" %in% CANONICAL_DATABASES)
})

test_that("SUPPORTED_OMICS_TYPES includes proteomics and rnaseq", {
  expect_true(all(c("proteomics", "rnaseq") %in% SUPPORTED_OMICS_TYPES))
})

test_that("SUPPORTED_DIFF_METHODS includes ttest", {
  expect_true("ttest" %in% SUPPORTED_DIFF_METHODS)
})

test_that("SUPPORTED_DIFF_ANALYSIS_TYPES includes group", {
  expect_true("group" %in% SUPPORTED_DIFF_ANALYSIS_TYPES)
})

test_that("SUPPORTED_ENRICH_TYPES includes ora and gsea", {
  expect_true(all(c("ora", "gsea") %in% SUPPORTED_ENRICH_TYPES))
})

test_that("SUPPORTED_INTEGRATION_METHODS includes concordance", {
  expect_true("concordance" %in% SUPPORTED_INTEGRATION_METHODS)
})
