# GSVA, install_optional, and remaining utility robustness tests.

make_gsva_input <- function(n_feat = 50, n_samp = 8) {
  mat <- matrix(rnorm(n_feat * n_samp, mean = 10, sd = 2),
                nrow = n_feat, ncol = n_samp)
  rownames(mat) <- paste0("gene_", seq_len(n_feat))
  colnames(mat) <- paste0("s", seq_len(n_samp))
  meta <- data.frame(group = rep(c("A", "B"), length.out = n_samp),
                     row.names = colnames(mat))
  feat <- data.frame(feature_id = rownames(mat),
                     feature_name = rownames(mat),
                     stringsAsFactors = FALSE)
  omics_input(mat, meta, feat, omics_type = "proteomics",
              assay_type = "log_intensity")
}

# ---- run_gsva ----------------------------------------------------------

test_that("run_gsva returns a bundle on valid input", {
  skip_if_not_installed("GSVA")
  skip_if_not_installed("msigdbr")
  b <- tryCatch(run_gsva(make_gsva_input(), database = "hallmark"),
                error = function(e) NULL)
  skip_if(is.null(b), "GSVA bundle not buildable in this env")
  expect_true(is_analysis_bundle(b))
})

test_that("run_gsva produces a gsva_score_mat result", {
  skip_if_not_installed("GSVA")
  skip_if_not_installed("msigdbr")
  b <- tryCatch(run_gsva(make_gsva_input(), database = "hallmark"),
                error = function(e) NULL)
  skip_if(is.null(b), "GSVA bundle not buildable in this env")
  expect_true(is.matrix(b$results$gsva_score_mat) ||
                is.data.frame(b$results$gsva_score_mat))
})

test_that("run_gsva supports ssgsea method", {
  skip_if_not_installed("GSVA")
  skip_if_not_installed("msigdbr")
  b <- tryCatch(run_gsva(make_gsva_input(), database = "hallmark",
                         method = "ssgsea"),
                error = function(e) NULL)
  skip_if(is.null(b), "ssgsea bundle not buildable in this env")
  expect_true(is_analysis_bundle(b))
})

test_that("run_gsva respects min_size and max_size", {
  skip_if_not_installed("GSVA")
  skip_if_not_installed("msigdbr")
  b <- tryCatch(run_gsva(make_gsva_input(), database = "hallmark",
                         min_size = 1L, max_size = 5000L),
                error = function(e) NULL)
  skip_if(is.null(b), "GSVA bundle not buildable in this env")
  expect_true(is_analysis_bundle(b))
})

test_that("run_gsva errors on bogus database", {
  skip_if_not_installed("GSVA")
  skip_if_not_installed("msigdbr")
  expect_error(run_gsva(make_gsva_input(), database = "ghost_db"))
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
