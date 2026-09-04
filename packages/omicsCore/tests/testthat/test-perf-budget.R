# Performance budgets, at the size of the data this app is used on.
#
# The follicle RNA-seq workbook is 63,241 genes by 258 samples. Every
# step below has been fast enough on it at some point, and nothing
# would say if one stopped being: the roadmap's Phase 4 latency target
# has no test, and a per-feature loop that quietly replaces a vectorised
# one only shows up when a user waits.
#
# Gated, because generating the matrix and writing it as text takes a
# minute. Run with
#
#   OMICSCORE_PERF_TESTS=1 Rscript -e 'testthat::test_file("tests/testthat/test-perf-budget.R")'
#
# The ceilings are loose -- about three times what a 2023 laptop needs
# -- so they catch a change of algorithm, not a busy machine.

skip_unless_perf <- function() {
  skip_if(!nzchar(Sys.getenv("OMICSCORE_PERF_TESTS", "")),
          "set OMICSCORE_PERF_TESTS=1 to run the performance budget")
}

# Runs `code`, asserts it finished inside `seconds`, and returns what it
# produced so the caller can check the result as well as the clock.
budget <- function(seconds, code) {
  started <- Sys.time()
  value <- force(code)
  elapsed <- as.numeric(Sys.time() - started, units = "secs")
  expect_lt(elapsed, seconds)
  value
}

big_counts <- local({
  cache <- NULL
  function(n_genes = 63241L, n_samples = 258L) {
    if (!is.null(cache)) return(cache)
    set.seed(1)
    mat <- matrix(stats::rpois(n_genes * n_samples, 30L), n_genes, n_samples)
    # A fifth of the genes never vary, as in the real file
    mat[seq_len(n_genes %/% 5L), ] <- 0L
    dimnames(mat) <- list(sprintf("ENSG%011d", seq_len(n_genes)),
                          sprintf("S%03d", seq_len(n_samples)))
    cache <<- mat
    mat
  }
})

test_that("a 63k x 258 tab-separated counts file imports within a minute", {
  skip_unless_perf()
  mat <- big_counts()
  path <- tempfile(fileext = ".xls")   # the pipeline's name for a TSV
  on.exit(unlink(path), add = TRUE)
  df <- data.frame(gene_id = rownames(mat), mat, check.names = FALSE)
  utils::write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)

  got <- budget(60, read_omics(path, omics_type = "rnaseq", assay_type = "raw_count"))
  expect_identical(dim(got$input$expr_mat), dim(mat))
})

test_that("PCA over 63k genes, a fifth of them constant, takes seconds not minutes", {
  skip_unless_perf()
  mat <- big_counts()
  budget(1, row_variance(mat))
  pca <- budget(20, pca_over_samples(log2(mat + 1)))
  expect_identical(attr(pca, "n_dropped"), nrow(mat) %/% 5L)
})

test_that("re-masking a 63k-row result at a new threshold is instant", {
  skip_unless_perf()
  n <- 63241L
  df <- data.frame(
    feature_id = sprintf("g%d", seq_len(n)), feature_symbol = sprintf("g%d", seq_len(n)),
    feature_type = "gene", omics_type = "rnaseq", method = "deseq2",
    analysis_type = "group", comparison = "B_vs_A",
    effect = stats::rnorm(n), effect_type = "log2FC",
    statistic = stats::rnorm(n), statistic_type = "wald",
    p_value = stats::runif(n), adj_p_value = stats::runif(n),
    direction = "ns", base_mean = stats::runif(n, 1, 1000),
    model_fit = NA_real_, is_significant = NA,
    stringsAsFactors = FALSE
  )
  # This is what a slider drag costs the user; the view debounces to
  # one call per 250 ms, so one call must be well inside that.
  budget(0.25, filter_diff_results(df, p_cutoff = 0.05, effect_cutoff = 0.5))
  budget(0.5, make_ranked_features(df))
})

test_that("limma over 20k features and 258 samples fits inside the progress bar", {
  skip_unless_perf()
  skip_if_not_installed("limma")
  set.seed(2)
  n_feat <- 20000L; n_samp <- 258L
  mat <- matrix(stats::rnorm(n_feat * n_samp, 20, 2), n_feat, n_samp,
                dimnames = list(sprintf("P%05d", seq_len(n_feat)),
                                sprintf("S%03d", seq_len(n_samp))))
  meta <- data.frame(group = rep(c("G1", "G2"), length.out = n_samp),
                     age = stats::runif(n_samp, 20, 80),
                     row.names = colnames(mat))
  inp <- omics_input(mat, meta, data.frame(feature_id = rownames(mat)),
                     omics_type = "proteomics", assay_type = "normalized_intensity")
  b <- budget(30, run_diff(inp, method = "limma", analysis_type = "group",
                           group_col = "group", control_group = "G1",
                           case_group = "G2", covariates = "age"))
  expect_identical(nrow(b$results$diff_result_df), n_feat)
})
