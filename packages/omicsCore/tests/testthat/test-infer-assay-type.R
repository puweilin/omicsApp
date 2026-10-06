# RNA-seq used to be called raw counts whatever it held, so a TPM or log
# table was offered DESeq2 and edgeR, which round it and report p-values
# for a model it never fitted. The scale is now read from the values.

rna_values <- function(kind, n = 2000, s = 6) {
  set.seed(21)
  counts <- matrix(stats::rnbinom(n * s, mu = exp(stats::rnorm(n, 5, 1.5)), size = 5), n,
                   dimnames = list(sprintf("ENSG%011d", seq_len(n)), paste0("s", seq_len(s))))
  len <- stats::runif(n, 500, 4000)
  switch(kind,
    counts = counts,
    tpm = {
      rate <- counts / len
      sweep(rate, 2L, colSums(rate), "/") * 1e6
    },
    fpkm = sweep(counts / (len / 1e3), 2L, colSums(counts) / 1e6, "/"),
    logcpm = log2(sweep(counts + 0.5, 2L, colSums(counts) + 1, "/") * 1e6),
    vst = log2(counts + 1) + stats::rnorm(n * s, 0, 0.01),
    estimated = counts + stats::runif(n * s) * (counts > 0)
  )
}

test_that("non-negative whole numbers are read counts", {
  expect_identical(infer_assay_type(rna_values("counts"), "rnaseq"), "raw_count")
  g <- infer_assay_type(rna_values("counts"), "rnaseq", explain = TRUE)
  expect_identical(g$assay_type, "raw_count")
  expect_match(g$reason, "whole number")
})

test_that("fractional values adding up to a million per sample are TPM", {
  tpm <- rna_values("tpm")
  expect_equal(unname(colSums(tpm)), rep(1e6, 6))
  expect_identical(infer_assay_type(tpm, "rnaseq"), "tpm")
  # Still TPM after the low genes were filtered out.
  kept <- tpm[rowMeans(tpm) > stats::quantile(rowMeans(tpm), 0.25), ]
  expect_lt(max(colSums(kept)), 1e6)
  expect_identical(infer_assay_type(kept, "rnaseq"), "tpm")
})

test_that("negative values, or fractional values below 30, are log-scale", {
  logcpm <- rna_values("logcpm")
  expect_true(any(logcpm < 0))
  expect_identical(infer_assay_type(logcpm, "rnaseq"), "logcpm")
  vst <- abs(rna_values("vst"))
  expect_lt(max(vst), 30)
  expect_identical(infer_assay_type(vst, "rnaseq"), "logcpm")
  expect_match(infer_assay_type(vst, "rnaseq", explain = TRUE)$reason, "never reach 30")
})

test_that("other fractional values are FPKM-like, and estimated counts stay counts", {
  expect_identical(infer_assay_type(rna_values("fpkm"), "rnaseq"), "fpkm")
  # Salmon / RSEM estimated reads: fractional, library-sized totals.
  est <- rna_values("estimated") * 20
  expect_gt(stats::median(colSums(est)), 2e6)
  expect_identical(infer_assay_type(est, "rnaseq"), "raw_count")
})

test_that("proteomics is still told apart by magnitude", {
  set.seed(2)
  expect_identical(infer_assay_type(matrix(2^stats::rnorm(40, 20, 2), 10), "proteomics"),
                   "raw_intensity")
  expect_identical(infer_assay_type(matrix(stats::rnorm(40, 20, 2), 10), "proteomics"),
                   "normalized_intensity")
  expect_identical(infer_assay_type(matrix(1, 2, 2), "metabolomics"), NA_character_)
})

test_that("only a counts layer is offered DESeq2 and edgeR", {
  mk <- function(kind, assay) {
    m <- rna_values(kind)
    suppressWarnings(omics_input(m, data.frame(group = rep(c("A", "B"), each = 3),
                                               row.names = colnames(m)),
                                 data.frame(feature_id = rownames(m)),
                                 omics_type = "rnaseq",
                                 assay_type = infer_assay_type(m, "rnaseq")))
  }
  counts <- mk("counts")
  expect_true(all(c("deseq2", "edger") %in% applicable_diff_methods(counts)))
  for (kind in c("tpm", "fpkm", "logcpm")) {
    inp <- mk(kind)
    methods <- applicable_diff_methods(inp)
    expect_false(any(c("deseq2", "edger") %in% methods), info = kind)
    expect_true("limma" %in% methods, info = kind)
  }
})
