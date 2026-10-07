# `run_diff(prefilter = )`: whether DESeq2 or edgeR set aside genes with
# too few counts before fitting. NULL keeps each engine's default (edgeR
# filters with filterByExpr(), DESeq2 does not); TRUE and FALSE turn the
# filter on or off for either. Genes set aside stay in the table as
# untested rows, and the bundle says how many.
#
# DESeq2's filter is off by default on the evidence: on a simulated
# 30,000 x 60 set (two groups of 30) its recommended rule removed 13,718
# genes, 12,383 of which had an adjusted p-value without it, and lost 386
# of 1,282 genes called at adjusted p < 0.05 -- a speed-up bought with
# results. See deseq2_prefilter() in diff-deseq2.R.

pf_counts <- function(per_group = 4L, seed = 3) {
  set.seed(seed)
  n <- 2L * per_group
  meta <- data.frame(group = rep(c("A", "B"), each = per_group),
                     dose = seq(0, 7, length.out = n),
                     row.names = sprintf("S%02d", seq_len(n)))
  mu <- c(exp(stats::runif(150, log(30), log(3000))), stats::runif(50, 0.05, 1.5))
  fc <- rep(1, 200)
  fc[1:15] <- 4
  m <- outer(mu, rep(1, n))
  m[, meta$group == "B"] <- m[, meta$group == "B"] * fc
  counts <- matrix(stats::rnbinom(200 * n, mu = m, size = 10), 200, n,
                   dimnames = list(sprintf("G%03d", 1:200), rownames(meta)))
  counts[200, ] <- 0L   # one gene with no reads at all
  omics_input(counts, meta, data.frame(feature_id = rownames(counts)),
              omics_type = "rnaseq", assay_type = "raw_count")
}

pf_group <- function(inp, method, ...) {
  suppressMessages(run_diff(inp, method = method, group_col = "group",
                            control_group = "A", case_group = "B", ...))
}

test_that("leaving prefilter unset changes nothing, and FALSE is DESeq2's default", {
  skip_if_not_installed("DESeq2")
  inp <- pf_counts()
  unset <- pf_group(inp, "deseq2")
  off <- pf_group(inp, "deseq2", prefilter = FALSE)
  expect_identical(off$results$diff_result_df, unset$results$diff_result_df)
  expect_identical(off$results$diff_raw_df, unset$results$diff_raw_df)
  expect_identical(off$warnings, unset$warnings)
  # A call without the option records nothing about it, so a script
  # exported from it is the call it always was.
  expect_false("prefilter" %in% names(unset$params))
  expect_identical(off$params$prefilter, FALSE)
})

test_that("DESeq2 with prefilter = TRUE fits the genes the rule keeps and lists the rest untested", {
  skip_if_not_installed("DESeq2")
  inp <- pf_counts()
  b <- pf_group(inp, "deseq2", prefilter = TRUE)
  df <- b$results$diff_result_df
  expect_identical(df$feature_id, rownames(inp$expr_mat))
  keep <- rowSums(inp$expr_mat >= 10) >= 4   # the smaller group has 4 samples
  expect_gt(sum(!keep), 30L)
  expect_true(all(is.na(df$p_value[!keep])))
  expect_true(all(is.na(df$adj_p_value[!keep])))
  expect_true(all(is.na(df$effect[!keep])))
  expect_true(all(is.na(df$signed_stat[!keep])))
  expect_true(all(df$direction[!keep] == "ns"))
  expect_false(anyNA(df$p_value[keep]))
  expect_identical(b$params$prefilter, TRUE)
  expect_true(any(grepl(sprintf("%d of 200 genes with too few counts to test were set aside before fitting",
                                sum(!keep)), b$warnings, fixed = TRUE)))
  expect_true(any(grepl("fewer than 10 reads in 4 or more samples", b$warnings, fixed = TRUE)))
  # The fitted genes are DESeq2 on those genes alone.
  cd <- data.frame(g = factor(inp$meta_df$group), row.names = rownames(inp$meta_df))
  dds <- suppressMessages(DESeq2::DESeq(
    DESeq2::DESeqDataSetFromMatrix(inp$expr_mat[keep, ], cd, ~ g), quiet = TRUE))
  ref <- as.data.frame(DESeq2::results(dds, contrast = c("g", "B", "A")))
  expect_equal(df$p_value[keep], ref$pvalue, tolerance = 1e-8)
  expect_equal(df$effect[keep], ref$log2FoldChange, tolerance = 1e-8)
  expect_equal(df$adj_p_value[keep], ref$padj, tolerance = 1e-8)
  # The planted changes are still found.
  expect_gte(sum(df$adj_p_value[1:15] < 0.05, na.rm = TRUE), 12L)
})

test_that("the smallest group sets DESeq2's rule, and designs without groups get an equivalent", {
  skip_if_not_installed("DESeq2")
  inp <- pf_counts()
  # Three B samples of four: the rule asks for 3 samples, not 4.
  inp3 <- subset_omics_samples(inp, rownames(inp$meta_df)[-8])
  b <- pf_group(inp3, "deseq2", prefilter = TRUE)
  keep <- unname(rowSums(inp3$expr_mat >= 10) >= 3)
  expect_identical(is.na(b$results$diff_result_df$p_value), !keep)
  expect_true(any(grepl("in 3 or more samples", b$warnings, fixed = TRUE)))
  # A continuous variable: one over the largest leverage, as edgeR's
  # filterByExpr() reads a design.
  ct <- suppressMessages(run_diff(inp, method = "deseq2", analysis_type = "continuous",
                                  continuous_col = "dose", prefilter = TRUE))
  X <- stats::model.matrix(~ dose, inp$meta_df)
  n_min <- ceiling(1 / max(stats::hat(X, intercept = FALSE)) - 1e-8)
  keep_c <- unname(rowSums(inp$expr_mat >= 10) >= n_min)
  expect_identical(is.na(ct$results$diff_result_df$p_value), !keep_c)
  expect_true(any(grepl(sprintf("in %d or more samples", n_min), ct$warnings, fixed = TRUE)))
  # And the global test.
  an <- suppressMessages(run_diff(inp, method = "deseq2", analysis_type = "anova",
                                  group_col = "group", prefilter = TRUE))
  expect_identical(nrow(an$results$diff_result_df), 200L)
  expect_identical(is.na(an$results$diff_result_df$p_value),
                   unname(rowSums(inp$expr_mat >= 10) < 4))
  expect_true(any(grepl("set aside before fitting", an$warnings, fixed = TRUE)))
})

test_that("edgeR filters by default, and prefilter = FALSE fits every gene", {
  skip_if_not_installed("edgeR")
  inp <- pf_counts()
  unset <- pf_group(inp, "edger")
  on <- pf_group(inp, "edger", prefilter = TRUE)
  expect_identical(on$results$diff_result_df, unset$results$diff_result_df)
  expect_true(any(grepl("filterByExpr", unset$warnings)))
  expect_true(anyNA(unset$results$diff_result_df$p_value))
  off <- pf_group(inp, "edger", prefilter = FALSE)
  expect_false(any(grepl("filterByExpr", off$warnings)))
  expect_false(anyNA(off$results$diff_result_df$p_value))
  expect_identical(nrow(off$results$diff_result_df), 200L)
  an <- suppressMessages(run_diff(inp, method = "edger", analysis_type = "anova",
                                  group_col = "group", prefilter = FALSE))
  expect_false(anyNA(an$results$diff_result_df$p_value))
})

test_that("prefilter is a count-model option: others ignore it and say so", {
  skip_if_not_installed("limma")
  set.seed(2)
  m <- matrix(stats::rnorm(40 * 6, 20), 40, 6,
              dimnames = list(paste0("P", 1:40), paste0("S", 1:6)))
  inp <- omics_input(m, data.frame(group = rep(c("A", "B"), each = 3), row.names = colnames(m)),
                     data.frame(feature_id = rownames(m)), omics_type = "proteomics",
                     assay_type = "normalized_intensity")
  plain <- pf_group(inp, "limma")
  expect_warning(asked <- pf_group(inp, "limma", prefilter = TRUE),
                 "method = 'limma' does not support `prefilter`; it was ignored.", fixed = TRUE)
  expect_identical(asked$results$diff_result_df, plain$results$diff_result_df)
  expect_false("prefilter" %in% names(asked$params))
  expect_error(pf_group(inp, "limma", prefilter = "yes"), "prefilter")
  expect_error(pf_group(inp, "limma", prefilter = NA), "prefilter")
})

test_that("an exported script repeats the prefilter choice", {
  expect_true("prefilter" %in% script_arg_names(run_diff))
})
