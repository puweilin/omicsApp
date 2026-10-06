# Gene-length offsets from tximport in the global test across groups.
#
# Counts summarised from transcripts carry each sample's average
# transcript length: a gene whose dominant isoform is four times longer
# in one group has four times the reads there at the same expression.
# The pairwise edgeR and DESeq2 fits corrected for it; the edgeR global
# test (QL F on every group coefficient) normalised by library size
# alone, and called such genes changed.

txg_input <- function(seed = 21) {
  set.seed(seed)
  g <- 300
  groups <- rep(c("A", "B", "C"), each = 4)
  ids <- paste0("G", seq_len(g))
  samp <- paste0("s", seq_along(groups))
  # Transcript length per gene and sample: genes 1-30 switch to an
  # isoform four times as long in group B, with expression unchanged.
  base_len <- stats::runif(g, 800, 3000)
  len <- matrix(base_len, g, length(groups), dimnames = list(ids, samp))
  len <- len * matrix(stats::runif(length(len), 0.95, 1.05), g)
  len[1:30, groups == "B"] <- len[1:30, groups == "B"] * 4
  # Expression: genes 31-45 really go up four-fold in group C.
  expr <- matrix(exp(stats::rnorm(g, log(0.3), 0.8)), g, length(groups))
  expr[31:45, groups == "C"] <- expr[31:45, groups == "C"] * 4
  m <- matrix(stats::rnbinom(length(expr), mu = expr * len, size = 30), g,
              dimnames = list(ids, samp))
  inp <- omics_input(m, data.frame(group = groups, row.names = samp),
                     data.frame(feature_id = ids), omics_type = "rnaseq",
                     assay_type = "raw_count")
  inp$misc <- list(tximport = list(length = len, counts_from_abundance = "no"))
  inp
}

txg_hits <- function(b, rows) {
  df <- b$results$diff_result_df
  df <- df[match(paste0("G", rows), df$feature_id), ]
  sum(!is.na(df$adj_p_value) & df$adj_p_value < 0.05)
}

test_that("the edgeR global test uses tximport's length offsets", {
  skip_if_not_installed("edgeR")
  with_len <- txg_input()
  without <- with_len
  without$misc <- list()
  run <- function(inp) run_diff(inp, method = "edger", analysis_type = "anova",
                                group_col = "group")
  b_off <- run(with_len)
  b_none <- run(without)
  # Without the lengths, the isoform switch reads as a change in expression.
  expect_gte(txg_hits(b_none, 1:30), 25L)
  # With them it does not ...
  expect_lte(txg_hits(b_off, 1:30), 2L)
  # ... and the genes whose expression did change are still found.
  expect_gte(txg_hits(b_off, 31:45), 13L)
})

test_that("the DESeq2 global test uses tximport's length offsets", {
  skip_if_not_installed("DESeq2")
  with_len <- txg_input()
  without <- with_len
  without$misc <- list()
  run <- function(inp) suppressMessages(run_diff(inp, method = "deseq2",
                                                 analysis_type = "anova",
                                                 group_col = "group"))
  b_off <- run(with_len)
  b_none <- run(without)
  expect_gte(txg_hits(b_none, 1:30), 25L)
  expect_lte(txg_hits(b_off, 1:30), 2L)
  expect_gte(txg_hits(b_off, 31:45), 13L)
})

test_that("the pairwise fits and the global test agree about lengths", {
  skip_if_not_installed("edgeR")
  skip_if_not_installed("DESeq2")
  inp <- txg_input()
  for (m in c("edger", "deseq2")) {
    b <- suppressMessages(run_diff(inp, method = m, group_col = "group",
                                   control_group = "A", case_group = "B"))
    expect_lte(txg_hits(b, 1:30), 2L)
  }
})

test_that("counts from abundance get no second length correction", {
  skip_if_not_installed("edgeR")
  inp <- txg_input()
  plain <- inp
  plain$misc <- list()
  inp$misc$tximport$counts_from_abundance <- "lengthScaledTPM"
  a <- run_diff(inp, method = "edger", analysis_type = "anova", group_col = "group")
  z <- run_diff(plain, method = "edger", analysis_type = "anova", group_col = "group")
  expect_identical(a$results$diff_result_df, z$results$diff_result_df)
})

test_that("without length metadata the edgeR global test is edgeR's own", {
  skip_if_not_installed("edgeR")
  inp <- txg_input()
  inp$misc <- list()
  b <- run_diff(inp, method = "edger", analysis_type = "anova", group_col = "group")
  # The same fit written out by hand: filterByExpr, TMM, QL F-test on
  # the two group coefficients.
  grp <- factor(inp$meta_df$group)
  design <- stats::model.matrix(~ grp)
  y <- edgeR::DGEList(counts = inp$expr_mat)
  keep <- edgeR::filterByExpr(y, design = design)
  y <- edgeR::calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
  y <- edgeR::estimateDisp(y, design = design)
  fit <- edgeR::glmQLFit(y, design = design)
  tt <- edgeR::topTags(edgeR::glmQLFTest(fit, coef = 2:3), n = Inf, sort.by = "none")$table
  df <- b$results$diff_result_df
  df <- df[match(rownames(tt), df$feature_id), ]
  expect_equal(df$statistic, tt$F)
  expect_equal(df$p_value, tt$PValue)
  expect_equal(df$adj_p_value, tt$FDR)
})
