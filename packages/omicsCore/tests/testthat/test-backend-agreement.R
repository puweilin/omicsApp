# Each differential engine, reproduced by hand from the same input.
#
# run_diff() wraps limma, DESeq2, edgeR, a Welch t-test and lm, and
# maps their output onto one schema. A wrapper can drift from the
# library it wraps -- a changed default, a contrast built the other way
# round, a column mapped to the wrong slot -- and produce a complete,
# plausible, wrong table. The check is the one this project already
# uses for normalization: two implementations of one transform have to
# agree exactly, not approximately.

by_feature <- function(df, ids) df[match(ids, df$feature_id), , drop = FALSE]

expect_direction_follows_effect <- function(df) {
  expect_identical(df$direction[df$effect > 0 & !is.na(df$effect)],
                   rep("up", sum(df$effect > 0, na.rm = TRUE)))
  expect_identical(df$direction[df$effect < 0 & !is.na(df$effect)],
                   rep("down", sum(df$effect < 0, na.rm = TRUE)))
}

expect_one_row_per_feature <- function(df, input) {
  expect_setequal(df$feature_id, rownames(input$expr_mat))
  expect_false(anyDuplicated(df$feature_id) > 0L)
  expect_true(all(is.na(df$is_significant)))
}

# ---- limma ----------------------------------------------------------------

test_that("limma group with a covariate matches lmFit / contrasts.fit / eBayes", {
  skip_if_not_installed("limma")
  inp <- realistic_input()
  b <- run_diff(inp, method = "limma", analysis_type = "group",
                group_col = "group", control_group = "G1", case_group = "G2",
                covariates = "age")
  got <- b$results$diff_result_df

  meta <- inp$meta_df
  meta$group <- factor(meta$group, levels = c("G1", "G2"))
  design <- stats::model.matrix(~ 0 + group + age, data = meta)
  colnames(design) <- sub("^group", "", colnames(design))
  fit <- limma::lmFit(inp$expr_mat[, rownames(meta)], design)
  contrast <- limma::makeContrasts(contrasts = "G2 - G1", levels = design)
  fit2 <- limma::eBayes(limma::contrasts.fit(fit, contrast))
  ref <- limma::topTable(fit2, number = Inf, sort.by = "none")

  got <- by_feature(got, rownames(ref))
  expect_equal(got$effect, ref$logFC, tolerance = 1e-12)
  expect_equal(got$statistic, ref$t, tolerance = 1e-12)
  expect_equal(got$p_value, ref$P.Value, tolerance = 1e-12)
  expect_equal(got$adj_p_value, ref$adj.P.Val, tolerance = 1e-12)
  expect_equal(got$base_mean, ref$AveExpr, tolerance = 1e-12)
  expect_identical(unique(got$effect_type), "log2FC")
  expect_identical(unique(got$statistic_type), "t")
  expect_identical(unique(got$method), "limma")
  expect_direction_follows_effect(got)
  expect_one_row_per_feature(got, inp)
})

test_that("limma continuous matches topTable(coef = 2) and a Spearman correlation", {
  skip_if_not_installed("limma")
  inp <- realistic_input()
  b <- run_diff(inp, method = "limma", analysis_type = "continuous",
                continuous_col = "age")
  got <- b$results$diff_result_df

  design <- stats::model.matrix(~ age, data = inp$meta_df)
  fit <- limma::eBayes(limma::lmFit(inp$expr_mat, design))
  ref <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "none")
  rho <- apply(inp$expr_mat, 1L, function(y) {
    suppressWarnings(stats::cor.test(y, inp$meta_df$age, method = "spearman",
                                     exact = FALSE)$estimate)
  })

  got <- by_feature(got, rownames(ref))
  expect_equal(got$effect, unname(rho[rownames(ref)]), tolerance = 1e-12)
  expect_equal(got$statistic, ref$t, tolerance = 1e-12)
  expect_equal(got$p_value, ref$P.Value, tolerance = 1e-12)
  expect_equal(got$adj_p_value, ref$adj.P.Val, tolerance = 1e-12)
  expect_identical(unique(got$effect_type), "correlation")
  expect_one_row_per_feature(got, inp)
})

# ---- t-test and lm --------------------------------------------------------

test_that("ttest matches per-feature Welch t.test and BH", {
  inp <- realistic_input()
  b <- run_diff(inp, method = "ttest", analysis_type = "group",
                group_col = "group", control_group = "G1", case_group = "G2")
  got <- b$results$diff_result_df

  ctrl <- rownames(inp$meta_df)[inp$meta_df$group == "G1"]
  case <- rownames(inp$meta_df)[inp$meta_df$group == "G2"]
  ref <- t(apply(inp$expr_mat, 1L, function(y) {
    tt <- stats::t.test(y[case], y[ctrl], var.equal = FALSE)
    c(effect = mean(y[case]) - mean(y[ctrl]), t = unname(tt$statistic),
      p = tt$p.value, base = (mean(y[case]) + mean(y[ctrl])) / 2)
  }))

  got <- by_feature(got, rownames(ref))
  expect_equal(got$effect, unname(ref[, "effect"]), tolerance = 1e-12)
  expect_equal(got$statistic, unname(ref[, "t"]), tolerance = 1e-12)
  expect_equal(got$p_value, unname(ref[, "p"]), tolerance = 1e-12)
  expect_equal(got$adj_p_value, unname(stats::p.adjust(ref[, "p"], "BH")), tolerance = 1e-12)
  expect_equal(got$base_mean, unname(ref[, "base"]), tolerance = 1e-12)
  expect_identical(unique(got$effect_type), "mean_diff")
  expect_direction_follows_effect(got)
  expect_one_row_per_feature(got, inp)
})

test_that("lm group with a covariate matches per-feature lm()", {
  inp <- realistic_input()
  b <- run_diff(inp, method = "lm", analysis_type = "group",
                group_col = "group", control_group = "G1", case_group = "G2",
                covariates = "age")
  got <- b$results$diff_result_df

  meta <- inp$meta_df
  meta$group <- factor(meta$group, levels = c("G1", "G2"))
  ref <- t(apply(inp$expr_mat[, rownames(meta)], 1L, function(y) {
    fit <- stats::lm(y ~ group + age, data = data.frame(y = y, meta))
    s <- summary(fit)
    c(beta = s$coefficients["groupG2", "Estimate"],
      t = s$coefficients["groupG2", "t value"],
      p = s$coefficients["groupG2", "Pr(>|t|)"],
      r2 = s$adj.r.squared, base = mean(y))
  }))

  got <- by_feature(got, rownames(ref))
  expect_equal(got$effect, unname(ref[, "beta"]), tolerance = 1e-12)
  expect_equal(got$statistic, unname(ref[, "t"]), tolerance = 1e-12)
  expect_equal(got$p_value, unname(ref[, "p"]), tolerance = 1e-12)
  expect_equal(got$adj_p_value, unname(stats::p.adjust(ref[, "p"], "BH")), tolerance = 1e-12)
  expect_equal(got$model_fit, unname(ref[, "r2"]), tolerance = 1e-12)
  expect_equal(got$base_mean, unname(ref[, "base"]), tolerance = 1e-12)
  expect_identical(unique(got$effect_type), "beta")
  expect_direction_follows_effect(got)
  expect_one_row_per_feature(got, inp)
})

# ---- the count engines ----------------------------------------------------

test_that("deseq2 group matches DESeqDataSetFromMatrix / DESeq / results", {
  skip_if_not_installed("DESeq2")
  inp <- realistic_input("rnaseq", n_per_group = 5L)
  b <- suppressMessages(run_diff(
    inp, method = "deseq2", analysis_type = "group",
    group_col = "group", control_group = "G1", case_group = "G2"
  ))
  got <- b$results$diff_result_df

  meta <- inp$meta_df
  meta$group <- stats::relevel(factor(meta$group), ref = "G1")
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(inp$expr_mat[, rownames(meta)]),
    colData = meta, design = ~ group
  )
  dds <- suppressMessages(DESeq2::DESeq(dds, quiet = TRUE))
  ref <- as.data.frame(DESeq2::results(dds, contrast = c("group", "G2", "G1")))

  got <- by_feature(got, rownames(ref))
  # 1e-8 rather than 1e-12: the two runs are the same code path, but
  # DESeq2's fit is iterative and sums in an order that can differ by
  # one ulp from the wrapper's, which subsets the matrix before fitting.
  expect_equal(got$effect, ref$log2FoldChange, tolerance = 1e-8)
  expect_equal(got$statistic, ref$stat, tolerance = 1e-8)
  expect_equal(got$p_value, ref$pvalue, tolerance = 1e-8)
  expect_equal(got$adj_p_value, ref$padj, tolerance = 1e-8)
  expect_equal(got$base_mean, ref$baseMean, tolerance = 1e-8)
  expect_identical(unique(got$effect_type), "log2FC")
  expect_identical(unique(got$statistic_type), "wald")
  expect_direction_follows_effect(got)
  expect_one_row_per_feature(got, inp)
  # The perturbed block is the one that moves. Relative, not absolute:
  # a quarter of the genes shifted one way is enough to pull the
  # library-size normalisation, so the untouched genes sit at about
  # -0.35 rather than 0 -- a property of median-of-ratios, not a fault.
  hit <- got$feature_symbol %in% REAL_GENE_SETS$G2M
  expect_gt(median(got$effect[hit]) - median(got$effect[!hit], na.rm = TRUE), 1)
})

test_that("edger group matches DGEList / estimateDisp / glmQLFit / glmQLFTest", {
  skip_if_not_installed("edgeR")
  inp <- realistic_input("rnaseq", n_per_group = 5L)
  b <- run_diff(inp, method = "edger", analysis_type = "group",
                group_col = "group", control_group = "G1", case_group = "G2")
  got <- b$results$diff_result_df

  meta <- inp$meta_df
  meta$group <- factor(meta$group, levels = c("G1", "G2"))
  design <- stats::model.matrix(~ group, data = meta)
  y <- edgeR::DGEList(counts = inp$expr_mat[, rownames(meta)])
  y <- edgeR::calcNormFactors(y)
  y <- edgeR::estimateDisp(y, design = design)
  fit <- edgeR::glmQLFit(y, design = design)
  qlf <- edgeR::glmQLFTest(fit, coef = "groupG2")
  ref <- as.data.frame(edgeR::topTags(qlf, n = Inf, sort.by = "none")$table)

  got <- by_feature(got, rownames(ref))
  expect_equal(got$effect, ref$logFC, tolerance = 1e-8)
  expect_equal(got$statistic, ref$F, tolerance = 1e-8)
  expect_equal(got$p_value, ref$PValue, tolerance = 1e-8)
  expect_equal(got$adj_p_value, ref$FDR, tolerance = 1e-8)
  expect_equal(got$base_mean, ref$logCPM, tolerance = 1e-8)
  expect_identical(unique(got$effect_type), "log2FC")
  expect_identical(unique(got$statistic_type), "F")
  expect_direction_follows_effect(got)
  expect_one_row_per_feature(got, inp)
})

test_that("the count engines agree with each other on what moved", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("edgeR")
  inp <- realistic_input("rnaseq", n_per_group = 5L)
  run <- function(m) suppressMessages(run_diff(
    inp, method = m, analysis_type = "group",
    group_col = "group", control_group = "G1", case_group = "G2"
  ))$results$diff_result_df
  a <- run("deseq2")
  z <- by_feature(run("edger"), a$feature_id)
  # Same sign on every feature with a clear effect, and highly
  # correlated fold changes: two estimators of the same quantity.
  clear <- abs(a$effect) > 0.5 & !is.na(a$effect)
  expect_true(all(sign(a$effect[clear]) == sign(z$effect[clear])))
  expect_gt(stats::cor(a$effect, z$effect, use = "complete.obs"), 0.95)
  top_a <- a$feature_id[order(a$adj_p_value)][1:20]
  top_z <- z$feature_id[order(z$adj_p_value)][1:20]
  expect_gt(length(intersect(top_a, top_z)), 12L)
})
