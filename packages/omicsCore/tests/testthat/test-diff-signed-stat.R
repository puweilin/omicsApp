# `signed_stat`: one signed test statistic for every engine.
#
# `statistic` is whatever the engine reported -- limma's t, DESeq2's Wald
# statistic, edgeR's unsigned QL F -- so a table from edgeR could not be
# ranked or compared with one from limma without knowing which. The
# optional `signed_stat` column puts all of them on the t scale with the
# sign of the effect, and is NA where a test has no direction.

ss_counts_input <- function(groups = rep(c("A", "B", "C"), each = 4), seed = 11) {
  set.seed(seed)
  g <- 200
  mu <- matrix(exp(stats::rnorm(g, log(300), 1)), g, length(groups))
  mu[1:20, groups == "B"] <- mu[1:20, groups == "B"] * 4
  mu[21:40, groups == "B"] <- mu[21:40, groups == "B"] / 4
  m <- matrix(stats::rnbinom(length(mu), mu = mu, size = 20), g,
              dimnames = list(paste0("G", seq_len(g)), paste0("s", seq_along(groups))))
  omics_input(m, data.frame(group = groups, dose = seq_along(groups) + stats::rnorm(length(groups)),
                            row.names = colnames(m)),
              data.frame(feature_id = rownames(m)), omics_type = "rnaseq",
              assay_type = "raw_count")
}

ss_prot_input <- function(seed = 12) {
  set.seed(seed)
  groups <- rep(c("A", "B", "C"), each = 5)
  m <- matrix(stats::rnorm(60 * 15, 10, 0.5), 60,
              dimnames = list(paste0("P", 1:60), paste0("s", 1:15)))
  m[1:10, groups == "B"] <- m[1:10, groups == "B"] + 2
  m[11:20, groups == "B"] <- m[11:20, groups == "B"] - 2
  dose <- seq_len(15) + stats::rnorm(15)
  m[21:30, ] <- m[21:30, ] + rep(0.3 * dose, each = 10)
  omics_input(m, data.frame(group = groups, dose = dose, row.names = colnames(m)),
              data.frame(feature_id = rownames(m)), omics_type = "proteomics",
              assay_type = "normalized_intensity")
}

# The sign of signed_stat agrees with the sign of the log2FC, wherever
# both are defined and non-zero.
expect_signed_like_effect <- function(df) {
  ok <- is.finite(df$signed_stat) & is.finite(df$effect) & df$effect != 0 &
    df$signed_stat != 0
  expect_gt(sum(ok), 10L)
  expect_identical(sign(df$signed_stat[ok]), sign(df$effect[ok]))
}

test_that("limma, the t-test and lm report their t as signed_stat", {
  skip_if_not_installed("limma")
  inp <- ss_prot_input()
  for (m in c("limma", "ttest", "lm")) {
    df <- run_diff(inp, method = m, group_col = "group", control_group = "A",
                   case_group = "B")$results$diff_result_df
    expect_true("signed_stat" %in% names(df), info = m)
    expect_identical(df$signed_stat, df$statistic, info = m)
    expect_signed_like_effect(df)
  }
})

test_that("continuous fits report the slope's t; a spline has no sign", {
  skip_if_not_installed("limma")
  inp <- ss_prot_input()
  for (m in c("limma", "lm")) {
    df <- run_diff(inp, method = m, analysis_type = "continuous",
                   continuous_col = "dose")$results$diff_result_df
    # The effect is Spearman's rho, which can disagree with the slope;
    # the signed statistic is the slope's t.
    expect_identical(df$signed_stat, df$statistic, info = m)
    expect_true(all(df$signed_stat[21:30] > 0), info = m)
  }
  sp <- run_diff(inp, method = "limma", analysis_type = "continuous",
                 continuous_col = "dose", model = "spline", df = 3)$results$diff_result_df
  expect_identical(unique(sp$statistic_type), "F")
  expect_true(all(is.na(sp$signed_stat)))
})

test_that("DESeq2 reports its Wald statistic as signed_stat", {
  skip_if_not_installed("DESeq2")
  inp <- ss_counts_input()
  df <- suppressMessages(run_diff(inp, method = "deseq2", group_col = "group",
                                  control_group = "A", case_group = "B"))$results$diff_result_df
  expect_identical(df$signed_stat, df$statistic)
  expect_signed_like_effect(df)
  ct <- suppressMessages(run_diff(inp, method = "deseq2", analysis_type = "continuous",
                                  continuous_col = "dose"))$results$diff_result_df
  expect_identical(ct$signed_stat, ct$statistic)
})

test_that("edgeR's signed_stat is sign(logFC) * sqrt(F), on the t scale", {
  skip_if_not_installed("edgeR")
  inp <- ss_counts_input()
  df <- run_diff(inp, method = "edger", group_col = "group", control_group = "A",
                 case_group = c("B", "C"))$results$diff_result_df
  expect_identical(unique(df$statistic_type), "F")
  expect_equal(df$signed_stat, sign(df$effect) * sqrt(df$statistic))
  expect_signed_like_effect(df)
  # Genes the filter set aside have no statistic, signed or not.
  expect_identical(is.na(df$signed_stat), is.na(df$statistic))
  # On the same scale as a t: the planted 4-fold changes stand well out.
  b <- df[df$comparison == "B_vs_A", ]
  expect_true(all(b$signed_stat[1:20] > 3, na.rm = TRUE))
  expect_true(all(b$signed_stat[21:40] < -3, na.rm = TRUE))
})

test_that("a global test has no signed statistic", {
  skip_if_not_installed("limma")
  skip_if_not_installed("edgeR")
  skip_if_not_installed("DESeq2")
  an <- run_diff(ss_prot_input(), method = "limma", analysis_type = "anova",
                 group_col = "group")$results$diff_result_df
  expect_true("signed_stat" %in% names(an))
  expect_true(all(is.na(an$signed_stat)))
  cnt <- ss_counts_input()
  for (m in c("edger", "deseq2")) {
    df <- suppressMessages(run_diff(cnt, method = m, analysis_type = "anova",
                                    group_col = "group"))$results$diff_result_df
    expect_true("signed_stat" %in% names(df), info = m)
    expect_true(all(is.na(df$signed_stat)), info = m)
  }
})

test_that("GSEA ranks by signed_stat, and the ranking is the one it was", {
  skip_if_not_installed("limma")
  skip_if_not_installed("edgeR")
  skip_if_not_installed("DESeq2")
  prot <- ss_prot_input()
  cnt <- ss_counts_input()
  runs <- list(
    limma = run_diff(prot, method = "limma", group_col = "group",
                     control_group = "A", case_group = "B"),
    ttest = run_diff(prot, method = "ttest", group_col = "group",
                     control_group = "A", case_group = "B"),
    deseq2 = suppressMessages(run_diff(cnt, method = "deseq2", group_col = "group",
                                       control_group = "A", case_group = "B")),
    edger = run_diff(cnt, method = "edger", group_col = "group",
                     control_group = "A", case_group = "B")
  )
  for (nm in names(runs)) {
    df <- runs[[nm]]$results$diff_result_df
    old <- df
    old$signed_stat <- NULL
    new_rank <- gsea_rank_vector(df, "feature_id")
    old_rank <- gsea_rank_vector(old, "feature_id")
    # The same genes, in the same order, with the same values: a bundle
    # saved before the column and one made after it give one GSEA.
    expect_identical(names(new_rank), names(old_rank), info = nm)
    expect_identical(as.numeric(new_rank), as.numeric(old_rank), info = nm)
    expect_identical(attr(new_rank, "metric"), attr(old_rank, "metric"), info = nm)
    # And it is signed_stat that is ranked.
    expect_equal(as.numeric(new_rank), unname(df$signed_stat[match(names(new_rank), df$feature_id)]),
                 info = nm)
  }
  # The column wins over the statistic when both are there.
  df <- runs$edger$results$diff_result_df
  df$signed_stat <- -df$signed_stat
  r <- gsea_rank_vector(df, "feature_id")
  expect_equal(as.numeric(r), unname(df$signed_stat[match(names(r), df$feature_id)]))
})

test_that("limma GSEA results are unchanged by the new column", {
  skip_if_not_installed("limma")
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("msigdbr")
  skip_if_not_installed("fgsea")
  b <- realistic_diff_bundle()
  old <- b
  old$results$diff_result_df$signed_stat <- NULL
  a <- suppressMessages(suppressWarnings(run_enrichment(b, type = "gsea", database = "hallmark")))
  z <- suppressMessages(suppressWarnings(run_enrichment(old, type = "gsea", database = "hallmark")))
  ea <- a$results$enrich_result_df
  ez <- z$results$enrich_result_df
  expect_gt(nrow(ea), 0L)
  expect_identical(ea$pathway_id, ez$pathway_id)
  expect_equal(ea$effect, ez$effect)
  expect_equal(ea$p_value, ez$p_value)
  expect_identical(a$params$rank_metric, z$params$rank_metric)
})

test_that("a bundle saved before signed_stat existed still works everywhere", {
  skip_if_not_installed("limma")
  inp <- ss_prot_input()
  b <- run_diff(inp, method = "limma", group_col = "group", control_group = "A",
                case_group = c("B", "C"))
  old <- b
  old$results$diff_result_df$signed_stat <- NULL
  expect_silent(check_diff_result_schema(old$results$diff_result_df))
  expect_identical(diff_comparisons(old), c("B_vs_A", "C_vs_A"))
  one <- select_comparison(old, "B_vs_A")
  expect_false("signed_stat" %in% names(one$results$diff_result_df))
  expect_gt(nrow(filter_diff_results(one$results$diff_result_df)), 0L)
  expect_identical(summarize_diff_contrasts(old), summarize_diff_contrasts(b))
  expect_s3_class(plot_volcano(one), "ggplot")
  expect_s3_class(plot_ma(one), "ggplot")
  expect_s3_class(plot_diff_contrasts(old), "ggplot")
  expect_identical(gsea_rank_metric(one), "signed test statistic")

  # An old bundle and a new one, side by side in an integration.
  p <- omics_project("q", list(a = inp, b = inp))
  int <- run_integration(p, "concordance", c("a", "b"),
                         diff_bundles = list(a = one, b = select_comparison(b, "B_vs_A")))
  expect_gt(nrow(int$results$integration_df), 0L)

  # And through a save and a load.
  p$bundles <- list(old = one)
  path <- withr::local_tempfile(fileext = ".rds")
  save_project(p, path)
  back <- load_project(path)
  expect_identical(back$bundles$old$results$diff_result_df, one$results$diff_result_df)
})

test_that("the empty template carries the optional column", {
  tpl <- new_diff_result_template()
  expect_true(all(DIFF_RESULT_OPTIONAL_COLS %in% names(tpl)))
  expect_type(tpl$signed_stat, "double")
  expect_false(any(DIFF_RESULT_OPTIONAL_COLS %in% DIFF_RESULT_REQUIRED_COLS))
})
