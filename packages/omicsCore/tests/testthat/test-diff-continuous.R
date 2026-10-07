# run_diff(analysis_type = "continuous") on raw counts runs DESeq2 with
# the variable as a slope: `~ [pair +] variable [+ covariates]`, Wald test
# on the variable's coefficient, effect in log2 fold change per unit. It
# is what `method = "auto"` picks for a counts layer, and it had no test
# of its own.

dc_counts <- function(n_genes = 400L, seed = 11) {
  set.seed(seed)
  n <- 16L
  meta <- data.frame(
    dose = rep(c(0, 1, 2.5, 5, 7.5, 10, 12.5, 15), 2),
    batch = rep(c("b1", "b2"), each = 8),
    # Each subject seen at two neighbouring doses, so a pairing block
    # leaves the within-subject change in dose to estimate the slope from.
    subject = rep(paste0("p", 1:8), each = 2),
    row.names = sprintf("S%02d", seq_len(n))
  )
  mu <- exp(stats::runif(n_genes, log(50), log(2000)))
  slope <- numeric(n_genes)
  slope[1:20] <- 0.15     # log2 fold change per unit of dose
  slope[21:40] <- -0.15
  batch_eff <- ifelse(meta$batch == "b2", 1.6, 1)
  m <- outer(mu, batch_eff) * 2^outer(slope, meta$dose)
  counts <- matrix(stats::rnbinom(n_genes * n, mu = m, size = 20), n_genes, n,
                   dimnames = list(sprintf("G%03d", seq_len(n_genes)), rownames(meta)))
  omics_input(counts, meta, data.frame(feature_id = rownames(counts)),
              omics_type = "rnaseq", assay_type = "raw_count")
}

dc_run <- function(inp, ...) {
  suppressMessages(run_diff(inp, method = "deseq2", analysis_type = "continuous",
                            continuous_col = "dose", ...))
}

# DESeq2 as a user would call it, on the same counts and design.
dc_direct <- function(inp, design, name = "dose", txi = NULL) {
  cd <- inp$meta_df
  for (v in intersect(c("batch", "subject"), names(cd))) cd[[v]] <- factor(cd[[v]])
  dds <- if (is.null(txi)) {
    DESeq2::DESeqDataSetFromMatrix(inp$expr_mat, cd, design)
  } else {
    DESeq2::DESeqDataSetFromTximport(txi, cd, design)
  }
  dds <- suppressMessages(DESeq2::DESeq(dds, quiet = TRUE))
  as.data.frame(DESeq2::results(dds, name = name))
}

test_that("DESeq2 continuous recovers a planted dose trend, with the right sign", {
  skip_if_not_installed("DESeq2")
  df <- dc_run(dc_counts())$results$diff_result_df
  up <- sprintf("G%03d", 1:20)
  down <- sprintf("G%03d", 21:40)
  sig <- df$feature_id[!is.na(df$adj_p_value) & df$adj_p_value < 0.05]
  expect_gte(length(intersect(sig, c(up, down))), 36L)
  expect_lte(length(setdiff(sig, c(up, down))), 6L)
  rownames(df) <- df$feature_id
  # The slope is on the planted scale: log2 fold change per unit of dose.
  expect_equal(stats::median(df[up, "effect"]), 0.15, tolerance = 0.2)
  expect_equal(stats::median(df[down, "effect"]), -0.15, tolerance = 0.2)
  expect_true(all(df[up, "direction"] == "positive"))
  expect_true(all(df[down, "direction"] == "negative"))
  expect_true(all(df[up, "signed_stat"] > 0))
  expect_true(all(df[down, "signed_stat"] < 0))
})

test_that("DESeq2 continuous returns the standard schema, one row per gene", {
  skip_if_not_installed("DESeq2")
  inp <- dc_counts()
  b <- dc_run(inp)
  df <- b$results$diff_result_df
  expect_silent(check_diff_result_schema(df))
  expect_identical(names(df), c(DIFF_RESULT_REQUIRED_COLS[1:11], "signed_stat",
                                DIFF_RESULT_REQUIRED_COLS[12:17]))
  expect_identical(df$feature_id, rownames(inp$expr_mat))
  expect_identical(unique(df$method), "deseq2")
  expect_identical(unique(df$analysis_type), "continuous_linear")
  expect_identical(unique(df$comparison), "dose")
  expect_identical(unique(df$effect_type), "log2FC_per_unit")
  expect_identical(unique(df$statistic_type), "wald")
  expect_identical(df$signed_stat, df$statistic)
  expect_true(all(is.na(df$is_significant)))
  ok <- !is.na(df$effect) & df$effect != 0 & !is.na(df$signed_stat)
  expect_identical(sign(df$signed_stat[ok]), sign(df$effect[ok]))
  raw <- b$results$diff_raw_df
  expect_identical(df$base_mean, raw$baseMean)
  expect_identical(df$p_value, raw$pvalue)
  expect_identical(b$params$analysis_type, "continuous")
  expect_identical(b$params$comparison, "dose")
})

test_that("DESeq2 continuous agrees with calling DESeq2 directly", {
  skip_if_not_installed("DESeq2")
  inp <- dc_counts()
  df <- dc_run(inp)$results$diff_result_df
  ref <- dc_direct(inp, ~ dose)
  expect_equal(df$effect, ref$log2FoldChange, tolerance = 1e-8)
  expect_equal(df$statistic, ref$stat, tolerance = 1e-8)
  expect_equal(df$p_value, ref$pvalue, tolerance = 1e-8)
  expect_equal(df$adj_p_value, ref$padj, tolerance = 1e-8)
})

test_that("with covariates it is the Wald test of the slope in ~ covariates + dose", {
  skip_if_not_installed("DESeq2")
  inp <- dc_counts()
  df <- dc_run(inp, covariates = "batch")$results$diff_result_df
  # The design written the way the DESeq2 vignette writes it, variable of
  # interest last; the order of terms does not change the slope's test.
  ref <- dc_direct(inp, ~ batch + dose)
  expect_equal(df$effect, ref$log2FoldChange, tolerance = 1e-6)
  expect_equal(df$p_value, ref$pvalue, tolerance = 1e-6)
  expect_equal(df$adj_p_value, ref$padj, tolerance = 1e-6)
  # The adjustment matters here: batch moves every gene 1.6-fold and is
  # correlated with nothing, so leaving it out costs power, not bias.
  plain <- dc_run(inp)$results$diff_result_df
  expect_false(isTRUE(all.equal(df$p_value, plain$p_value)))
})

test_that("a pairing column enters DESeq2 as a block", {
  skip_if_not_installed("DESeq2")
  inp <- dc_counts()
  df <- dc_run(inp, paired_col = "subject")$results$diff_result_df
  ref <- dc_direct(inp, ~ subject + dose)
  expect_equal(df$p_value, ref$pvalue, tolerance = 1e-6)
  expect_equal(df$effect, ref$log2FoldChange, tolerance = 1e-6)
})

test_that("a variable that does not change within pairs is refused in words", {
  skip_if_not_installed("DESeq2")
  inp <- dc_counts(n_genes = 60L)
  # One dose per subject: with the subject as a block there is no
  # within-subject change left. DESeq2 stopped on "the model matrix is
  # not full rank".
  inp$meta_df$dose <- rep(c(0, 2.5, 5, 7.5, 10, 12.5, 15, 20), each = 2)
  expect_error(dc_run(inp, paired_col = "subject"),
               "`dose` does not vary within the blocks of `subject`", fixed = TRUE)
  # limma models the pairs as a correlation when it has to, and runs.
  skip_if_not_installed("limma")
  expect_no_error(suppressWarnings(run_diff(inp, method = "limma", analysis_type = "continuous",
                                            continuous_col = "dose", paired_col = "subject")))
  # The same for subjects nested in the groups of a global test.
  inp$meta_df$group <- rep(c("A", "B"), each = 8)
  expect_error(suppressMessages(run_diff(inp, method = "deseq2", analysis_type = "anova",
                                         group_col = "group", paired_col = "subject")),
               "`group` does not vary within the blocks of `subject`", fixed = TRUE)
})

test_that("reversing or rescaling the variable moves the slope, not the test", {
  skip_if_not_installed("DESeq2")
  inp <- dc_counts()
  base <- dc_run(inp)$results$diff_result_df
  neg <- inp
  neg$meta_df$dose <- -inp$meta_df$dose
  flipped <- dc_run(neg)$results$diff_result_df
  expect_equal(flipped$effect, -base$effect, tolerance = 1e-6)
  expect_equal(flipped$signed_stat, -base$signed_stat, tolerance = 1e-6)
  expect_equal(flipped$p_value, base$p_value, tolerance = 1e-6)
  ug <- inp
  ug$meta_df$dose <- inp$meta_df$dose * 1000   # mg -> ug
  scaled <- dc_run(ug)$results$diff_result_df
  expect_equal(scaled$effect * 1000, base$effect, tolerance = 1e-5)
  expect_equal(scaled$p_value, base$p_value, tolerance = 1e-5)
})

test_that("a variable stored as text or as a factor of numbers is read as numbers", {
  skip_if_not_installed("DESeq2")
  inp <- dc_counts(n_genes = 120L)
  base <- dc_run(inp)$results$diff_result_df
  as_text <- inp
  as_text$meta_df$dose <- format(inp$meta_df$dose)
  expect_equal(dc_run(as_text)$results$diff_result_df$p_value, base$p_value)
  as_factor <- inp
  as_factor$meta_df$dose <- factor(inp$meta_df$dose)
  expect_equal(dc_run(as_factor)$results$diff_result_df$p_value, base$p_value)
})

test_that("a variable that is not numbers stops with an error naming the values", {
  inp <- dc_counts(n_genes = 60L)
  words <- inp
  words$meta_df$dose <- rep(c("low", "mid", "high", "none"), 4)
  expect_error(dc_run(words),
               "`dose` must hold numbers to be used as a continuous variable; 16 sample(s) have values that are not numbers: 'low', 'mid', 'high', ...",
               fixed = TRUE)
  # One typo among numbers is named as such -- it used to be reported as
  # a missing value, in a sample whose value was there.
  typo <- inp
  typo$meta_df$dose <- as.character(inp$meta_df$dose)
  typo$meta_df$dose[3] <- "10 mg"
  expect_error(dc_run(typo), "1 sample(s) have values that are not numbers: '10 mg'.",
               fixed = TRUE)
  # A blank is a missing value, and said to be one.
  blank <- inp
  blank$meta_df$dose <- as.character(inp$meta_df$dose)
  blank$meta_df$dose[5] <- ""
  expect_error(dc_run(blank), "`dose` has missing values in 1 sample(s): S05.", fixed = TRUE)
  expect_error(dc_run(inp, continuous_col = NULL), "continuous_col")
})

test_that("tximport gene lengths enter the continuous fit as offsets", {
  skip_if_not_installed("DESeq2")
  inp <- dc_counts()
  # Genes 41-60 do not change with dose, but their effective length does
  # (a longer isoform takes over), so their raw counts rise with dose.
  len <- matrix(1500, nrow(inp$expr_mat), ncol(inp$expr_mat),
                dimnames = dimnames(inp$expr_mat))
  stretch <- 2^(0.12 * inp$meta_df$dose)
  len[41:60, ] <- outer(rep(1500, 20), stretch)
  inp$expr_mat[41:60, ] <- round(sweep(inp$expr_mat[41:60, ], 2L, stretch, "*"))
  txi_inp <- inp
  txi_inp$misc$tximport <- list(length = len, counts_from_abundance = "no")

  lengthy <- sprintf("G%03d", 41:60)
  without <- dc_run(inp)$results$diff_result_df
  with <- dc_run(txi_inp)$results$diff_result_df
  rownames(without) <- without$feature_id
  rownames(with) <- with$feature_id
  # Read as a dose trend without the lengths, and not with them.
  expect_gte(sum(without[lengthy, "adj_p_value"] < 0.05, na.rm = TRUE), 15L)
  expect_lte(sum(with[lengthy, "adj_p_value"] < 0.05, na.rm = TRUE), 2L)
  # The genuine trends are still found.
  expect_gte(sum(with[sprintf("G%03d", 1:40), "adj_p_value"] < 0.05, na.rm = TRUE), 36L)
  # And it is what DESeqDataSetFromTximport() gives directly.
  ref <- dc_direct(inp, ~ dose, txi = list(counts = inp$expr_mat, length = len,
                                          countsFromAbundance = "no"))
  expect_equal(with$p_value, ref$pvalue, tolerance = 1e-8)
  expect_equal(with$effect, ref$log2FoldChange, tolerance = 1e-8)
})

test_that("raw counts offer DESeq2 for a continuous variable, and auto picks it", {
  skip_if_not_installed("DESeq2")
  inp <- dc_counts(n_genes = 60L)
  m <- applicable_diff_methods(inp, analysis_type = "continuous")
  expect_true("deseq2" %in% m)
  # edgeR's backend here has no continuous mode.
  expect_false("edger" %in% m)
  expect_true(auto_select_diff_method(inp, "continuous") %in% m)
  b <- suppressMessages(run_diff(inp, analysis_type = "continuous", continuous_col = "dose"))
  expect_identical(b$params$method, "deseq2")
})

test_that("limma's continuous fit with a covariate or pairing survives missing values", {
  skip_if_not_installed("limma")
  set.seed(5)
  n <- 12L
  meta <- data.frame(age = seq(20, 75, length.out = n), sex = rep(c("F", "M"), 6),
                     subject = rep(paste0("p", 1:6), each = 2),
                     row.names = sprintf("S%02d", seq_len(n)))
  m <- matrix(stats::rnorm(50 * n, 20), 50, n,
              dimnames = list(sprintf("P%02d", 1:50), rownames(meta)))
  m[1:10, ] <- m[1:10, ] + outer(rep(0.05, 10), meta$age)
  m[cbind(c(3, 15, 15, 40), c(2, 5, 9, 11))] <- NA
  inp <- omics_input(m, meta, data.frame(feature_id = rownames(m)),
                     omics_type = "proteomics", assay_type = "normalized_intensity")
  gappy <- c("P03", "P15", "P40")
  for (args in list(list(covariates = "sex"), list(paired_col = "subject"))) {
    # It stopped with "NA/NaN/Inf in foreign function call (arg 5)".
    df <- do.call(run_diff, c(list(inp, method = "limma", analysis_type = "continuous",
                                   continuous_col = "age"), args))$results$diff_result_df
    expect_identical(nrow(df), 50L)
    # The partial rank correlation needs every sample; the test does not.
    expect_true(all(is.na(df$effect[df$feature_id %in% gappy])))
    expect_false(anyNA(df$effect[!df$feature_id %in% gappy]))
    expect_false(anyNA(df$p_value))
    # Each complete feature's rho is its own, whatever the other rows hold.
    full <- do.call(run_diff, c(list(subset_omics_features(inp, setdiff(rownames(m), gappy)),
                                     method = "limma", analysis_type = "continuous",
                                     continuous_col = "age"), args))$results$diff_result_df
    expect_equal(df$effect[!df$feature_id %in% gappy], full$effect)
    if (!is.null(args$covariates)) expect_true(all(df$effect[1:10] > 0, na.rm = TRUE))
  }
})
