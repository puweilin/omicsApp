# Boundary shapes the plot layer has to survive. A view renders whatever
# the analysis produced, and an analysis can legitimately produce
# nothing: every feature filtered out by QC, a contrast with no
# variance, a database that matched no pathway. Erroring there takes the
# whole view down instead of showing an empty panel.

pb_input <- function(n_features = 6L, values = NULL) {
  n <- n_features
  mat <- matrix(values %||% as.numeric(seq_len(n * 4L)), nrow = n,
                dimnames = list(paste0("g", seq_len(n)), paste0("s", 1:4)))
  meta <- data.frame(group = c("A", "A", "B", "B"),
                     row.names = paste0("s", 1:4))
  feat <- data.frame(feature_id = paste0("g", seq_len(n)),
                     feature_symbol = paste0("SYM", seq_len(n)))
  omics_input(mat, meta, feat, omics_type = "proteomics",
              assay_type = "normalized_intensity")
}

pb_diff <- function(input = pb_input()) {
  run_diff(input, method = "ttest", analysis_type = "group",
           group_col = "group", control_group = "A", case_group = "B")
}

# ---- empty and near-empty results -------------------------------------

test_that("a diff result with no rows still draws a volcano", {
  b <- pb_diff()
  b$results$diff_result_df <- b$results$diff_result_df[0, , drop = FALSE]
  # The axis label is read off the first row; without a guard this was
  # a subscript error rather than an empty panel.
  expect_s3_class(plot_volcano(b), "ggplot")
  expect_s3_class(plot_volcano(b, p_threshold = 0.05), "ggplot")
})

test_that("a single-feature result draws", {
  b <- pb_diff()
  b$results$diff_result_df <- b$results$diff_result_df[1, , drop = FALSE]
  expect_s3_class(plot_volcano(b), "ggplot")
})

test_that("an all-NA effect column does not stop the volcano", {
  b <- pb_diff()
  b$results$diff_result_df$effect <- NA_real_
  expect_s3_class(plot_volcano(b), "ggplot")
})

test_that("an all-NA p column does not stop the volcano", {
  b <- pb_diff()
  b$results$diff_result_df$adj_p_value <- NA_real_
  expect_s3_class(plot_volcano(b, p_threshold = 0.05), "ggplot")
})

# ---- threshold handling ------------------------------------------------

test_that("explicit NULL thresholds mean 'draw no rule', not an error", {
  b <- pb_diff()
  # An earlier version reached log10(NULL) and died.
  expect_s3_class(plot_volcano(b, p_threshold = NULL), "ggplot")
  expect_s3_class(plot_volcano(b, effect_threshold = NULL), "ggplot")
  expect_s3_class(plot_volcano(b, p_threshold = NULL,
                               effect_threshold = NULL), "ggplot")
})

test_that("NULL thresholds fall back to the analysis's own mask", {
  b <- pb_diff()
  b$results$diff_result_df$is_significant <- rep(c(TRUE, FALSE), length.out =
                                                   nrow(b$results$diff_result_df))
  built <- ggplot2::ggplot_build(
    plot_volcano(b, p_threshold = NULL, effect_threshold = NULL))$data[[1]]
  # Both classes present: the stored mask was honoured rather than every
  # point being swept into one colour.
  expect_equal(length(unique(built$colour)), 2L)
})

test_that("a threshold nothing passes colours nothing significant", {
  b <- pb_diff()
  built <- ggplot2::ggplot_build(
    plot_volcano(b, p_threshold = 1e-300, effect_threshold = 1e6))$data[[1]]
  expect_equal(length(unique(built$colour)), 1L)
})

# ---- labels ------------------------------------------------------------

test_that("top_n larger than the result is not an error", {
  b <- pb_diff()
  expect_s3_class(plot_volcano(b, top_n = 10000L), "ggplot")
})

test_that("top_n of zero draws no labels", {
  b <- pb_diff()
  expect_s3_class(plot_volcano(b, top_n = 0L), "ggplot")
})

test_that("label_features naming nothing present is not an error", {
  b <- pb_diff()
  expect_s3_class(plot_volcano(b, label_features = "NOT_A_GENE"), "ggplot")
})

# ---- enrichment --------------------------------------------------------

test_that("an enrichment result without an effect column falls back", {
  df <- data.frame(
    database = "hallmark", result_type = "ora", comparison = "B_vs_A",
    pathway_id = paste0("P", 1:3), pathway_name = paste0("Pathway ", 1:3),
    effect = NA_real_, effect_type = "log2_odds",
    direction = "up", p_value = c(0.001, 0.01, 0.05),
    adj_p_value = c(0.002, 0.02, 0.1), q_value = c(0.002, 0.02, 0.1),
    gene_set_size = c(100L, 80L, 60L), overlap_size = c(10L, 8L, 6L),
    overlap_features = NA_character_, leading_features = NA_character_,
    source_label = "test", stringsAsFactors = FALSE
  )
  b <- new_analysis_bundle("run_enrichment", list(omics_type = "proteomics"),
                           params = list(type = "ora"),
                           results = list(enrich_result_df = df))
  # An ORA result can carry no usable effect; the dot plot has to fall
  # back to overlap size rather than plotting an all-NA axis.
  expect_s3_class(plot_enrichment(b, view = "dot"), "ggplot")
})

test_that("an enrichment result with no rows returns an empty plot", {
  df <- new_enrich_result_template()
  b <- new_analysis_bundle("run_enrichment", list(omics_type = "proteomics"),
                           params = list(type = "ora"),
                           results = list(enrich_result_df = df))
  expect_s3_class(plot_enrichment(b, view = "dot"), "ggplot")
})

# ---- integration -------------------------------------------------------

test_that("the integration views reject a method they cannot draw", {
  df <- new_integration_result_template()
  b <- new_analysis_bundle("run_integration",
                           list(omics_type = c("a", "b")),
                           params = list(method = "correlation"),
                           results = list(integration_df = df))
  # Failing loudly here is right: a quadrant plot of a correlation run
  # would be a picture of a column that does not exist.
  expect_error(plot_integration(b, view = "quadrant"), "concordance")
  expect_error(plot_integration(b, view = "effect_pair"), "concordance")
})

test_that("an integration bundle with no rows returns an empty plot", {
  df <- new_integration_result_template()
  b <- new_analysis_bundle("run_integration",
                           list(omics_type = c("a", "b")),
                           params = list(method = "concordance",
                                         experiments = c("a", "b")),
                           results = list(integration_df = df))
  expect_s3_class(plot_integration(b, view = "effect_pair"), "ggplot")
})

# ---- where "significant" comes from -----------------------------------
#
# Every run_diff() backend writes is_significant = FALSE and leaves the
# judgement to whoever filters the result. A volcano that deferred to
# that column drew one colour over data where a fifth of the features
# cleared 0.05 -- and in a report that reads as a finding rather than as
# a missing step.

pb_split_diff <- function() {
  set.seed(42)
  n_sig <- 20L; n_null <- 80L; ns <- 20L
  mat <- rbind(
    matrix(c(rnorm(n_sig * ns / 2, 8, 1), rnorm(n_sig * ns / 2, 11, 1)),
           nrow = n_sig),
    matrix(rnorm(n_null * ns, 8, 1), nrow = n_null)
  )
  dimnames(mat) <- list(paste0("g", seq_len(n_sig + n_null)),
                        paste0("s", seq_len(ns)))
  meta <- data.frame(group = rep(c("A", "B"), each = ns / 2),
                     row.names = colnames(mat))
  feat <- data.frame(feature_id = rownames(mat),
                     feature_symbol = rownames(mat))
  inp <- omics_input(mat, meta, feat, omics_type = "proteomics",
                     assay_type = "normalized_intensity")
  run_diff(inp, method = "limma", analysis_type = "group",
           group_col = "group", control_group = "A", case_group = "B")
}

sig_counts <- function(p) {
  table(factor(ggplot2::ggplot_build(p)$plot$data$.sig,
               levels = c("ns", "significant")))
}

test_that("the default cut separates the features that clear it", {
  skip_if_not_installed("limma")
  b <- pb_split_diff()
  df <- b$results$diff_result_df
  # The premise of the regression: the bundle answers nothing, so the
  # colouring has to come from the threshold this call was given.
  expect_true(all(is.na(df$is_significant)))
  expect_gt(sum(df$adj_p_value < 0.05, na.rm = TRUE), 0L)

  counts <- sig_counts(plot_volcano(b))
  expect_equal(unname(counts[["significant"]]),
               sum(df$adj_p_value < 0.05, na.rm = TRUE))
  expect_gt(unname(counts[["ns"]]), 0L)
})

test_that("an effect threshold narrows the set further", {
  skip_if_not_installed("limma")
  b <- pb_split_diff()
  loose <- sig_counts(plot_volcano(b))
  tight <- sig_counts(plot_volcano(b, effect_threshold = 10))
  expect_lt(unname(tight[["significant"]]), unname(loose[["significant"]]))
})

test_that("both thresholds NULL defers to the stored column", {
  skip_if_not_installed("limma")
  b <- pb_split_diff()
  counts <- sig_counts(plot_volcano(b, p_threshold = NULL,
                                    effect_threshold = NULL))
  # Which run_diff() leaves FALSE -- "make no distinction" is a request
  # the caller is allowed to make, and this is what it looks like.
  expect_equal(unname(counts[["significant"]]), 0L)
})

test_that("the figure states the cut it was drawn at", {
  skip_if_not_installed("limma")
  b <- pb_split_diff()
  cap <- function(p) ggplot2::ggplot_build(p)$plot$labels$caption
  # A screenshot outlives the session that produced it, so the cut has
  # to travel with the picture.
  expect_match(cap(plot_volcano(b)), "adj_p_value < 0.05", fixed = TRUE)
  expect_match(cap(plot_volcano(b, effect_threshold = 1.5)),
               "|effect| > 1.5", fixed = TRUE)
  expect_match(cap(plot_volcano(b, p_threshold = NULL,
                                effect_threshold = NULL)),
               "as recorded", fixed = TRUE)
})

test_that("run_diff leaves significance unanswered rather than answering FALSE", {
  skip_if_not_installed("limma")
  b <- pb_split_diff()
  df <- b$results$diff_result_df
  # NA is the truth: no cutoff was supplied, so the question is open.
  # FALSE was a claim -- "this feature is not significant" -- made about
  # features with adj.P near zero, and two plot functions believed it.
  expect_true(all(is.na(df$is_significant)))
  expect_gt(sum(df$adj_p_value < 0.05, na.rm = TRUE), 0L)
})

test_that("filter_diff_results still answers TRUE on what it keeps", {
  skip_if_not_installed("limma")
  df <- pb_split_diff()$results$diff_result_df
  filt <- filter_diff_results(df, p_cutoff = 0.05)
  # Applying a cutoff is what turns the open question into an answer.
  expect_true(all(filt$is_significant))
  expect_lt(nrow(filt), nrow(df))
})

test_that("the MA plot decides significance the same way the volcano does", {
  skip_if_not_installed("limma")
  b <- pb_split_diff()
  n_sig <- sum(b$results$diff_result_df$adj_p_value < 0.05, na.rm = TRUE)
  counts <- table(factor(ggplot2::ggplot_build(plot_ma(b))$plot$data$.sig,
                         levels = c("ns", "significant")))
  # It read the same always-FALSE column and drew the same single
  # colour; fixing one and not the other would have left the pair
  # disagreeing about the same data.
  expect_equal(unname(counts[["significant"]]), n_sig)
  expect_match(ggplot2::ggplot_build(plot_ma(b))$plot$labels$caption,
               "adj_p_value < 0.05", fixed = TRUE)
})

test_that("the MA plot honours a supplied threshold", {
  skip_if_not_installed("limma")
  b <- pb_split_diff()
  loose <- table(factor(ggplot2::ggplot_build(plot_ma(b))$plot$data$.sig,
                        levels = c("ns", "significant")))
  tight <- table(factor(ggplot2::ggplot_build(
    plot_ma(b, effect_threshold = 10))$plot$data$.sig,
    levels = c("ns", "significant")))
  expect_lt(unname(tight[["significant"]]), unname(loose[["significant"]]))
})
