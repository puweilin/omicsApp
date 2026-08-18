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
              assay_type = "intensity")
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
