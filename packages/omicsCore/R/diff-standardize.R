# Standardizers convert backend-native diff outputs (limma, DESeq2, edgeR,
# t-test, lm) into the unified DIFF_RESULT_REQUIRED_COLS schema declared in
# `analysis-result-schema.R`. All standardizers are internal — they are called
# by the corresponding backend in `diff-*.R`, which in turn is dispatched by
# `run_diff()` in `run-diff.R`.
#
# Each one is a mapping from an engine's column names to the standard
# columns, handed to `standardize_diff_columns()`. They used to be eight
# copies of the same join-and-transmute, which is how a new column
# (`signed_stat`) had to be added eight times and could be missed in one.

prep_feature_df_for_standardize <- function(feature_df) {
  if (!"feature_id" %in% colnames(feature_df)) {
    stop("`feature_df` must contain `feature_id`.")
  }
  if (!"feature_symbol" %in% colnames(feature_df)) {
    feature_df$feature_symbol <- feature_df$feature_id
  }
  if (!"feature_type" %in% colnames(feature_df)) {
    feature_df$feature_type <- "feature"
  }
  feature_df
}

# The shared body of every standardizer: check the engine's required
# columns, join the feature annotation, and lay the standard columns out
# in schema order.
#
# `effect`, `statistic`, `signed_stat`, `p`, `padj`, `base_mean` and
# `model_fit` each name a column of the engine's table, or are a
# function of the joined table (for a value computed from several
# columns), or NULL for a column the engine does not have (NA). A named
# column the engine left out is NA as well -- limma without `AveExpr`,
# DESeq2 without `stat` -- as each standardizer used to fill by hand.
# `direction` is "group" (up / down by the effect's sign), "continuous"
# (positive / negative) or "none" (a global test: always "ns").
standardize_diff_columns <- function(raw_df, feature_df, comparison, omics_type,
                                     required, method, analysis_type,
                                     effect, effect_type,
                                     statistic, statistic_type, signed_stat,
                                     p, padj, base_mean = NULL, model_fit = NULL,
                                     direction = c("group", "continuous", "none")) {
  direction <- match.arg(direction)
  check_required_cols(raw_df, required, object_name = "raw_df")

  feature_df <- prep_feature_df_for_standardize(feature_df)
  out <- dplyr::left_join(raw_df, feature_df, by = "feature_id")
  n <- nrow(out)
  pick <- function(x) {
    if (is.null(x)) return(rep(NA_real_, n))
    if (is.function(x)) return(x(out))
    if (!x %in% colnames(out)) return(rep(NA_real_, n))
    out[[x]]
  }
  eff <- pick(effect)
  dir <- switch(
    direction,
    group = dplyr::case_when(eff > 0 ~ "up", eff < 0 ~ "down", TRUE ~ "ns"),
    continuous = dplyr::case_when(eff > 0 ~ "positive", eff < 0 ~ "negative", TRUE ~ "ns"),
    none = rep("ns", n)
  )

  res <- data.frame(
    feature_id = out$feature_id,
    feature_symbol = out$feature_symbol,
    feature_type = out$feature_type,
    omics_type = rep(omics_type, n),
    method = rep(method, n),
    analysis_type = rep(analysis_type, n),
    comparison = rep(comparison, n),
    effect = eff,
    effect_type = rep(effect_type, n),
    statistic = pick(statistic),
    statistic_type = rep(statistic_type, n),
    signed_stat = pick(signed_stat),
    p_value = pick(p),
    adj_p_value = pick(padj),
    direction = dir,
    base_mean = pick(base_mean),
    model_fit = pick(model_fit),
    is_significant = rep(NA, n),
    stringsAsFactors = FALSE
  )

  check_diff_result_schema(res)
  res
}

#' Standardize limma continuous results
#'
#' @param raw_df Raw limma topTable result with `feature_id`, `P.Value`,
#'   `adj.P.Val`, `spearman_rho`, `adj_r_squared`.
#' @param feature_df Feature metadata.
#' @param comparison Continuous variable name.
#' @param omics_type Omics modality label.
#'
#' @return Standardized diff result data frame.
#' @keywords internal
standardize_limma_continuous_results <- function(
  raw_df,
  feature_df,
  comparison,
  omics_type = "proteomics"
) {
  # A spline fit reports one F across its basis columns, which has no
  # direction; a linear fit reports the slope's t.
  is_spline <- "F" %in% colnames(raw_df)
  standardize_diff_columns(
    raw_df, feature_df, comparison, omics_type,
    required = c("feature_id", "P.Value", "adj.P.Val", "spearman_rho", "adj_r_squared"),
    method = "limma",
    analysis_type = if (is_spline) "continuous_spline" else "continuous_linear",
    effect = "spearman_rho", effect_type = "correlation",
    statistic = if (is_spline) "F" else "t",
    statistic_type = if (is_spline) "F" else "t",
    signed_stat = if (is_spline) NULL else "t",
    p = "P.Value", padj = "adj.P.Val",
    base_mean = "AveExpr", model_fit = "adj_r_squared",
    direction = "continuous"
  )
}

#' Standardize limma two-group results
#'
#' @param raw_df Raw limma topTable result with `feature_id`, `logFC`,
#'   `P.Value`, `adj.P.Val`.
#' @param feature_df Feature metadata.
#' @param comparison Comparison label.
#' @param omics_type Omics modality label.
#'
#' @return Standardized diff result data frame.
#' @keywords internal
standardize_limma_group_results <- function(
  raw_df,
  feature_df,
  comparison,
  omics_type = "proteomics"
) {
  standardize_diff_columns(
    raw_df, feature_df, comparison, omics_type,
    required = c("feature_id", "logFC", "P.Value", "adj.P.Val"),
    method = "limma", analysis_type = "group",
    effect = "logFC", effect_type = "log2FC",
    statistic = "t", statistic_type = "t", signed_stat = "t",
    p = "P.Value", padj = "adj.P.Val",
    base_mean = "AveExpr",
    direction = "group"
  )
}

#' Standardize edgeR two-group results
#'
#' @param raw_df Raw edgeR topTags result with `feature_id`, `logFC`, `PValue`,
#'   `FDR`.
#' @param feature_df Feature metadata.
#' @param comparison Comparison label.
#' @param omics_type Omics modality label.
#'
#' @return Standardized diff result data frame.
#' @keywords internal
standardize_edger_group_results <- function(
  raw_df,
  feature_df,
  comparison,
  omics_type = "rnaseq"
) {
  standardize_diff_columns(
    raw_df, feature_df, comparison, omics_type,
    required = c("feature_id", "logFC", "PValue", "FDR"),
    method = "edger", analysis_type = "group",
    effect = "logFC", effect_type = "log2FC",
    statistic = "F", statistic_type = "F",
    # A QL F on one degree of freedom is a squared t; its root, with
    # the sign of the fold change, is on the other engines' scale.
    signed_stat = function(out) {
      f <- if ("F" %in% colnames(out)) out$F else NA_real_
      sign(out$logFC) * sqrt(f)
    },
    p = "PValue", padj = "FDR",
    base_mean = "logCPM",
    direction = "group"
  )
}

#' Standardize DESeq2 two-group results
#'
#' @param raw_df Raw DESeq2 result with `feature_id`, `log2FoldChange`,
#'   `pvalue`, `padj`.
#' @param feature_df Feature metadata.
#' @param comparison Comparison label.
#' @param omics_type Omics modality label.
#'
#' @return Standardized diff result data frame.
#' @keywords internal
standardize_deseq2_group_results <- function(
  raw_df,
  feature_df,
  comparison,
  omics_type = "rnaseq"
) {
  standardize_diff_columns(
    raw_df, feature_df, comparison, omics_type,
    required = c("feature_id", "log2FoldChange", "pvalue", "padj"),
    method = "deseq2", analysis_type = "group",
    effect = "log2FoldChange", effect_type = "log2FC",
    statistic = "stat", statistic_type = "wald", signed_stat = "stat",
    p = "pvalue", padj = "padj",
    base_mean = "baseMean",
    direction = "group"
  )
}

#' Standardize DESeq2 continuous results
#'
#' @param raw_df Raw DESeq2 result with `feature_id`, `log2FoldChange`,
#'   `pvalue`, `padj`.
#' @param feature_df Feature metadata.
#' @param comparison Continuous variable name.
#' @param omics_type Omics modality label.
#'
#' @return Standardized diff result data frame.
#' @keywords internal
standardize_deseq2_continuous_results <- function(
  raw_df,
  feature_df,
  comparison,
  omics_type = "rnaseq"
) {
  standardize_diff_columns(
    raw_df, feature_df, comparison, omics_type,
    required = c("feature_id", "log2FoldChange", "pvalue", "padj"),
    method = "deseq2", analysis_type = "continuous_linear",
    effect = "log2FoldChange", effect_type = "log2FC_per_unit",
    statistic = "stat", statistic_type = "wald", signed_stat = "stat",
    p = "pvalue", padj = "padj",
    base_mean = "baseMean",
    direction = "continuous"
  )
}

#' Standardize per-feature t-test group results
#'
#' @param raw_df Per-feature t-test table with `feature_id`, `mean_diff`,
#'   `t_stat`, `p_value`, `adj_p_value`, optionally `mean_ctrl`, `mean_case`.
#' @param feature_df Feature metadata.
#' @param comparison Comparison label.
#' @param omics_type Omics modality label.
#'
#' @return Standardized diff result data frame.
#' @keywords internal
standardize_ttest_group_results <- function(
  raw_df,
  feature_df,
  comparison,
  omics_type = "proteomics"
) {
  has_means <- all(c("mean_ctrl", "mean_case") %in% colnames(raw_df))
  standardize_diff_columns(
    raw_df, feature_df, comparison, omics_type,
    required = c("feature_id", "mean_diff", "t_stat", "p_value", "adj_p_value"),
    method = "ttest", analysis_type = "group",
    effect = "mean_diff", effect_type = "mean_diff",
    statistic = "t_stat", statistic_type = "t", signed_stat = "t_stat",
    p = "p_value", padj = "adj_p_value",
    base_mean = if (has_means) function(out) (out$mean_ctrl + out$mean_case) / 2,
    direction = "group"
  )
}

#' Standardize per-feature lm group results
#'
#' @param raw_df Per-feature lm table with `feature_id`, `beta`, `t_stat`,
#'   `p_value`, `adj_p_value`, `adj_r_squared`, `base_mean`.
#' @param feature_df Feature metadata.
#' @param comparison Comparison label.
#' @param omics_type Omics modality label.
#'
#' @return Standardized diff result data frame.
#' @keywords internal
standardize_lm_group_results <- function(
  raw_df,
  feature_df,
  comparison,
  omics_type = "proteomics"
) {
  standardize_diff_columns(
    raw_df, feature_df, comparison, omics_type,
    required = c("feature_id", "beta", "t_stat", "p_value", "adj_p_value",
                 "adj_r_squared", "base_mean"),
    method = "lm", analysis_type = "group",
    effect = "beta", effect_type = "beta",
    statistic = "t_stat", statistic_type = "t", signed_stat = "t_stat",
    p = "p_value", padj = "adj_p_value",
    base_mean = "base_mean", model_fit = "adj_r_squared",
    direction = "group"
  )
}

#' Standardize per-feature lm continuous results
#'
#' @param raw_df Per-feature lm table with `feature_id`, `spearman_rho`,
#'   `t_stat`, `p_value`, `adj_p_value`, `adj_r_squared`, `base_mean`.
#' @param feature_df Feature metadata.
#' @param comparison Continuous variable name.
#' @param omics_type Omics modality label.
#'
#' @return Standardized diff result data frame.
#' @keywords internal
standardize_lm_continuous_results <- function(
  raw_df,
  feature_df,
  comparison,
  omics_type = "proteomics"
) {
  standardize_diff_columns(
    raw_df, feature_df, comparison, omics_type,
    required = c("feature_id", "spearman_rho", "t_stat", "p_value",
                 "adj_p_value", "adj_r_squared", "base_mean"),
    method = "lm", analysis_type = "continuous_linear",
    effect = "spearman_rho", effect_type = "correlation",
    statistic = "t_stat", statistic_type = "t", signed_stat = "t_stat",
    p = "p_value", padj = "adj_p_value",
    base_mean = "base_mean", model_fit = "adj_r_squared",
    direction = "continuous"
  )
}

# A global test across groups (limma's moderated F, edgeR's QL F,
# DESeq2's LRT): the test statistic stands in for the effect, and there
# is no direction and so no signed statistic.
standardize_global_test_results <- function(raw_df, feature_df, comparison, omics_type,
                                            method, stat, stat_type, p, padj,
                                            base_mean) {
  standardize_diff_columns(
    raw_df, feature_df, comparison, omics_type,
    required = c("feature_id", stat, p, padj),
    method = method, analysis_type = "anova",
    effect = stat, effect_type = paste0(stat_type, "_statistic"),
    statistic = stat, statistic_type = stat_type, signed_stat = NULL,
    p = p, padj = padj, base_mean = base_mean,
    direction = "none"
  )
}
