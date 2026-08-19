# Standardizers convert backend-native diff outputs (limma, DESeq2, edgeR,
# t-test, lm) into the unified DIFF_RESULT_REQUIRED_COLS schema declared in
# `analysis-result-schema.R`. All standardizers are internal — they are called
# by the corresponding backend in `diff-*.R`, which in turn is dispatched by
# `run_diff()` in `run-diff.R`.

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
  required_cols <- c("feature_id", "P.Value", "adj.P.Val", "spearman_rho", "adj_r_squared")
  check_required_cols(raw_df, required_cols, object_name = "raw_df")

  feature_df <- prep_feature_df_for_standardize(feature_df)
  out <- dplyr::left_join(raw_df, feature_df, by = "feature_id")

  if (!"AveExpr" %in% colnames(out)) out$AveExpr <- NA_real_
  if (!"t" %in% colnames(out)) out$t <- NA_real_

  is_spline <- "F" %in% colnames(raw_df)
  out <- out |>
    dplyr::transmute(
      feature_id = .data$feature_id,
      feature_symbol = .data$feature_symbol,
      feature_type = .data$feature_type,
      omics_type = omics_type,
      method = "limma",
      analysis_type = if (is_spline) "continuous_spline" else "continuous_linear",
      comparison = comparison,
      effect = .data$spearman_rho,
      effect_type = "correlation",
      statistic = if (is_spline) .data$F else .data$t,
      statistic_type = if (is_spline) "F" else "t",
      p_value = .data$P.Value,
      adj_p_value = .data$adj.P.Val,
      direction = dplyr::case_when(
        .data$spearman_rho > 0 ~ "positive",
        .data$spearman_rho < 0 ~ "negative",
        TRUE ~ "ns"
      ),
      base_mean = .data$AveExpr,
      model_fit = .data$adj_r_squared,
      is_significant = NA
    )

  check_diff_result_schema(out)
  out
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
  required_cols <- c("feature_id", "logFC", "P.Value", "adj.P.Val")
  check_required_cols(raw_df, required_cols, object_name = "raw_df")

  feature_df <- prep_feature_df_for_standardize(feature_df)
  out <- dplyr::left_join(raw_df, feature_df, by = "feature_id")

  if (!"AveExpr" %in% colnames(out)) out$AveExpr <- NA_real_
  if (!"t" %in% colnames(out)) out$t <- NA_real_

  out <- out |>
    dplyr::transmute(
      feature_id = .data$feature_id,
      feature_symbol = .data$feature_symbol,
      feature_type = .data$feature_type,
      omics_type = omics_type,
      method = "limma",
      analysis_type = "group",
      comparison = comparison,
      effect = .data$logFC,
      effect_type = "log2FC",
      statistic = .data$t,
      statistic_type = "t",
      p_value = .data$P.Value,
      adj_p_value = .data$adj.P.Val,
      direction = dplyr::case_when(
        .data$logFC > 0 ~ "up",
        .data$logFC < 0 ~ "down",
        TRUE ~ "ns"
      ),
      base_mean = .data$AveExpr,
      model_fit = NA_real_,
      is_significant = NA
    )

  check_diff_result_schema(out)
  out
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
  required_cols <- c("feature_id", "logFC", "PValue", "FDR")
  check_required_cols(raw_df, required_cols, object_name = "raw_df")

  feature_df <- prep_feature_df_for_standardize(feature_df)
  out <- dplyr::left_join(raw_df, feature_df, by = "feature_id")

  if (!"logCPM" %in% colnames(out)) out$logCPM <- NA_real_
  if (!"F" %in% colnames(out)) out$F <- NA_real_

  out <- out |>
    dplyr::transmute(
      feature_id = .data$feature_id,
      feature_symbol = .data$feature_symbol,
      feature_type = .data$feature_type,
      omics_type = omics_type,
      method = "edger",
      analysis_type = "group",
      comparison = comparison,
      effect = .data$logFC,
      effect_type = "log2FC",
      statistic = .data$F,
      statistic_type = "F",
      p_value = .data$PValue,
      adj_p_value = .data$FDR,
      direction = dplyr::case_when(
        .data$logFC > 0 ~ "up",
        .data$logFC < 0 ~ "down",
        TRUE ~ "ns"
      ),
      base_mean = .data$logCPM,
      model_fit = NA_real_,
      is_significant = NA
    )

  check_diff_result_schema(out)
  out
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
  required_cols <- c("feature_id", "log2FoldChange", "pvalue", "padj")
  check_required_cols(raw_df, required_cols, object_name = "raw_df")

  feature_df <- prep_feature_df_for_standardize(feature_df)
  out <- dplyr::left_join(raw_df, feature_df, by = "feature_id")

  if (!"baseMean" %in% colnames(out)) out$baseMean <- NA_real_
  if (!"stat" %in% colnames(out)) out$stat <- NA_real_

  out <- out |>
    dplyr::transmute(
      feature_id = .data$feature_id,
      feature_symbol = .data$feature_symbol,
      feature_type = .data$feature_type,
      omics_type = omics_type,
      method = "deseq2",
      analysis_type = "group",
      comparison = comparison,
      effect = .data$log2FoldChange,
      effect_type = "log2FC",
      statistic = .data$stat,
      statistic_type = "wald",
      p_value = .data$pvalue,
      adj_p_value = .data$padj,
      direction = dplyr::case_when(
        .data$log2FoldChange > 0 ~ "up",
        .data$log2FoldChange < 0 ~ "down",
        TRUE ~ "ns"
      ),
      base_mean = .data$baseMean,
      model_fit = NA_real_,
      is_significant = NA
    )

  check_diff_result_schema(out)
  out
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
  required_cols <- c("feature_id", "log2FoldChange", "pvalue", "padj")
  check_required_cols(raw_df, required_cols, object_name = "raw_df")

  feature_df <- prep_feature_df_for_standardize(feature_df)
  out <- dplyr::left_join(raw_df, feature_df, by = "feature_id")

  if (!"baseMean" %in% colnames(out)) out$baseMean <- NA_real_
  if (!"stat" %in% colnames(out)) out$stat <- NA_real_

  out <- out |>
    dplyr::transmute(
      feature_id = .data$feature_id,
      feature_symbol = .data$feature_symbol,
      feature_type = .data$feature_type,
      omics_type = omics_type,
      method = "deseq2",
      analysis_type = "continuous_linear",
      comparison = comparison,
      effect = .data$log2FoldChange,
      effect_type = "log2FC_per_unit",
      statistic = .data$stat,
      statistic_type = "wald",
      p_value = .data$pvalue,
      adj_p_value = .data$padj,
      direction = dplyr::case_when(
        .data$log2FoldChange > 0 ~ "positive",
        .data$log2FoldChange < 0 ~ "negative",
        TRUE ~ "ns"
      ),
      base_mean = .data$baseMean,
      model_fit = NA_real_,
      is_significant = NA
    )

  check_diff_result_schema(out)
  out
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
  required_cols <- c("feature_id", "mean_diff", "t_stat", "p_value", "adj_p_value")
  check_required_cols(raw_df, required_cols, object_name = "raw_df")

  feature_df <- prep_feature_df_for_standardize(feature_df)
  out <- dplyr::left_join(raw_df, feature_df, by = "feature_id")

  has_means <- all(c("mean_ctrl", "mean_case") %in% colnames(raw_df))
  out <- out |>
    dplyr::transmute(
      feature_id = .data$feature_id,
      feature_symbol = .data$feature_symbol,
      feature_type = .data$feature_type,
      omics_type = omics_type,
      method = "ttest",
      analysis_type = "group",
      comparison = comparison,
      effect = .data$mean_diff,
      effect_type = "mean_diff",
      statistic = .data$t_stat,
      statistic_type = "t",
      p_value = .data$p_value,
      adj_p_value = .data$adj_p_value,
      direction = dplyr::case_when(
        .data$mean_diff > 0 ~ "up",
        .data$mean_diff < 0 ~ "down",
        TRUE ~ "ns"
      ),
      base_mean = if (has_means) (.data$mean_ctrl + .data$mean_case) / 2 else NA_real_,
      model_fit = NA_real_,
      is_significant = NA
    )

  check_diff_result_schema(out)
  out
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
  required_cols <- c("feature_id", "beta", "t_stat", "p_value", "adj_p_value",
                     "adj_r_squared", "base_mean")
  check_required_cols(raw_df, required_cols, object_name = "raw_df")

  feature_df <- prep_feature_df_for_standardize(feature_df)
  out <- dplyr::left_join(raw_df, feature_df, by = "feature_id")

  out <- out |>
    dplyr::transmute(
      feature_id = .data$feature_id,
      feature_symbol = .data$feature_symbol,
      feature_type = .data$feature_type,
      omics_type = omics_type,
      method = "lm",
      analysis_type = "group",
      comparison = comparison,
      effect = .data$beta,
      effect_type = "beta",
      statistic = .data$t_stat,
      statistic_type = "t",
      p_value = .data$p_value,
      adj_p_value = .data$adj_p_value,
      direction = dplyr::case_when(
        .data$beta > 0 ~ "up",
        .data$beta < 0 ~ "down",
        TRUE ~ "ns"
      ),
      base_mean = .data$base_mean,
      model_fit = .data$adj_r_squared,
      is_significant = NA
    )

  check_diff_result_schema(out)
  out
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
  required_cols <- c("feature_id", "spearman_rho", "t_stat", "p_value",
                     "adj_p_value", "adj_r_squared", "base_mean")
  check_required_cols(raw_df, required_cols, object_name = "raw_df")

  feature_df <- prep_feature_df_for_standardize(feature_df)
  out <- dplyr::left_join(raw_df, feature_df, by = "feature_id")

  out <- out |>
    dplyr::transmute(
      feature_id = .data$feature_id,
      feature_symbol = .data$feature_symbol,
      feature_type = .data$feature_type,
      omics_type = omics_type,
      method = "lm",
      analysis_type = "continuous_linear",
      comparison = comparison,
      effect = .data$spearman_rho,
      effect_type = "correlation",
      statistic = .data$t_stat,
      statistic_type = "t",
      p_value = .data$p_value,
      adj_p_value = .data$adj_p_value,
      direction = dplyr::case_when(
        .data$spearman_rho > 0 ~ "positive",
        .data$spearman_rho < 0 ~ "negative",
        TRUE ~ "ns"
      ),
      base_mean = .data$base_mean,
      model_fit = .data$adj_r_squared,
      is_significant = NA
    )

  check_diff_result_schema(out)
  out
}
