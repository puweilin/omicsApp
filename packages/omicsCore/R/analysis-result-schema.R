#' Standardized differential analysis result columns
#'
#' All differential analysis backends (DESeq2, edgeR, limma, t-test, lm)
#' return a `data.frame` with these columns so that downstream visualization
#' and integration functions can operate uniformly.
#'
#' @keywords internal
DIFF_RESULT_REQUIRED_COLS <- c(
  "feature_id",
  "feature_symbol",
  "feature_type",
  "omics_type",
  "method",
  "analysis_type",
  "comparison",
  "effect",
  "effect_type",
  "statistic",
  "statistic_type",
  "p_value",
  "adj_p_value",
  "direction",
  "base_mean",
  "model_fit",
  # `NA` until something applies a threshold. Significance is not a
  # property of a feature, it is a property of a (feature, cutoff) pair,
  # and `run_diff()` is never told a cutoff -- p < 0.05, adj.P < 0.05,
  # with or without a fold-change floor are all defensible choices that
  # belong to whoever is asking. `filter_diff_results()` sets it TRUE on
  # the rows it keeps; a figure decides for itself and prints the cut it
  # used. It used to be written FALSE, which asserted "not significant"
  # about features with adj.P near zero and produced monochrome volcano
  # and MA plots that read as findings.
  "is_significant"
)

#' Create an empty diff result template
#'
#' @return Empty `data.frame` following the standardized diff-result schema.
#' @keywords internal
new_diff_result_template <- function() {
  out <- as.data.frame(
    stats::setNames(
      replicate(length(DIFF_RESULT_REQUIRED_COLS), logical(0), simplify = FALSE),
      DIFF_RESULT_REQUIRED_COLS
    )
  )

  out$feature_id <- character(0)
  out$feature_symbol <- character(0)
  out$feature_type <- character(0)
  out$omics_type <- character(0)
  out$method <- character(0)
  out$analysis_type <- character(0)
  out$comparison <- character(0)
  out$effect <- numeric(0)
  out$effect_type <- character(0)
  out$statistic <- numeric(0)
  out$statistic_type <- character(0)
  out$p_value <- numeric(0)
  out$adj_p_value <- numeric(0)
  out$direction <- character(0)
  out$base_mean <- numeric(0)
  out$model_fit <- numeric(0)
  out$is_significant <- logical(0)

  out
}

#' Validate the standard diff result schema
#'
#' @param result_df Data frame to validate.
#'
#' @return Invisibly returns `TRUE` on success.
#' @keywords internal
check_diff_result_schema <- function(result_df) {
  missing_cols <- setdiff(DIFF_RESULT_REQUIRED_COLS, colnames(result_df))
  if (length(missing_cols) > 0) {
    stop("Missing required diff result columns: ", paste(missing_cols, collapse = ", "))
  }
  invisible(TRUE)
}

#' Standardized enrichment result columns
#'
#' All enrichment backends (ORA, GSEA, GSVA) return a `data.frame` with
#' these columns.
#'
#' @keywords internal
ENRICH_RESULT_REQUIRED_COLS <- c(
  "database",
  "result_type",
  "comparison",
  "pathway_id",
  "pathway_name",
  "effect",
  "effect_type",
  "direction",
  "p_value",
  "adj_p_value",
  "q_value",
  "gene_set_size",
  "overlap_size",
  "overlap_features",
  "leading_features",
  "source_label"
)

#' Create an empty enrichment result template
#'
#' @return Empty `data.frame` following the standardized enrichment schema.
#' @keywords internal
new_enrich_result_template <- function() {
  data.frame(
    database = character(0),
    result_type = character(0),
    comparison = character(0),
    pathway_id = character(0),
    pathway_name = character(0),
    effect = numeric(0),
    effect_type = character(0),
    direction = character(0),
    p_value = numeric(0),
    adj_p_value = numeric(0),
    q_value = numeric(0),
    gene_set_size = numeric(0),
    overlap_size = numeric(0),
    overlap_features = character(0),
    leading_features = character(0),
    source_label = character(0),
    stringsAsFactors = FALSE
  )
}

#' Validate the standard enrichment result schema
#'
#' @param result_df Data frame to validate.
#'
#' @return Invisibly returns `TRUE` on success.
#' @keywords internal
check_enrich_result_schema <- function(result_df) {
  missing_cols <- setdiff(ENRICH_RESULT_REQUIRED_COLS, colnames(result_df))
  if (length(missing_cols) > 0) {
    stop("Missing required enrich result columns: ", paste(missing_cols, collapse = ", "))
  }
  invisible(TRUE)
}

#' Standardized integration result columns
#'
#' All [run_integration()] methods return a `data.frame` with these columns
#' so downstream visualisation and reporting code can treat correlation,
#' concordance, and pathway integration uniformly.
#'
#' `effect` carries the method-specific score (Pearson/Spearman r, the
#' difference of effect sizes for concordance, or the combined -log10(p)
#' for ActivePathways) and `effect_type` records which one. `direction`
#' is `"concordant"` / `"discordant"` for concordance, `"positive"` /
#' `"negative"` for correlation, and `"shared"` / `"unique"` for pathway
#' integration. `quadrant` is filled for concordance (e.g. `"up_up"`).
#'
#' @keywords internal
INTEGRATION_RESULT_REQUIRED_COLS <- c(
  "feature_id",
  "feature_symbol",
  "result_type",
  "experiments",
  "comparison",
  "effect",
  "effect_type",
  "statistic",
  "statistic_type",
  "p_value",
  "adj_p_value",
  "direction",
  "quadrant",
  "is_significant",
  "source_label"
)

#' Create an empty integration result template
#'
#' @return Empty `data.frame` following the standardized integration schema.
#' @keywords internal
new_integration_result_template <- function() {
  data.frame(
    feature_id = character(0),
    feature_symbol = character(0),
    result_type = character(0),
    experiments = character(0),
    comparison = character(0),
    effect = numeric(0),
    effect_type = character(0),
    statistic = numeric(0),
    statistic_type = character(0),
    p_value = numeric(0),
    adj_p_value = numeric(0),
    direction = character(0),
    quadrant = character(0),
    is_significant = logical(0),
    source_label = character(0),
    stringsAsFactors = FALSE
  )
}

#' Validate the standard integration result schema
#'
#' @param result_df Data frame to validate.
#'
#' @return Invisibly returns `TRUE` on success.
#' @keywords internal
check_integration_result_schema <- function(result_df) {
  missing_cols <- setdiff(INTEGRATION_RESULT_REQUIRED_COLS, colnames(result_df))
  if (length(missing_cols) > 0) {
    stop("Missing required integration result columns: ",
         paste(missing_cols, collapse = ", "))
  }
  invisible(TRUE)
}
