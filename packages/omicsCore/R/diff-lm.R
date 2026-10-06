#' Per-feature linear regression for a two-group comparison
#'
#' Fits `expr ~ group_col [+ covariates]` per feature using `stats::lm`,
#' extracting the case-group coefficient as the effect estimate. RNA-seq raw
#' counts are automatically log2(x+1) transformed. Internal — call via
#' [run_diff()].
#'
#' @param input A validated `omics_input`.
#' @param group_col Group column in sample metadata.
#' @param control_group Control-group label.
#' @param case_group Case-group label.
#' @param covariates Optional character vector of covariate column names.
#'
#' @return List with `results_raw`, `results_std`, `model_object` (`NULL`),
#'   and `analysis_info`.
#' @keywords internal
run_lm_group <- function(
  input,
  group_col,
  control_group,
  case_group,
  covariates = NULL
) {
  validate_omics_input(input)

  expr_mat <- input$expr_mat
  meta_df <- input$meta_df
  feature_df <- input$feature_df

  if (identical(input$assay_type, "raw_count")) {
    expr_mat <- log2(expr_mat + 1)
  }
  if (!group_col %in% colnames(meta_df)) {
    stop("`group_col` not found in `meta_df`: ", group_col)
  }

  meta_df[[group_col]] <- factor(meta_df[[group_col]])
  target_meta <- meta_df[meta_df[[group_col]] %in% c(control_group, case_group), , drop = FALSE]
  target_meta[[group_col]] <- factor(target_meta[[group_col]], levels = c(control_group, case_group))

  if (!is.null(covariates)) {
    missing_cov <- setdiff(covariates, colnames(target_meta))
    if (length(missing_cov) > 0L) {
      stop("Missing covariates: ", paste(missing_cov, collapse = ", "))
    }
  }

  keep_samples <- rownames(target_meta)
  expr_sub <- expr_mat[, keep_samples, drop = FALSE]

  # The model sees placeholder column names, so a column called
  # "Treatment Group" or "Body mass" cannot break the formula.
  cov_names <- if (length(covariates)) paste0(".cov", seq_along(covariates)) else character(0)
  design_df <- stats::setNames(
    data.frame(target_meta[, c(group_col, covariates), drop = FALSE],
               check.names = FALSE, stringsAsFactors = FALSE),
    c(".grp", cov_names))
  design_df <- droplevels(design_df)
  rhs <- paste(c(".grp", cov_names), collapse = " + ")
  formula_obj <- stats::as.formula(paste("y ~", rhs))
  coef_name <- paste0(".grp", case_group)

  feature_ids <- rownames(expr_sub)
  n_features <- length(feature_ids)

  # One QR per pattern of missing values instead of an lm() per feature
  # (lm_rows(), fast-rows.R): the same estimates, 9x faster with missing
  # values and ~200x without.
  X <- stats::model.matrix(stats::delete.response(stats::terms(formula_obj)), data = design_df)
  fit <- lm_rows(expr_sub, X, coef_name, fac = design_factors(design_df))
  beta <- fit$beta
  t_stat <- fit$t
  p_value <- fit$p
  adj_r_squared <- fit$adj_r2
  base_mean <- unname(rowMeans(expr_sub, na.rm = TRUE))

  raw_df <- data.frame(
    feature_id = feature_ids,
    beta = beta,
    t_stat = t_stat,
    p_value = p_value,
    adj_p_value = stats::p.adjust(p_value, method = "BH"),
    adj_r_squared = adj_r_squared,
    base_mean = base_mean,
    stringsAsFactors = FALSE
  )

  comparison <- paste0(case_group, "_vs_", control_group)
  results_std <- standardize_lm_group_results(
    raw_df = raw_df,
    feature_df = feature_df,
    comparison = comparison,
    omics_type = input$omics_type
  )

  list(
    results_raw = raw_df,
    results_std = results_std,
    model_object = NULL,
    analysis_info = list(
      omics_type = input$omics_type,
      method = "lm",
      analysis_type = "group",
      comparison = comparison,
      covariates = covariates
    )
  )
}

#' Per-feature linear regression against a continuous variable
#'
#' Fits `expr ~ continuous_col [+ covariates]` per feature, returning the
#' coefficient on `continuous_col` as the effect, plus a (partial) Spearman
#' rank correlation. RNA-seq raw counts are automatically log2(x+1)
#' transformed. Internal — call via [run_diff()].
#'
#' @param input A validated `omics_input`.
#' @param continuous_col Continuous metadata column name.
#' @param covariates Optional character vector of covariate column names.
#'
#' @return List with `results_raw`, `results_std`, `model_object` (`NULL`),
#'   and `analysis_info`.
#' @keywords internal
run_lm_continuous <- function(
  input,
  continuous_col,
  covariates = NULL
) {
  validate_omics_input(input)

  expr_mat <- input$expr_mat
  meta_df <- input$meta_df
  feature_df <- input$feature_df

  if (identical(input$assay_type, "raw_count")) {
    expr_mat <- log2(expr_mat + 1)
  }
  if (!continuous_col %in% colnames(meta_df)) {
    stop("`continuous_col` not found in `meta_df`: ", continuous_col)
  }

  cont_vals <- coerce_continuous_col(meta_df[[continuous_col]], continuous_col)
  meta_df[[continuous_col]] <- cont_vals

  if (!is.null(covariates)) {
    missing_cov <- setdiff(covariates, colnames(meta_df))
    if (length(missing_cov) > 0L) {
      stop("Missing covariates: ", paste(missing_cov, collapse = ", "))
    }
  }

  cov_names <- if (length(covariates)) paste0(".cov", seq_along(covariates)) else character(0)
  design_df <- stats::setNames(
    data.frame(meta_df[, c(continuous_col, covariates), drop = FALSE],
               check.names = FALSE, stringsAsFactors = FALSE),
    c(".cont", cov_names))
  rhs <- paste(c(".cont", cov_names), collapse = " + ")
  formula_obj <- stats::as.formula(paste("y ~", rhs))
  adjustment_terms <- cov_names

  feature_ids <- rownames(expr_mat)
  n_features <- length(feature_ids)

  X <- stats::model.matrix(stats::delete.response(stats::terms(formula_obj)), data = design_df)
  fit <- lm_rows(expr_mat, X, ".cont", fac = design_factors(design_df))
  beta <- fit$beta
  t_stat <- fit$t
  p_value <- fit$p
  adj_r_squared <- fit$adj_r2
  base_mean <- unname(rowMeans(expr_mat, na.rm = TRUE))
  if (length(adjustment_terms) == 0L) {
    spearman_rho <- spearman_rows(expr_mat, cont_vals)
  } else {
    # Partial rank correlation: both sides residualised on the
    # adjustment terms, on the features observed in every sample.
    Xc <- stats::model.matrix(stats::as.formula(
      paste("~", paste(adjustment_terms, collapse = " + "))), data = design_df)
    cont_res <- unname(stats::lm.fit(Xc, cont_vals)$residuals)
    full <- rowSums(is.na(expr_mat)) == 0L
    spearman_rho <- rep(NA_real_, nrow(expr_mat))
    if (any(full)) {
      spearman_rho[full] <- spearman_rows(
        t(qr.resid(qr(Xc), t(expr_mat[full, , drop = FALSE]))), cont_res)
    }
  }

  raw_df <- data.frame(
    feature_id = feature_ids,
    beta = beta,
    t_stat = t_stat,
    p_value = p_value,
    adj_p_value = stats::p.adjust(p_value, method = "BH"),
    adj_r_squared = adj_r_squared,
    spearman_rho = spearman_rho,
    base_mean = base_mean,
    stringsAsFactors = FALSE
  )

  results_std <- standardize_lm_continuous_results(
    raw_df = raw_df,
    feature_df = feature_df,
    comparison = continuous_col,
    omics_type = input$omics_type
  )

  list(
    results_raw = raw_df,
    results_std = results_std,
    model_object = NULL,
    analysis_info = list(
      omics_type = input$omics_type,
      method = "lm",
      analysis_type = "continuous_linear",
      comparison = continuous_col,
      covariates = covariates
    )
  )
}
