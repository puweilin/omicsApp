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

  # Terms and levels renamed once; the loop below models them through `~ .`.
  # Besides keeping the uploader's strings away from the parser, this fixes a
  # standing breakage: `data.frame()` runs `make.names()` over its columns, so
  # a legitimate column called `年龄 (岁)` became `年龄..岁.` while the pasted
  # formula still asked for the original, and the fit failed for every feature.
  sf <- safe_frame(target_meta, c(group_col, covariates))
  group_id <- sf$ids[[1L]]
  formula_obj <- stats::as.formula("y ~ .")
  coef_name <- design_coef_name(sf, group_id, case_group)

  feature_ids <- rownames(expr_sub)
  n_features <- length(feature_ids)

  beta <- numeric(n_features)
  t_stat <- numeric(n_features)
  p_value <- numeric(n_features)
  adj_r_squared <- numeric(n_features)
  base_mean <- numeric(n_features)

  for (i in seq_len(n_features)) {
    model_df <- data.frame(y = as.numeric(expr_sub[i, ]), sf$df)
    fit <- tryCatch(stats::lm(formula_obj, data = model_df), error = function(e) NULL)
    base_mean[i] <- mean(model_df$y, na.rm = TRUE)

    if (!is.null(fit)) {
      s <- summary(fit)
      coefs <- s$coefficients
      if (coef_name %in% rownames(coefs)) {
        beta[i] <- coefs[coef_name, "Estimate"]
        t_stat[i] <- coefs[coef_name, "t value"]
        p_value[i] <- coefs[coef_name, "Pr(>|t|)"]
      } else {
        beta[i] <- NA_real_
        t_stat[i] <- NA_real_
        p_value[i] <- NA_real_
      }
      adj_r_squared[i] <- s$adj.r.squared
    } else {
      beta[i] <- NA_real_
      t_stat[i] <- NA_real_
      p_value[i] <- NA_real_
      adj_r_squared[i] <- NA_real_
    }
  }

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

  # As in run_lm_group(): generated names, `~ .`, no parsing of uploaded
  # strings. `cont` is the continuous term, `adj` the covariates.
  sf <- safe_frame(meta_df, covariates)
  formula_obj <- stats::as.formula("y ~ .")
  cont_coef <- "cont"
  adjustment_terms <- if (!is.null(covariates)) covariates else character(0)

  feature_ids <- rownames(expr_mat)
  n_features <- length(feature_ids)

  beta <- numeric(n_features)
  t_stat <- numeric(n_features)
  p_value <- numeric(n_features)
  adj_r_squared <- numeric(n_features)
  spearman_rho <- numeric(n_features)
  base_mean <- numeric(n_features)

  for (i in seq_len(n_features)) {
    y <- as.numeric(expr_mat[i, ])
    model_df <- data.frame(y = y, cont = cont_vals, sf$df)
    fit <- tryCatch(stats::lm(formula_obj, data = model_df), error = function(e) NULL)
    base_mean[i] <- mean(y, na.rm = TRUE)

    if (!is.null(fit)) {
      s <- summary(fit)
      coefs <- s$coefficients
      if (cont_coef %in% rownames(coefs)) {
        beta[i] <- coefs[cont_coef, "Estimate"]
        t_stat[i] <- coefs[cont_coef, "t value"]
        p_value[i] <- coefs[cont_coef, "Pr(>|t|)"]
      } else {
        beta[i] <- NA_real_
        t_stat[i] <- NA_real_
        p_value[i] <- NA_real_
      }
      adj_r_squared[i] <- s$adj.r.squared
    } else {
      beta[i] <- NA_real_
      t_stat[i] <- NA_real_
      p_value[i] <- NA_real_
      adj_r_squared[i] <- NA_real_
    }

    rho_val <- if (length(adjustment_terms) == 0L) {
      tryCatch(
        suppressWarnings(stats::cor.test(y, cont_vals, method = "spearman", exact = FALSE)$estimate),
        error = function(e) NA_real_
      )
    } else {
      tryCatch({
        y_resid <- stats::residuals(
          stats::lm(y ~ ., data = data.frame(y = y, sf$df)))
        cont_resid <- stats::residuals(
          stats::lm(cont ~ ., data = data.frame(cont = cont_vals, sf$df)))
        suppressWarnings(stats::cor.test(y_resid, cont_resid, method = "spearman", exact = FALSE)$estimate)
      }, error = function(e) NA_real_)
    }
    spearman_rho[i] <- unname(rho_val)
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
