# All limma backends are Suggests-gated. When the user calls run_diff() with
# method = "limma" but limma is not installed, ensure_limma() emits a hint
# pointing at install_optional("proteomics").

ensure_limma <- function() {
  if (!is_installed("limma")) {
    stop(
      "Package 'limma' is required for the limma differential backend. ",
      "Install with: omicsCore::install_optional('proteomics').",
      call. = FALSE
    )
  }
}

#' Limma two-group differential test
#'
#' Fits `~ 0 + group_col [+ covariates]` and contrasts `case_group -
#' control_group`. Optionally accommodates a `paired_col` via
#' `limma::duplicateCorrelation`. Internal — call via [run_diff()].
#'
#' @param input A validated proteomics-style `omics_input`.
#' @param group_col Group column in sample metadata.
#' @param control_group Control-group label.
#' @param case_group Case-group label.
#' @param covariates Optional covariate column names.
#' @param paired_col Optional pairing/block column.
#'
#' @return List with `results_raw`, `results_std`, `model_object`, and
#'   `analysis_info`.
#' @keywords internal
run_limma_group <- function(
  input,
  group_col,
  control_group,
  case_group,
  covariates = NULL,
  paired_col = NULL
) {
  validate_omics_input(input)
  ensure_limma()

  expr_mat <- input$expr_mat
  meta_df <- input$meta_df
  feature_df <- input$feature_df

  if (!group_col %in% colnames(meta_df)) {
    stop("`group_col` not found in `meta_df`: ", group_col)
  }
  check_paired_col(meta_df, paired_col, object_name = "meta_df")

  meta_df[[group_col]] <- factor(meta_df[[group_col]])
  target_meta <- meta_df[meta_df[[group_col]] %in% c(control_group, case_group), , drop = FALSE]
  target_meta[[group_col]] <- factor(target_meta[[group_col]], levels = c(control_group, case_group))
  validate_two_group_pairing(
    target_meta,
    group_col = group_col,
    paired_col = paired_col,
    control_group = control_group,
    case_group = case_group,
    object_name = "target_meta"
  )

  keep_samples <- rownames(target_meta)
  expr_sub <- expr_mat[, keep_samples, drop = FALSE]

  formula_str <- paste0("~ 0 + ", group_col)
  if (!is.null(covariates)) {
    missing_cov <- setdiff(covariates, colnames(target_meta))
    if (length(missing_cov) > 0L) {
      stop("Missing covariates: ", paste(missing_cov, collapse = ", "))
    }
    formula_str <- paste0(formula_str, " + ", paste(covariates, collapse = " + "))
  }
  design <- stats::model.matrix(stats::as.formula(formula_str), data = target_meta)
  colnames(design) <- gsub(paste0("^", group_col), "", colnames(design))

  if (!is.null(paired_col)) {
    if (any(is.na(target_meta[[paired_col]]))) {
      stop("`paired_col` contains missing values after group filtering: ", paired_col)
    }
    block_var <- factor(target_meta[[paired_col]])
    corfit <- limma::duplicateCorrelation(expr_sub, design, block = block_var)
    fit <- limma::lmFit(expr_sub, design, block = block_var, correlation = corfit$consensus)
  } else {
    fit <- limma::lmFit(expr_sub, design)
  }

  contrast_str <- paste0(make.names(case_group), " - ", make.names(control_group))
  contrast_matrix <- limma::makeContrasts(contrasts = contrast_str, levels = design)
  fit2 <- limma::eBayes(limma::contrasts.fit(fit, contrast_matrix))

  raw_df <- limma::topTable(fit2, number = Inf, sort.by = "none")
  raw_df <- tibble::rownames_to_column(raw_df, "feature_id")

  comparison <- paste0(case_group, "_vs_", control_group)
  results_std <- standardize_limma_group_results(
    raw_df = raw_df,
    feature_df = feature_df,
    comparison = comparison,
    omics_type = input$omics_type
  )

  list(
    results_raw = raw_df,
    results_std = results_std,
    model_object = fit2,
    analysis_info = list(
      omics_type = input$omics_type,
      method = "limma",
      analysis_type = "group",
      comparison = comparison,
      covariates = covariates,
      paired_col = paired_col
    )
  )
}

#' Limma continuous-variable differential test
#'
#' Either linear or spline (natural cubic via `splines::ns`) modeling of
#' `continuous_col`, returning per-feature limma statistics plus a (partial)
#' Spearman rho and adjusted R^2. Internal — call via [run_diff()].
#'
#' @param input A validated `omics_input`.
#' @param continuous_col Continuous metadata column.
#' @param method Either `"linear"` or `"spline"`.
#' @param df Degrees of freedom for spline fits.
#' @param covariates Optional covariate column names.
#' @param paired_col Optional pairing/block column.
#'
#' @return List with `results_raw`, `results_std`, `model_object`, and
#'   `analysis_info`.
#' @keywords internal
run_limma_continuous <- function(
  input,
  continuous_col,
  method = c("linear", "spline"),
  df = 3,
  covariates = NULL,
  paired_col = NULL
) {
  validate_omics_input(input)
  ensure_limma()
  method <- match.arg(method)

  expr_mat <- input$expr_mat
  meta_df <- input$meta_df
  feature_df <- input$feature_df

  if (!continuous_col %in% colnames(meta_df)) {
    stop("`continuous_col` not found in `meta_df`: ", continuous_col)
  }
  check_paired_col(meta_df, paired_col, object_name = "meta_df")
  validate_continuous_pairing(meta_df, paired_col, object_name = "meta_df")

  cont_vals <- coerce_continuous_col(meta_df[[continuous_col]], continuous_col)
  meta_df[[continuous_col]] <- cont_vals

  if (method == "spline") {
    design_df <- as.data.frame(splines::ns(cont_vals, df = df))
    colnames(design_df) <- paste0("ns", seq_len(ncol(design_df)))
  } else {
    design_df <- data.frame(x = cont_vals)
    colnames(design_df) <- continuous_col
  }

  if (!is.null(covariates)) {
    missing_cov <- setdiff(covariates, colnames(meta_df))
    if (length(missing_cov) > 0L) {
      stop("Missing covariates: ", paste(missing_cov, collapse = ", "))
    }
    for (cov in covariates) {
      design_df[[cov]] <- meta_df[[cov]]
    }
  }

  design <- stats::model.matrix(~ ., data = design_df)
  if (!is.null(paired_col)) {
    if (any(is.na(meta_df[[paired_col]]))) {
      stop("`paired_col` contains missing values: ", paired_col)
    }
    block_var <- factor(meta_df[[paired_col]])
    corfit <- limma::duplicateCorrelation(expr_mat, design, block = block_var)
    fit <- limma::eBayes(limma::lmFit(
      expr_mat, design, block = block_var, correlation = corfit$consensus
    ))
  } else {
    fit <- limma::eBayes(limma::lmFit(expr_mat, design))
  }

  if (method == "spline") {
    coef_idx <- 2:(1L + df)
    raw_df <- limma::topTable(fit, coef = coef_idx, number = Inf, sort.by = "none")
  } else {
    raw_df <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "none")
  }
  raw_df <- tibble::rownames_to_column(raw_df, "feature_id")

  adjustment_terms <- unique(c(if (!is.null(paired_col)) paired_col, covariates))

  adj_r2 <- apply(expr_mat, 1L, function(y) {
    model_df <- data.frame(y = y, cont = cont_vals, meta_df[, adjustment_terms, drop = FALSE])
    if (method == "spline") {
      f <- if (length(adjustment_terms) > 0L) {
        stats::as.formula(paste(
          "y ~ splines::ns(cont, df =", df, ") +",
          paste(adjustment_terms, collapse = " + ")
        ))
      } else {
        stats::as.formula(paste("y ~ splines::ns(cont, df =", df, ")"))
      }
    } else {
      f <- if (length(adjustment_terms) > 0L) {
        stats::as.formula(paste("y ~ cont +", paste(adjustment_terms, collapse = " + ")))
      } else {
        stats::as.formula("y ~ cont")
      }
    }
    summary(stats::lm(f, data = model_df))$adj.r.squared
  })

  rho <- apply(expr_mat, 1L, function(y) {
    if (length(adjustment_terms) == 0L) {
      return(suppressWarnings(
        stats::cor.test(y, cont_vals, method = "spearman", exact = FALSE)$estimate
      ))
    }
    model_df <- data.frame(y = y, cont = cont_vals, meta_df[, adjustment_terms, drop = FALSE])
    y_resid <- stats::residuals(stats::lm(
      stats::as.formula(paste("y ~", paste(adjustment_terms, collapse = " + "))),
      data = model_df
    ))
    cont_resid <- stats::residuals(stats::lm(
      stats::as.formula(paste("cont ~", paste(adjustment_terms, collapse = " + "))),
      data = model_df
    ))
    suppressWarnings(stats::cor.test(y_resid, cont_resid, method = "spearman", exact = FALSE)$estimate)
  })

  raw_df$adj_r_squared <- unname(adj_r2[raw_df$feature_id])
  raw_df$spearman_rho <- unname(rho[raw_df$feature_id])

  results_std <- standardize_limma_continuous_results(
    raw_df = raw_df,
    feature_df = feature_df,
    comparison = continuous_col,
    omics_type = input$omics_type
  )

  list(
    results_raw = raw_df,
    results_std = results_std,
    model_object = fit,
    analysis_info = list(
      omics_type = input$omics_type,
      method = "limma",
      analysis_type = if (method == "spline") "continuous_spline" else "continuous_linear",
      comparison = continuous_col,
      covariates = covariates,
      paired_col = paired_col
    )
  )
}

#' Limma multi-group ANOVA-style test
#'
#' Builds the standardized diff schema in place (no shared standardizer) so
#' that F-statistic and AveExpr fields are passed through cleanly. Internal —
#' call via [run_diff()] with `analysis_type = "anova"`.
#'
#' @param input A validated `omics_input`.
#' @param group_col Grouping column in sample metadata.
#' @param covariates Optional covariate column names.
#' @param selected_groups Optional subset of groups to retain.
#' @param paired_col Optional pairing/block column.
#'
#' @return List with `results_raw`, `results_std`, `model_object`, and
#'   `analysis_info`.
#' @keywords internal
run_limma_anova <- function(
  input,
  group_col,
  covariates = NULL,
  selected_groups = NULL,
  paired_col = NULL
) {
  validate_omics_input(input)
  ensure_limma()

  expr_mat <- input$expr_mat
  meta_df <- input$meta_df
  feature_df <- input$feature_df

  if (!group_col %in% colnames(meta_df)) {
    stop("`group_col` not found in `meta_df`: ", group_col)
  }
  check_paired_col(meta_df, paired_col, object_name = "meta_df")

  target_meta <- meta_df
  if (!is.null(selected_groups)) {
    target_meta <- target_meta[target_meta[[group_col]] %in% selected_groups, , drop = FALSE]
  }
  target_meta[[group_col]] <- factor(target_meta[[group_col]])

  keep_samples <- rownames(target_meta)
  expr_sub <- expr_mat[, keep_samples, drop = FALSE]

  formula_str <- paste0("~ 0 + ", group_col)
  if (!is.null(covariates)) {
    missing_cov <- setdiff(covariates, colnames(target_meta))
    if (length(missing_cov) > 0L) {
      stop("Missing covariates: ", paste(missing_cov, collapse = ", "))
    }
    formula_str <- paste0(formula_str, " + ", paste(covariates, collapse = " + "))
  }
  design <- stats::model.matrix(stats::as.formula(formula_str), data = target_meta)

  if (!is.null(paired_col)) {
    if (any(is.na(target_meta[[paired_col]]))) {
      stop("`paired_col` contains missing values after group filtering: ", paired_col)
    }
    block_var <- factor(target_meta[[paired_col]])
    corfit <- limma::duplicateCorrelation(expr_sub, design, block = block_var)
    fit <- limma::eBayes(limma::lmFit(
      expr_sub, design, block = block_var, correlation = corfit$consensus
    ))
  } else {
    fit <- limma::eBayes(limma::lmFit(expr_sub, design))
  }

  coef_idx <- seq_len(nlevels(target_meta[[group_col]]))
  raw_df <- limma::topTable(fit, coef = coef_idx, number = Inf, sort.by = "none")
  raw_df <- tibble::rownames_to_column(raw_df, "feature_id")

  feature_df <- prep_feature_df_for_standardize(feature_df)
  results_std <- dplyr::left_join(raw_df, feature_df, by = "feature_id") |>
    dplyr::transmute(
      feature_id = .data$feature_id,
      feature_symbol = .data$feature_symbol,
      feature_type = .data$feature_type,
      omics_type = input$omics_type,
      method = "limma",
      analysis_type = "anova",
      comparison = group_col,
      effect = .data$F,
      effect_type = "F_statistic",
      statistic = .data$F,
      statistic_type = "F",
      p_value = .data$P.Value,
      adj_p_value = .data$adj.P.Val,
      direction = "ns",
      base_mean = if ("AveExpr" %in% colnames(raw_df)) .data$AveExpr else NA_real_,
      model_fit = NA_real_,
      is_significant = FALSE
    )
  check_diff_result_schema(results_std)

  list(
    results_raw = raw_df,
    results_std = results_std,
    model_object = fit,
    analysis_info = list(
      omics_type = input$omics_type,
      method = "limma",
      analysis_type = "anova",
      comparison = group_col,
      covariates = covariates,
      selected_groups = selected_groups,
      paired_col = paired_col
    )
  )
}
