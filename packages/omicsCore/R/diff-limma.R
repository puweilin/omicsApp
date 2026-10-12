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
#' @param contrasts Parsed contrast specs (from `run_diff(contrasts = )`); when
#'   given, every group they name is fitted and each contrast read off the fit.
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
  paired_col = NULL,
  contrasts = NULL
) {
  validate_omics_input(input)
  ensure_limma()

  expr_mat <- input$expr_mat
  meta_df <- input$meta_df
  feature_df <- input$feature_df

  # `case_group` may name several groups, or `contrasts` any comparisons
  # between groups. Every group involved is fitted in one model, so every
  # contrast shares one residual variance and one eBayes prior -- which
  # is the point of fitting them together rather than as separate
  # two-group runs.
  tm <- contrast_target_meta(meta_df, group_col, control_group, case_group,
                             contrasts = contrasts, paired_col = paired_col)
  specs <- tm$specs
  group_levels <- tm$group_levels
  target_meta <- tm$target_meta

  keep_samples <- rownames(target_meta)
  expr_sub <- expr_mat[, keep_samples, drop = FALSE]

  formula_str <- paste0("~ 0 + `", group_col, "`")
  if (!is.null(covariates)) {
    missing_cov <- setdiff(covariates, colnames(target_meta))
    if (length(missing_cov) > 0L) {
      stop("Missing covariates: ", paste(missing_cov, collapse = ", "))
    }
    formula_str <- paste0(formula_str, " + ",
                          paste0("`", covariates, "`", collapse = " + "))
  }
  design <- stats::model.matrix(stats::as.formula(formula_str), data = target_meta)
  # The group columns come first, one per level in `group_levels` order.
  # They are renamed to placeholders rather than to the level names:
  # makeContrasts() only accepts syntactic names, and group labels such
  # as "24h", "Drug A" or "KO-1" are not -- that used to stop the run.
  n_lv <- length(group_levels)
  grp_names <- paste0(".grp", seq_len(n_lv))
  colnames(design)[seq_len(n_lv)] <- grp_names
  colnames(design) <- make.names(colnames(design), unique = TRUE)

  if (!is.null(paired_col)) {
    if (any(is.na(target_meta[[paired_col]]))) {
      stop("`paired_col` contains missing values after group filtering: ", paired_col)
    }
    bf <- limma_blocked_fit(expr_sub, design, target_meta[[paired_col]])
    fit <- bf$fit
    design <- bf$design
  } else {
    bf <- NULL
    fit <- limma::lmFit(expr_sub, design)
  }

  # The contrast matrix is built from the weights directly, one row per
  # design column (zero for covariates), rather than parsed from strings.
  cw <- contrast_matrix_from_specs(specs, group_levels)
  contrast_matrix <- matrix(0, nrow = ncol(design), ncol = ncol(cw),
                            dimnames = list(colnames(design), colnames(cw)))
  contrast_matrix[grp_names, ] <- cw
  fit2 <- limma::eBayes(limma::contrasts.fit(fit, contrast_matrix))

  st <- stack_contrast_results(
    specs,
    raw_for = function(i) {
      raw_df <- limma::topTable(fit2, coef = i, number = Inf, sort.by = "none")
      tibble::rownames_to_column(raw_df, "feature_id")
    },
    standardize = function(raw_df, comparison) {
      standardize_limma_group_results(
        raw_df = raw_df,
        feature_df = feature_df,
        comparison = comparison,
        omics_type = input$omics_type
      )
    }
  )
  raw_df <- st$raw
  results_std <- st$std
  comparisons <- st$comparisons

  list(
    results_raw = raw_df,
    results_std = results_std,
    model_object = fit2,
    analysis_info = list(
      omics_type = input$omics_type,
      method = "limma",
      analysis_type = "group",
      comparison = comparisons,
      covariates = covariates,
      paired_col = paired_col,
      pairing = bf$how,
      warnings = bf$note
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
#' @param model Either `"linear"` or `"spline"`. (Called `model` rather
#'   than `method` because [run_diff()] has a `method` of its own, and the
#'   clash made the spline unreachable: `run_diff_continuous(...,
#'   method = "spline")` failed, and `df = 3` alone ran a linear fit.)
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
  model = c("linear", "spline"),
  df = 3,
  covariates = NULL,
  paired_col = NULL
) {
  validate_omics_input(input)
  ensure_limma()
  method <- match.arg(model)

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
    bf <- limma_blocked_fit(expr_mat, design, meta_df[[paired_col]])
    fit <- limma::eBayes(bf$fit)
  } else {
    bf <- NULL
    fit <- limma::eBayes(limma::lmFit(expr_mat, design))
  }

  if (method == "spline") {
    coef_idx <- 2:(1L + df)
    raw_df <- limma::topTable(fit, coef = coef_idx, number = Inf, sort.by = "none")
  } else {
    raw_df <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "none")
  }
  raw_df <- tibble::rownames_to_column(raw_df, "feature_id")

  # The adjustment terms under placeholder names, with the pairing column
  # as a factor: numeric subject ids 1..6 used to enter lm() as a slope,
  # which gave adjusted R^2 = 0.245 where the blocked fit gives -0.242,
  # and a covariate called "Body mass" broke the formula outright.
  adj_cols <- unique(c(if (!is.null(paired_col)) paired_col, covariates))
  adj_df <- meta_df[, adj_cols, drop = FALSE]
  if (!is.null(paired_col)) adj_df[[paired_col]] <- factor(adj_df[[paired_col]])
  adjustment_terms <- if (length(adj_cols)) paste0(".adj", seq_along(adj_cols)) else character(0)
  names(adj_df) <- adjustment_terms

  # Model fit and rank correlation for every feature at once. The
  # apply() loops this replaces had no guard: one all-missing protein
  # stopped the whole run with "0 (non-NA) cases".
  ex <- limma_continuous_extras(expr_mat, cont_vals, adj_df, method, df)
  adj_r2 <- ex$adj_r2
  rho <- ex$rho

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
      paired_col = paired_col,
      pairing = bf$how,
      warnings = bf$note
    )
  )
}

#' Limma multi-group ANOVA-style test
#'
#' Moderated F-test on all group coefficients together, standardized with
#' the other global tests (the F stands in for the effect). Internal —
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
  target_meta <- target_meta[!is.na(target_meta[[group_col]]), , drop = FALSE]
  target_meta[[group_col]] <- factor(target_meta[[group_col]])

  keep_samples <- rownames(target_meta)
  expr_sub <- expr_mat[, keep_samples, drop = FALSE]

  if (nlevels(target_meta[[group_col]]) < 2L) {
    stop("ANOVA needs at least two groups in `", group_col, "`.", call. = FALSE)
  }

  # Intercept + treatment coding: the group coefficients are then the
  # differences from the first level, and the F-test on all of them
  # together asks "do the group means differ?". The cell-means design
  # (`~ 0 + group`) this used to fit made the same F-test ask "are all
  # group means zero?" -- true of almost nothing on a log-intensity
  # scale, so nearly every feature came out significant.
  formula_str <- paste0("~ `", group_col, "`")
  if (!is.null(covariates)) {
    missing_cov <- setdiff(covariates, colnames(target_meta))
    if (length(missing_cov) > 0L) {
      stop("Missing covariates: ", paste(missing_cov, collapse = ", "))
    }
    formula_str <- paste0(formula_str, " + ",
                          paste0("`", covariates, "`", collapse = " + "))
  }
  design <- stats::model.matrix(stats::as.formula(formula_str), data = target_meta)

  if (!is.null(paired_col)) {
    if (any(is.na(target_meta[[paired_col]]))) {
      stop("`paired_col` contains missing values after group filtering: ", paired_col)
    }
    bf <- limma_blocked_fit(expr_sub, design, target_meta[[paired_col]])
    fit <- limma::eBayes(bf$fit)
  } else {
    bf <- NULL
    fit <- limma::eBayes(limma::lmFit(expr_sub, design))
  }

  coef_idx <- 1L + seq_len(nlevels(target_meta[[group_col]]) - 1L)
  raw_df <- limma::topTable(fit, coef = coef_idx, number = Inf, sort.by = "none")
  # Two groups make one coefficient, and topTable() then reports t, not
  # F; the run stopped on "Column `F` not found". F on one numerator
  # degree of freedom is t squared.
  if (!"F" %in% names(raw_df) && "t" %in% names(raw_df)) raw_df$F <- raw_df$t^2
  raw_df <- tibble::rownames_to_column(raw_df, "feature_id")

  results_std <- standardize_global_test_results(
    raw_df, feature_df, group_col, input$omics_type, method = "limma",
    stat = "F", stat_type = "F", p = "P.Value", padj = "adj.P.Val", base_mean = "AveExpr")

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
      paired_col = paired_col,
      pairing = bf$how,
      warnings = bf$note
    )
  )
}

# Pairs (or repeated measures) as limma models them.
#
# As a fixed block -- the pair as a factor in the design, limma's own
# recommendation for paired samples -- whenever that design can be fitted:
# every pair has at least two samples, and no covariate is constant
# within pairs. duplicateCorrelation(), which this always used, assumes
# one within-pair correlation shared by every feature; on data where it
# varies (most data), the paired test was anticonservative: 1.8% of null
# features under p < 0.01 instead of 1%, and an observed FDR of 0.46 at
# n = 4 pairs against a nominal 0.05. The random-effect fit remains for
# the designs only it can estimate (a subject-level covariate such as
# sex, or pairs with one sample), and says so.
limma_blocked_fit <- function(expr, design, block) {
  block <- factor(block)
  if (nlevels(block) >= 2L && all(table(block) >= 2L)) {
    bm <- stats::model.matrix(~ block)[, -1L, drop = FALSE]
    colnames(bm) <- paste0(".blk", seq_len(ncol(bm)))
    full <- cbind(design, bm)
    if (qr(full)$rank == ncol(full)) {
      return(list(fit = limma::lmFit(expr, full), design = full,
                  how = "fixed_block", note = NULL))
    }
  }
  corfit <- limma::duplicateCorrelation(expr, design, block = block)
  why <- if (any(table(block) < 2L)) "some pairs have a single sample"
         else "a covariate does not vary within pairs"
  list(
    fit = limma::lmFit(expr, design, block = block, correlation = corfit$consensus),
    design = design, how = "random_block",
    note = sprintf(paste(
      "Pairs entered limma as a random effect (duplicateCorrelation, consensus",
      "correlation %.2f) because %s; this assumes one within-pair correlation",
      "for every feature, and p-values can be optimistic."),
      corfit$consensus, why))
}
