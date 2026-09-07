# edgeR QL F-test backend, Suggests-gated. Uses an optional tximport length
# matrix (from input$misc$tximport$length) to scale offsets when available;
# otherwise falls back to standard library-size normalization.

ensure_edger <- function() {
  if (!is_installed("edgeR")) {
    stop(
      "Package 'edgeR' is required for the edgeR differential backend. ",
      "Install with: omicsCore::install_optional('rnaseq').",
      call. = FALSE
    )
  }
}

#' edgeR QL F-test for a two-group comparison
#'
#' @param input A validated RNA-seq `omics_input` built from raw counts.
#' @param group_col Group column in sample metadata.
#' @param control_group Control-group label.
#' @param case_group Case-group label.
#' @param covariates Optional covariate column names.
#' @param paired_col Optional pairing column.
#'
#' @return List with `results_raw`, `results_std`, `model_object`
#'   (`DGEGLM`), and `analysis_info`.
#' @keywords internal
run_edger_group <- function(
  input,
  group_col,
  control_group,
  case_group,
  covariates = NULL,
  paired_col = NULL
) {
  validate_omics_input(input)
  if (input$omics_type != "rnaseq") {
    stop("`run_edger_group()` requires `omics_type = 'rnaseq'`.")
  }
  if (!identical(input$assay_type, "raw_count")) {
    stop("`run_edger_group()` requires `assay_type = 'raw_count'`.")
  }
  ensure_edger()

  count_mat <- as.matrix(input$expr_mat)
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
  count_sub <- count_mat[, keep_samples, drop = FALSE]

  design_terms <- c(group_col)
  if (!is.null(paired_col)) {
    if (any(is.na(target_meta[[paired_col]]))) {
      stop("`paired_col` contains missing values after group filtering: ", paired_col)
    }
    target_meta[[paired_col]] <- factor(target_meta[[paired_col]])
    design_terms <- c(paired_col, design_terms)
  }
  if (!is.null(covariates)) {
    missing_cov <- setdiff(covariates, colnames(target_meta))
    if (length(missing_cov) > 0L) {
      stop("Missing covariates: ", paste(missing_cov, collapse = ", "))
    }
    design_terms <- c(design_terms, covariates)
  }

  design_mat <- build_design(target_meta, design_terms, primary = group_col)

  y <- edgeR::DGEList(counts = count_sub)
  txi_info <- get_tximport_info(input)
  if (!is.null(txi_info$length)) {
    length_sub <- txi_info$length[rownames(count_sub), colnames(count_sub), drop = FALSE]
    lib_sizes <- colSums(count_sub)
    log_length <- log(length_sub + 1)
    log_lib <- matrix(
      log(lib_sizes),
      nrow = nrow(count_sub),
      ncol = ncol(count_sub),
      byrow = TRUE
    )
    offsets <- log_length + log_lib
    y <- edgeR::scaleOffset(y, offset = offsets)
  } else {
    y <- edgeR::calcNormFactors(y)
  }

  y <- edgeR::estimateDisp(y, design = design_mat)
  fit <- edgeR::glmQLFit(y, design = design_mat)

  # The coefficient is addressed by its generated name. Searching for it with
  # `grep(group_col, ...)` treated the column name as a regex and the group
  # label as a substring, so a covariate whose name contained the group's
  # could be picked instead.
  group_coef <- design_coef_name(design_mat, attr(design_mat, "primary_id"), case_group)
  if (!group_coef %in% colnames(design_mat)) {
    stop("Could not locate group coefficient in design matrix for: ", case_group)
  }

  qlf <- edgeR::glmQLFTest(fit, coef = group_coef)
  tt <- edgeR::topTags(qlf, n = Inf, sort.by = "none")
  raw_df <- as.data.frame(tt$table) |>
    tibble::rownames_to_column("feature_id")

  comparison <- paste0(case_group, "_vs_", control_group)
  results_std <- standardize_edger_group_results(
    raw_df = raw_df,
    feature_df = feature_df,
    comparison = comparison,
    omics_type = input$omics_type
  )

  list(
    results_raw = raw_df,
    results_std = results_std,
    model_object = fit,
    analysis_info = list(
      omics_type = input$omics_type,
      method = "edger",
      analysis_type = "group",
      comparison = comparison,
      covariates = covariates,
      paired_col = paired_col
    )
  )
}
