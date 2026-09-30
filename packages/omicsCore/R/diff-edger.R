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

  group_levels <- c(control_group, case_group)
  meta_df[[group_col]] <- as.character(meta_df[[group_col]])
  target_meta <- meta_df[!is.na(meta_df[[group_col]]) &
                           meta_df[[group_col]] %in% group_levels, , drop = FALSE]
  target_meta[[group_col]] <- factor(target_meta[[group_col]], levels = group_levels)
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

  design_formula <- stats::as.formula(paste("~", paste(design_terms, collapse = " + ")))
  design_mat <- stats::model.matrix(design_formula, data = target_meta)

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

  # model.matrix() names a factor's columns `<column><level>` verbatim,
  # so the coefficient is found by exact name. It used to be found by
  # substring, which picked "B" out of "groupAB" when both were present.
  comparisons <- paste0(case_group, "_vs_", control_group)
  per <- lapply(seq_along(case_group), function(i) {
    group_coef <- paste0(group_col, case_group[[i]])
    if (!group_coef %in% colnames(design_mat)) {
      stop("Could not locate group coefficient in design matrix for: ",
           case_group[[i]])
    }
    qlf <- edgeR::glmQLFTest(fit, coef = group_coef)
    tt <- edgeR::topTags(qlf, n = Inf, sort.by = "none")
    raw_df <- as.data.frame(tt$table) |>
      tibble::rownames_to_column("feature_id")
    std <- standardize_edger_group_results(
      raw_df = raw_df,
      feature_df = feature_df,
      comparison = comparisons[[i]],
      omics_type = input$omics_type
    )
    raw_df$comparison <- comparisons[[i]]
    list(raw = raw_df, std = std)
  })
  raw_df <- do.call(rbind, lapply(per, `[[`, "raw"))
  if (length(case_group) == 1L) raw_df$comparison <- NULL
  results_std <- do.call(rbind, lapply(per, `[[`, "std"))
  rownames(raw_df) <- NULL
  rownames(results_std) <- NULL
  comparison <- comparisons

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
