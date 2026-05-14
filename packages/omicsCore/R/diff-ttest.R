#' Per-feature two-group t-test
#'
#' Welch (default) or paired t-test per feature, comparing `case_group`
#' against `control_group`. RNA-seq raw counts are automatically log2(x+1)
#' transformed before testing. Internal — call via [run_diff()].
#'
#' @param input A validated `omics_input`.
#' @param group_col Group column in sample metadata.
#' @param control_group Control-group label.
#' @param case_group Case-group label.
#' @param var_equal If `TRUE`, use equal-variance t-test. Default `FALSE`
#'   (Welch).
#' @param paired_col Optional pairing column for paired t-test.
#'
#' @return List with `results_raw`, `results_std`, `model_object` (`NULL`),
#'   and `analysis_info`.
#' @keywords internal
run_ttest_group <- function(
  input,
  group_col,
  control_group,
  case_group,
  var_equal = FALSE,
  paired_col = NULL
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

  ctrl_samples <- rownames(target_meta[target_meta[[group_col]] == control_group, , drop = FALSE])
  case_samples <- rownames(target_meta[target_meta[[group_col]] == case_group, , drop = FALSE])
  if (length(ctrl_samples) < 2L || length(case_samples) < 2L) {
    stop("Each group must have at least 2 samples for t-test.")
  }

  is_paired <- !is.null(paired_col)
  if (is_paired) {
    check_paired_col(meta_df, paired_col, object_name = "meta_df")
    validate_two_group_pairing(
      target_meta,
      group_col = group_col,
      paired_col = paired_col,
      control_group = control_group,
      case_group = case_group,
      object_name = "target_meta"
    )
    ctrl_order <- target_meta[ctrl_samples, paired_col]
    case_order <- target_meta[case_samples, paired_col]
    case_samples <- case_samples[match(ctrl_order, case_order)]
  }

  expr_sub <- expr_mat[, c(ctrl_samples, case_samples), drop = FALSE]
  feature_ids <- rownames(expr_sub)
  n_features <- length(feature_ids)

  mean_ctrl <- numeric(n_features)
  mean_case <- numeric(n_features)
  mean_diff <- numeric(n_features)
  t_stat <- numeric(n_features)
  p_value <- numeric(n_features)

  for (i in seq_len(n_features)) {
    x_ctrl <- as.numeric(expr_sub[i, ctrl_samples])
    x_case <- as.numeric(expr_sub[i, case_samples])

    tt <- tryCatch(
      stats::t.test(x_case, x_ctrl, paired = is_paired, var.equal = var_equal),
      error = function(e) NULL
    )
    mean_ctrl[i] <- mean(x_ctrl, na.rm = TRUE)
    mean_case[i] <- mean(x_case, na.rm = TRUE)
    mean_diff[i] <- mean_case[i] - mean_ctrl[i]
    if (!is.null(tt)) {
      t_stat[i] <- tt$statistic
      p_value[i] <- tt$p.value
    } else {
      t_stat[i] <- NA_real_
      p_value[i] <- NA_real_
    }
  }

  raw_df <- data.frame(
    feature_id = feature_ids,
    mean_ctrl = mean_ctrl,
    mean_case = mean_case,
    mean_diff = mean_diff,
    t_stat = t_stat,
    p_value = p_value,
    adj_p_value = stats::p.adjust(p_value, method = "BH"),
    stringsAsFactors = FALSE
  )

  comparison <- paste0(case_group, "_vs_", control_group)
  results_std <- standardize_ttest_group_results(
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
      method = "ttest",
      analysis_type = "group",
      comparison = comparison,
      var_equal = var_equal,
      paired_col = paired_col
    )
  )
}
