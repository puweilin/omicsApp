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
#' @param incomplete_pairs What a paired test does with a pair that lacks a
#'   sample in one of the two groups (or has two in one): `"error"` (the
#'   default, for a single comparison) refuses the design; `"drop"` leaves
#'   such pairs out and tests the complete ones, which is how [run_diff()]
#'   runs each comparison when several treatment groups share a control and
#'   not every subject received every treatment. At least two complete
#'   pairs are needed either way.
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
  paired_col = NULL,
  incomplete_pairs = c("error", "drop")
) {
  incomplete_pairs <- match.arg(incomplete_pairs)
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
  is_paired <- !is.null(paired_col)
  # Pairing per comparison counts complete pairs instead, and says so in
  # its own words.
  if (!(is_paired && identical(incomplete_pairs, "drop")) &&
      (length(ctrl_samples) < 2L || length(case_samples) < 2L)) {
    stop("Each group must have at least 2 samples for t-test.")
  }

  pairs <- NULL
  if (is_paired) {
    check_paired_col(meta_df, paired_col, object_name = "meta_df")
    if (identical(incomplete_pairs, "drop")) {
      pairs <- complete_pairs(target_meta, group_col, paired_col,
                              control_group, case_group)
      if (length(pairs$used) < 2L) {
        stop(sprintf(paste(
          "The paired t-test of %s against %s needs at least 2 complete pairs",
          "(a `%s` with one sample in each group); %s."),
          case_group, control_group, paired_col,
          if (length(pairs$used)) "there is 1" else "there are none"), call. = FALSE)
      }
      keep <- !is.na(target_meta[[paired_col]]) &
        as.character(target_meta[[paired_col]]) %in% pairs$used
      target_meta <- target_meta[keep, , drop = FALSE]
      ctrl_samples <- rownames(target_meta[target_meta[[group_col]] == control_group, , drop = FALSE])
      case_samples <- rownames(target_meta[target_meta[[group_col]] == case_group, , drop = FALSE])
    } else {
      validate_two_group_pairing(
        target_meta,
        group_col = group_col,
        paired_col = paired_col,
        control_group = control_group,
        case_group = case_group,
        object_name = "target_meta"
      )
    }
    ctrl_order <- target_meta[ctrl_samples, paired_col]
    case_order <- target_meta[case_samples, paired_col]
    case_samples <- case_samples[match(ctrl_order, case_order)]
  }

  expr_sub <- expr_mat[, c(ctrl_samples, case_samples), drop = FALSE]
  feature_ids <- rownames(expr_sub)
  n_features <- length(feature_ids)

  # Every feature at once (ttest_rows(), fast-rows.R) rather than one
  # t.test() per feature: the same arithmetic, 18x faster on a pairwise
  # run over 8,000 proteins.
  ctrl_m <- expr_sub[, ctrl_samples, drop = FALSE]
  case_m <- expr_sub[, case_samples, drop = FALSE]
  tt <- ttest_rows(case_m, ctrl_m, paired = is_paired, var_equal = var_equal)
  if (is_paired) {
    # The effect over the pairs the test used. Means over every sample
    # included the halves of incomplete pairs, and could point the other
    # way from t: -3.1 "down" beside t = +12.
    ok <- !is.na(ctrl_m) & !is.na(case_m)
    mean_ctrl <- unname(rowSums(ifelse(ok, ctrl_m, 0)) / rowSums(ok))
    mean_case <- unname(rowSums(ifelse(ok, case_m, 0)) / rowSums(ok))
  } else {
    mean_ctrl <- unname(rowMeans(ctrl_m, na.rm = TRUE))
    mean_case <- unname(rowMeans(case_m, na.rm = TRUE))
  }
  mean_diff <- mean_case - mean_ctrl
  t_stat <- tt$t
  p_value <- tt$p

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
      paired_col = paired_col,
      pairs_used = if (!is.null(pairs)) length(pairs$used),
      pairs_left_out = if (!is.null(pairs)) pairs$dropped
    )
  )
}

# The pairs of one comparison: the values of `paired_col` with exactly
# one sample in each of the two groups (`used`), and those present in
# either group without being such a pair (`dropped`) -- a subject who
# missed one of the two treatments, or was sampled twice under one.
complete_pairs <- function(meta, group_col, paired_col, control_group, case_group) {
  g <- as.character(meta[[group_col]])
  p <- as.character(meta[[paired_col]])
  ok <- !is.na(g) & !is.na(p) & g %in% c(control_group, case_group)
  g <- g[ok]
  p <- p[ok]
  ids <- sort(unique(p))
  n_ctrl <- vapply(ids, function(i) sum(p == i & g == control_group), integer(1))
  n_case <- vapply(ids, function(i) sum(p == i & g == case_group), integer(1))
  complete <- n_ctrl == 1L & n_case == 1L
  list(used = ids[complete], dropped = ids[!complete])
}
