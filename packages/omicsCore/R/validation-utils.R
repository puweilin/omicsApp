#' Check that required columns exist
#'
#' @param data Data frame to validate.
#' @param required_cols Required column names.
#' @param object_name Object label used in error messages.
#'
#' @return Invisibly returns `TRUE` on success.
#' @keywords internal
check_required_cols <- function(data, required_cols, object_name = "data") {
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in ", object_name, ": ",
      paste(missing_cols, collapse = ", ")
    )
  }
  invisible(TRUE)
}

#' Check that sample metadata matches an expression matrix
#'
#' @param expr_mat Expression matrix.
#' @param meta_df Sample metadata.
#'
#' @return Invisibly returns `TRUE` on success.
#' @keywords internal
check_meta_matches_expr <- function(expr_mat, meta_df) {
  if (is.null(colnames(expr_mat))) {
    stop("`expr_mat` must have sample column names.")
  }
  if (is.null(rownames(meta_df))) {
    stop("`meta_df` must have rownames.")
  }
  if (!all(colnames(expr_mat) %in% rownames(meta_df))) {
    stop("Not all samples in `expr_mat` are present in `rownames(meta_df)`.")
  }
  invisible(TRUE)
}

#' Check that feature metadata matches an expression matrix
#'
#' @param expr_mat Expression matrix.
#' @param feature_df Feature metadata.
#' @param feature_id_col Feature identifier column in `feature_df`.
#'
#' @return Invisibly returns `TRUE` on success.
#' @keywords internal
check_feature_matches_expr <- function(expr_mat, feature_df, feature_id_col = "feature_id") {
  if (is.null(rownames(expr_mat))) {
    stop("`expr_mat` must have feature row names.")
  }
  if (!feature_id_col %in% colnames(feature_df)) {
    stop("`feature_df` must contain `", feature_id_col, "`.")
  }
  if (!all(rownames(expr_mat) %in% feature_df[[feature_id_col]])) {
    stop("Not all features in `expr_mat` are present in `feature_df`.")
  }
  invisible(TRUE)
}

#' Check paired column validity in sample metadata
#'
#' @param meta_df Sample metadata.
#' @param paired_col Pairing column name.
#' @param object_name Object label used in error messages.
#'
#' @return Invisibly returns `TRUE` on success.
#' @keywords internal
check_paired_col <- function(meta_df, paired_col, object_name = "meta_df") {
  if (is.null(paired_col)) {
    return(invisible(TRUE))
  }

  if (!paired_col %in% colnames(meta_df)) {
    stop("`paired_col` not found in `", object_name, "`: ", paired_col)
  }

  pair_vals <- meta_df[[paired_col]]
  if (all(is.na(pair_vals))) {
    stop("`paired_col` contains only missing values: ", paired_col)
  }

  invisible(TRUE)
}

#' Validate two-group paired design completeness
#'
#' Ensures every pair has exactly one sample in each of the two target groups.
#'
#' @param meta_df Sample metadata already filtered to target samples.
#' @param group_col Group column name.
#' @param paired_col Pairing column name.
#' @param control_group Control-group label.
#' @param case_group Case-group label.
#' @param object_name Object label used in error messages.
#'
#' @return Invisibly returns `TRUE` on success.
#' @keywords internal
validate_two_group_pairing <- function(
  meta_df,
  group_col,
  paired_col,
  control_group,
  case_group,
  object_name = "meta_df"
) {
  if (is.null(paired_col)) {
    return(invisible(TRUE))
  }

  check_required_cols(meta_df, c(group_col, paired_col), object_name = object_name)

  pair_df <- data.frame(
    pair_id = meta_df[[paired_col]],
    group_id = meta_df[[group_col]],
    stringsAsFactors = FALSE
  )
  pair_df <- pair_df[!is.na(pair_df$pair_id) & !is.na(pair_df$group_id), , drop = FALSE]

  pair_levels <- split(pair_df$group_id, pair_df$pair_id)
  bad_pairs <- names(Filter(function(x) {
    length(x) != 2L ||
      !setequal(as.character(x), c(control_group, case_group))
  }, pair_levels))

  if (length(bad_pairs) > 0) {
    stop(
      "Invalid paired design in `", object_name, "` for `paired_col = '", paired_col, "'`. ",
      "Each pair must contain exactly one `", control_group, "` and one `", case_group, "` sample. ",
      "Problematic pairs: ", paste(bad_pairs, collapse = ", ")
    )
  }

  invisible(TRUE)
}

#' Validate continuous paired design minimum replication
#'
#' Ensures every pair has at least two observations.
#'
#' @param meta_df Sample metadata.
#' @param paired_col Pairing column name.
#' @param object_name Object label used in error messages.
#'
#' @return Invisibly returns `TRUE` on success.
#' @keywords internal
validate_continuous_pairing <- function(
  meta_df,
  paired_col,
  object_name = "meta_df"
) {
  if (is.null(paired_col)) {
    return(invisible(TRUE))
  }

  check_required_cols(meta_df, paired_col, object_name = object_name)
  pair_counts <- table(meta_df[[paired_col]], useNA = "no")
  bad_pairs <- names(pair_counts)[pair_counts < 2L]
  if (length(bad_pairs) > 0) {
    stop(
      "Invalid paired design in `", object_name, "` for `paired_col = '", paired_col, "'`. ",
      "Each pair must contain at least 2 samples. Problematic pairs: ",
      paste(bad_pairs, collapse = ", ")
    )
  }

  invisible(TRUE)
}
