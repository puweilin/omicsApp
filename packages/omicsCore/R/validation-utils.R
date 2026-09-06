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


# ---- argument contracts ------------------------------------------------
#
# Every public entry point answers a wrong argument in one voice: the
# argument by name, what it had to be, and what arrived. Without these
# a threshold given as a string reached `<` and answered "non-numeric
# argument to binary operator", and a column name given as a list
# reached `%in%` and answered "argument is of length zero" -- messages
# that name nothing the caller typed. Found by calling every exported
# function with one wrong argument at a time.

describe_value <- function(x) {
  if (is.null(x)) return("NULL")
  if (is.function(x)) return("a function")
  if (is.data.frame(x)) return(sprintf("a data.frame with %d row(s)", nrow(x)))
  if (is.matrix(x)) return(sprintf("a %d x %d matrix", nrow(x), ncol(x)))
  if (is.environment(x)) return("an environment")
  if (isS4(x)) return(paste0("an S4 object of class ", class(x)[[1L]]))
  cls <- setdiff(class(x), c("numeric", "character", "logical", "integer",
                             "double", "list", "complex"))
  if (length(cls)) return(paste0("an object of class ", cls[[1L]]))
  if (is.list(x)) return(sprintf("a list of length %d", length(x)))
  if (length(x) == 0L) return(sprintf("an empty %s vector", typeof(x)))
  if (length(x) > 1L) return(sprintf("a %s vector of length %d", typeof(x), length(x)))
  if (is.na(x)) return("NA")
  if (is.character(x)) return(sprintf("'%s'", x))
  format(x)
}

arg_stop <- function(arg, expected, x) {
  stop(sprintf("`%s` must be %s, not %s.", arg, expected, describe_value(x)),
       call. = FALSE)
}

assert_string <- function(x, arg, allow_null = FALSE, allow_empty = FALSE) {
  if (is.null(x) && allow_null) return(invisible(NULL))
  ok <- is.character(x) && length(x) == 1L && !is.na(x) && (allow_empty || nzchar(x))
  if (!ok) {
    arg_stop(arg, if (allow_empty) "a single string" else "a single non-empty string", x)
  }
  invisible(x)
}

# A group level is whatever a metadata column holds, so a number is as
# valid as a string.
assert_label <- function(x, arg, allow_null = FALSE) {
  if (is.null(x) && allow_null) return(invisible(NULL))
  ok <- (is.character(x) || is.numeric(x) || is.factor(x)) &&
    length(x) == 1L && !is.na(x)
  if (!ok) arg_stop(arg, "a single value (a string or a number)", x)
  invisible(x)
}

# Any number of names, including none.
assert_character <- function(x, arg, allow_null = FALSE) {
  if (is.null(x) && allow_null) return(invisible(NULL))
  if (is.factor(x)) x <- as.character(x)
  ok <- is.character(x) && !anyNA(x)
  if (!ok) arg_stop(arg, "a character vector", x)
  invisible(x)
}

# At least one name.
assert_names <- function(x, arg, allow_null = FALSE) {
  if (is.null(x) && allow_null) return(invisible(NULL))
  if (is.factor(x)) x <- as.character(x)
  ok <- is.character(x) && length(x) >= 1L && !anyNA(x) && all(nzchar(x))
  if (!ok) arg_stop(arg, "a character vector of one or more non-empty names", x)
  invisible(x)
}

assert_choice <- function(x, arg, choices, allow_null = FALSE) {
  if (is.null(x) && allow_null) return(invisible(NULL))
  ok <- is.character(x) && length(x) == 1L && !is.na(x) && x %in% choices
  if (!ok) {
    arg_stop(arg, paste0("one of ", paste0("'", choices, "'", collapse = ", ")), x)
  }
  invisible(x)
}

assert_subset <- function(x, arg, choices, allow_null = FALSE) {
  if (is.null(x) && allow_null) return(invisible(NULL))
  ok <- is.character(x) && length(x) >= 1L && !anyNA(x) && all(x %in% choices)
  if (!ok) {
    arg_stop(arg, paste0("one or more of ",
                         paste0("'", choices, "'", collapse = ", ")), x)
  }
  invisible(x)
}

assert_number <- function(x, arg, lower = -Inf, upper = Inf, allow_null = FALSE) {
  if (is.null(x) && allow_null) return(invisible(NULL))
  range_txt <- if (is.finite(lower) && is.finite(upper)) {
    sprintf(" between %s and %s", format(lower), format(upper))
  } else if (is.finite(lower)) {
    sprintf(" of at least %s", format(lower))
  } else if (is.finite(upper)) {
    sprintf(" of at most %s", format(upper))
  } else ""
  ok <- is.numeric(x) && length(x) == 1L && !is.na(x) && x >= lower && x <= upper
  if (!ok) arg_stop(arg, paste0("a single number", range_txt), x)
  invisible(x)
}

assert_count <- function(x, arg, lower = 0L, allow_null = FALSE) {
  if (is.null(x) && allow_null) return(invisible(NULL))
  ok <- is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x) &&
    x >= lower && x == floor(x)
  if (!ok) arg_stop(arg, sprintf("a single whole number of at least %d", as.integer(lower)), x)
  invisible(x)
}

assert_flag <- function(x, arg, allow_null = FALSE) {
  if (is.null(x) && allow_null) return(invisible(NULL))
  ok <- is.logical(x) && length(x) == 1L && !is.na(x)
  if (!ok) arg_stop(arg, "TRUE or FALSE", x)
  invisible(x)
}

assert_data_frame <- function(x, arg, allow_null = FALSE) {
  if (is.null(x) && allow_null) return(invisible(NULL))
  if (!is.data.frame(x)) arg_stop(arg, "a data.frame", x)
  invisible(x)
}

assert_list <- function(x, arg, allow_null = FALSE) {
  if (is.null(x) && allow_null) return(invisible(NULL))
  if (!is.list(x) || is.data.frame(x)) arg_stop(arg, "a list", x)
  invisible(x)
}

# Returns the matrix: a data.frame of numbers is accepted and handed
# back as one.
assert_numeric_matrix <- function(x, arg) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x) || !is.numeric(x)) arg_stop(arg, "a numeric matrix", x)
  x
}

assert_project <- function(x, arg = "project") {
  if (!is_omics_project(x)) arg_stop(arg, "an `omics_project`", x)
  invisible(x)
}
