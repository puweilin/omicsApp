#' Single-omics input container
#'
#' Constructs the canonical input object consumed by all `omicsCore` analysis
#' functions. The argument order is `expr_mat`, `meta_df`, `feature_df`,
#' `omics_type`, `assay_type` to match the contract in
#' `docs/export-manifest.md`.
#'
#' `expr_mat` must be a numeric matrix with sample column names and feature
#' row names. `meta_df` must have rownames matching the sample column names of
#' `expr_mat`. `feature_df` must contain a `feature_id` column.
#'
#' @param expr_mat Main expression matrix; rows are features, columns are
#'   samples. A `data.frame` is coerced to a matrix.
#' @param meta_df Sample metadata; `rownames(meta_df)` must include
#'   `colnames(expr_mat)`.
#' @param feature_df Feature metadata containing a `feature_id` column.
#' @param omics_type Omics modality, one of [SUPPORTED_OMICS_TYPES].
#' @param assay_type Optional assay semantic label, e.g. `"normalized_intensity"`,
#'   `"raw_counts"`.
#' @param raw_mat Optional raw matrix carried alongside `expr_mat`.
#' @param normalized_mat Optional normalized matrix carried alongside `expr_mat`.
#' @param raw_object Optional upstream raw object (e.g. a `DESeqDataSet` or
#'   `SummarizedExperiment`) preserved for reference.
#'
#' @return An object of class `omics_input`.
#' @export
#' @family omics_input
#' @examples
#' expr <- matrix(rnorm(20), nrow = 4, dimnames = list(
#'   paste0("g", 1:4), paste0("s", 1:5)
#' ))
#' meta <- data.frame(group = c("A", "A", "B", "B", "B"),
#'                    row.names = paste0("s", 1:5))
#' feat <- data.frame(feature_id = paste0("g", 1:4),
#'                    row.names = paste0("g", 1:4))
#' omics_input(expr, meta, feat, omics_type = "proteomics")
omics_input <- function(
  expr_mat,
  meta_df,
  feature_df,
  omics_type,
  assay_type = NULL,
  raw_mat = NULL,
  normalized_mat = NULL,
  raw_object = NULL
) {
  x <- new_omics_input(
    omics_type = omics_type,
    assay_type = assay_type,
    expr_mat = expr_mat,
    meta_df = meta_df,
    feature_df = feature_df,
    raw_mat = raw_mat,
    normalized_mat = normalized_mat,
    raw_object = raw_object
  )
  validate_omics_input(x)
  x
}

#' Low-level `omics_input` constructor
#'
#' Builds an `omics_input` without running [validate_omics_input()]. Prefer
#' [omics_input()] in user-facing code; this constructor is used internally
#' when validation is deferred (e.g. inside subset operations).
#'
#' @inheritParams omics_input
#'
#' @return An object of class `omics_input`.
#' @keywords internal
new_omics_input <- function(
  omics_type,
  assay_type,
  expr_mat,
  meta_df,
  feature_df,
  raw_mat = NULL,
  normalized_mat = NULL,
  raw_object = NULL
) {
  if (!is.matrix(expr_mat) && !is.data.frame(expr_mat)) {
    stop("`expr_mat` must be a matrix or data.frame.")
  }
  if (!is.data.frame(meta_df)) {
    stop("`meta_df` must be a data.frame.")
  }
  if (!is.data.frame(feature_df)) {
    stop("`feature_df` must be a data.frame.")
  }

  expr_mat <- as.matrix(expr_mat)

  structure(
    list(
      omics_type = omics_type,
      assay_type = assay_type,
      expr_mat = expr_mat,
      raw_mat = raw_mat,
      normalized_mat = normalized_mat,
      meta_df = meta_df,
      feature_df = feature_df,
      raw_object = raw_object
    ),
    class = "omics_input"
  )
}

#' Test whether an object is an `omics_input`
#'
#' @param x Object to test.
#'
#' @return Logical scalar.
#' @export
#' @family omics_input
is_omics_input <- function(x) {
  inherits(x, "omics_input")
}

#' Validate an `omics_input`
#'
#' Checks structural consistency between matrices and metadata: dimension
#' matches, presence of row/column names, presence of `feature_id`, and
#' membership of sample IDs in `meta_df`.
#'
#' `omics_type` is checked against [SUPPORTED_OMICS_TYPES] but only
#' *warns* when it does not match. omicsCore is a general engine and
#' callers may legitimately drive it with a modality the shipped
#' analysis dispatchers do not know about yet. The warning exists to
#' catch typos (`"proteomcis"`), not to gate extensibility.
#'
#' @param x Object expected to inherit from `omics_input`.
#'
#' @return Invisibly returns `TRUE` on success; otherwise raises an error.
#' @export
#' @family omics_input
validate_omics_input <- function(x) {
  if (!is_omics_input(x)) {
    stop("Object is not an `omics_input`.")
  }

  required_fields <- c(
    "omics_type", "assay_type", "expr_mat", "meta_df", "feature_df"
  )
  missing_fields <- setdiff(required_fields, names(x))
  if (length(missing_fields) > 0) {
    stop("Missing fields in `omics_input`: ", paste(missing_fields, collapse = ", "))
  }

  omics_type <- x$omics_type
  if (!is.character(omics_type) || length(omics_type) != 1L ||
      is.na(omics_type) || !nzchar(omics_type)) {
    stop("`omics_type` must be a non-empty single string.")
  }
  if (!omics_type %in% SUPPORTED_OMICS_TYPES) {
    warning(
      "Unrecognised `omics_type`: '", omics_type, "'. Known types are ",
      paste(sprintf("'%s'", SUPPORTED_OMICS_TYPES), collapse = ", "),
      ". Analysis dispatchers may fall back to generic defaults.",
      call. = FALSE
    )
  }

  expr_mat <- x$expr_mat
  meta_df <- x$meta_df
  feature_df <- x$feature_df

  if (ncol(expr_mat) != nrow(meta_df)) {
    stop("`ncol(expr_mat)` must equal `nrow(meta_df)`.")
  }
  if (nrow(expr_mat) != nrow(feature_df)) {
    stop("`nrow(expr_mat)` must equal `nrow(feature_df)`.")
  }

  if (is.null(colnames(expr_mat))) {
    stop("`expr_mat` must have sample column names.")
  }
  if (is.null(rownames(expr_mat))) {
    stop("`expr_mat` must have feature row names.")
  }
  if (is.null(rownames(meta_df))) {
    stop("`meta_df` must have rownames matching sample IDs.")
  }

  if (!all(colnames(expr_mat) %in% rownames(meta_df))) {
    stop("All `expr_mat` columns must be present in `rownames(meta_df)`.")
  }

  if (!"feature_id" %in% colnames(feature_df)) {
    stop("`feature_df` must contain `feature_id`.")
  }

  invisible(TRUE)
}

#' Print method for `omics_input`
#'
#' @param x An `omics_input` object.
#' @param ... Unused.
#'
#' @return Invisibly returns `x`.
#' @export
print.omics_input <- function(x, ...) {
  cat("<omics_input>\n")
  cat("  omics_type :", x$omics_type %||% "<unset>", "\n")
  cat("  assay_type :", x$assay_type %||% "<unset>", "\n")
  cat("  features   :", nrow(x$expr_mat), "\n")
  cat("  samples    :", ncol(x$expr_mat), "\n")
  cat("  missing %  :",
      sprintf("%.2f", mean(is.na(x$expr_mat)) * 100), "\n")
  invisible(x)
}

# Local null-coalescing helper. Kept package-local so we do not depend on
# rlang's `%||%` re-export through `@importFrom rlang %||%`.
`%||%` <- function(a, b) if (is.null(a)) b else a
