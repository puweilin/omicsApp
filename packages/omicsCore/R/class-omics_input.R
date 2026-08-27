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
#' @param assay_type Assay semantic label, ideally one of the values listed for
#'   this `omics_type` in [SUPPORTED_ASSAY_TYPES] (e.g. `"raw_intensity"`,
#'   `"normalized_intensity"`, `"raw_count"`). This is the only record of what
#'   scale `expr_mat` is on -- nothing downstream re-derives it -- so an
#'   inaccurate label silently changes what the analysis backends do. Superseded
#'   spellings in [DEPRECATED_ASSAY_TYPE_ALIASES] are rewritten with a warning.
#' @param raw_mat Optional raw matrix carried alongside `expr_mat`.
#' @param normalized_mat Optional normalized matrix carried alongside `expr_mat`.
#' @param raw_object Optional upstream raw object (e.g. a `DESeqDataSet` or
#'   `SummarizedExperiment`) preserved for reference.
#' @param source_fingerprint Optional string identifying the file this
#'   input was parsed from. Callers that re-import into a live project use
#'   it to tell "the user picked the same file again" from "the data
#'   changed"; analyses computed on the previous data are only valid for
#'   the former. `NULL` for inputs built directly from matrices.
#' @param source_path Optional path to the archived copy of that file.
#'   [export_script()] points its `read_omics()` line at it, which is
#'   what makes an exported script runnable rather than illustrative.
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
#' omics_input(expr, meta, feat, omics_type = "proteomics",
#'             assay_type = "normalized_intensity")
omics_input <- function(
  expr_mat,
  meta_df,
  feature_df,
  omics_type,
  assay_type = NULL,
  raw_mat = NULL,
  normalized_mat = NULL,
  raw_object = NULL,
  source_fingerprint = NULL,
  source_path = NULL
) {
  assay_type <- canonical_assay_type(assay_type)

  x <- new_omics_input(
    omics_type = omics_type,
    assay_type = assay_type,
    expr_mat = expr_mat,
    meta_df = meta_df,
    feature_df = feature_df,
    raw_mat = raw_mat,
    normalized_mat = normalized_mat,
    raw_object = raw_object,
    source_fingerprint = source_fingerprint,
    source_path = source_path
  )
  validate_omics_input(x)
  # Once, at construction -- see check_assay_scale() for why not in validate
  check_assay_scale(x)
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
  raw_object = NULL,
  source_fingerprint = NULL,
  source_path = NULL
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
      raw_object = raw_object,
      source_fingerprint = source_fingerprint,
      source_path = source_path
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

#' Warn when an assay type is unknown for this modality
#'
#' Mirrors how [validate_omics_input()] treats `omics_type`: a warning, not an
#' error, so a modality omicsCore ships no dispatcher for still works.
#'
#' @param assay_type Assay type label.
#' @param omics_type Omics modality.
#'
#' @return Invisibly `TRUE`.
#' @keywords internal
validate_assay_type <- function(assay_type, omics_type) {
  if (is.null(assay_type) || !is.character(assay_type) ||
      length(assay_type) != 1L || is.na(assay_type) || !nzchar(assay_type)) {
    warning(
      "`assay_type` is missing. It is the only record of what scale ",
      "`expr_mat` is on; analyses read it to decide whether to transform. ",
      "Set one of: ",
      paste(sprintf("'%s'", unlist(SUPPORTED_ASSAY_TYPES, use.names = FALSE)),
            collapse = ", "),
      call. = FALSE
    )
    return(invisible(TRUE))
  }

  if (assay_type %in% names(DEPRECATED_ASSAY_TYPE_ALIASES)) {
    warning(
      "`assay_type = \"", assay_type, "\"` is superseded by \"",
      unname(DEPRECATED_ASSAY_TYPE_ALIASES[[assay_type]]),
      "\". Rebuild this input with `omics_input()` to update it.",
      call. = FALSE
    )
    return(invisible(TRUE))
  }

  expected <- SUPPORTED_ASSAY_TYPES[[omics_type]]
  if (!is.null(expected) && !assay_type %in% expected) {
    warning(
      "Unrecognised `assay_type` for ", omics_type, ": '", assay_type,
      "'. Known types are ", paste(sprintf("'%s'", expected), collapse = ", "),
      ". Analyses will treat the values as-is without transforming them.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Rewrite a superseded assay-type spelling
#'
#' Maps the names in [DEPRECATED_ASSAY_TYPE_ALIASES] onto their replacements,
#' warning so the caller can update. Anything else passes through untouched,
#' including labels outside [SUPPORTED_ASSAY_TYPES].
#'
#' @param assay_type Assay type label, or `NULL`.
#'
#' @return The canonical label, or `assay_type` unchanged.
#' @keywords internal
canonical_assay_type <- function(assay_type) {
  if (is.null(assay_type) || !is.character(assay_type) ||
      length(assay_type) != 1L || is.na(assay_type)) {
    return(assay_type)
  }
  if (assay_type %in% names(DEPRECATED_ASSAY_TYPE_ALIASES)) {
    replacement <- unname(DEPRECATED_ASSAY_TYPE_ALIASES[[assay_type]])
    warning(
      "`assay_type = \"", assay_type, "\"` is superseded; using \"",
      replacement, "\". Update the caller to the new spelling.",
      call. = FALSE
    )
    return(replacement)
  }
  assay_type
}

#' Guess an assay type from the values
#'
#' A starting point for an import wizard, not a verdict: it separates linear
#' proteomics intensities from already-transformed ones by magnitude, the same
#' heuristic [check_assay_scale()] uses, and assumes counts for RNA-seq because
#' that is what gets uploaded. The caller is expected to show the guess and let
#' the user correct it -- nothing else recovers the scale if this is wrong.
#'
#' @param expr_mat Expression matrix.
#' @param omics_type Omics modality.
#'
#' @return A single assay type from [SUPPORTED_ASSAY_TYPES], or `NA_character_`
#'   for a modality with no vocabulary.
#' @export
#' @family omics_input
#' @examples
#' infer_assay_type(matrix(2^rnorm(40, 20, 2), nrow = 10), "proteomics")
#' infer_assay_type(matrix(rpois(40, 200), nrow = 10), "rnaseq")
infer_assay_type <- function(expr_mat, omics_type) {
  if (identical(omics_type, "rnaseq")) return("raw_count")
  if (!identical(omics_type, "proteomics")) return(NA_character_)

  max_value <- suppressWarnings(max(as.matrix(expr_mat), na.rm = TRUE))
  if (!is.finite(max_value)) return(NA_character_)

  if (max_value > MAX_PLAUSIBLE_LOG_SCALE_VALUE) {
    "raw_intensity"
  } else {
    "normalized_intensity"
  }
}

#' Check that the assay values look like the declared scale
#'
#' `assay_type` is the only record of what scale `expr_mat` is on, and the
#' analysis backends act on it without re-deriving anything: limma runs
#' straight on `expr_mat`, while `log2(x + 1)` is applied only for
#' `"raw_count"`. Mislabelling therefore changes results silently -- feeding
#' vsn-normalised values to a `raw_intensity` path, or running limma on linear
#' intensities, produces plausible numbers that are wrong.
#'
#' The check is a magnitude heuristic: log2 proteomics intensities top out
#' around 35, linear intensities run to millions. It warns rather than fails,
#' because a small or unusual matrix can legitimately land either side of
#' [MAX_PLAUSIBLE_LOG_SCALE_VALUE].
#'
#' Not called from [validate_omics_input()] on purpose: validation runs at
#' every analysis entry point, and this scans the matrix. [omics_input()] runs
#' it once at construction, which is where a wrong label enters.
#'
#' @param x An `omics_input`.
#'
#' @return Invisibly, `TRUE` when the values match the declared scale,
#'   `FALSE` when they do not.
#' @export
#' @family omics_input
check_assay_scale <- function(x) {
  if (!is_omics_input(x)) stop("Object is not an `omics_input`.")

  assay_type <- x$assay_type
  if (is.null(assay_type) || !is.character(assay_type) ||
      length(assay_type) != 1L || is.na(assay_type) ||
      !assay_type %in% SCALE_CHECKED_ASSAY_TYPES) {
    # Either the label's scale is unknown, or magnitude cannot resolve it --
    # see SCALE_CHECKED_ASSAY_TYPES for why counts are excluded
    return(invisible(TRUE))
  }

  max_value <- suppressWarnings(max(x$expr_mat, na.rm = TRUE))
  if (!is.finite(max_value)) return(invisible(TRUE))

  expects_log <- assay_type %in% LOG_SCALE_ASSAY_TYPES
  looks_log <- max_value <= MAX_PLAUSIBLE_LOG_SCALE_VALUE

  if (expects_log && !looks_log) {
    warning(
      "`assay_type = \"", assay_type, "\"` implies log-scale values, but the ",
      "matrix reaches ", format(max_value, digits = 4), ". If these are linear ",
      "intensities, normalize them first (see `normalize_omics()`); analyses ",
      "will otherwise run on untransformed data.",
      call. = FALSE
    )
    return(invisible(FALSE))
  }
  if (!expects_log && looks_log) {
    warning(
      "`assay_type = \"", assay_type, "\"` implies linear values, but the ",
      "matrix only reaches ", format(max_value, digits = 4), ", which looks ",
      "already log-transformed. Transforming again would compress the dynamic ",
      "range and shrink every fold change.",
      call. = FALSE
    )
    return(invisible(FALSE))
  }
  invisible(TRUE)
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

  validate_assay_type(x$assay_type, omics_type)

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
