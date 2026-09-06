# The `ImportReport` is the structured side-output of `read_omics()`. The
# input file may have multiple sheets/tabs (or columns, in the CSV case);
# each is classified as "matrix" / "metadata" / "feature_annot" / "unknown"
# with a confidence score and a free-form note. The report also carries:
#   * `warnings`  - human-readable strings the UI surfaces verbatim
#   * `suggested_input` - the names of the sheets we picked, the chosen
#                          orientation, and the ID column we found, so the
#                          Shiny wizard can pre-populate its controls.

#' Sheet role labels recognised by [read_omics()] and the classifier
#'
#' @export
#' @keywords internal
IMPORT_REPORT_ROLES <- c("matrix", "metadata", "feature_annot", "unknown")

#' Construct an `ImportReport`
#'
#' Side-output object emitted by [read_omics()]. Bundles a per-sheet
#' classification table, free-form warnings, and a `suggested_input` slot
#' that captures the assignment used to (try to) build an `omics_input`.
#'
#' @param sheets A `data.frame` (or `tibble`) with one row per sheet/tab
#'   detected in the input file. Required columns: `name`, `role`,
#'   `n_rows`, `n_cols`, `confidence`, `orientation`, `notes`. The
#'   constructor coerces partial inputs and fills in missing columns.
#' @param warnings Character vector of warnings to surface to the user.
#' @param suggested_input Named list describing the assignment used to
#'   build the returned `omics_input`. Conventional fields:
#'   `matrix_sheet`, `metadata_sheet`, `feature_sheet`, `orientation`,
#'   `id_column`, `omics_type`, `assay_type`. May be empty.
#' @param source Optional file path or other origin label.
#'
#' @return An object of class `ImportReport`.
#' @export
#' @family io
new_import_report <- function(
  sheets = NULL,
  warnings = character(0),
  suggested_input = list(),
  source = NA_character_
) {
  sheets <- coerce_sheet_table(sheets)
  if (!is.null(warnings) && !is.atomic(warnings)) {
    arg_stop("warnings", "a character vector", warnings)
  }
  if (!is.character(warnings)) warnings <- as.character(warnings)
  if (is.null(suggested_input)) suggested_input <- list()
  if (!is.list(suggested_input)) {
    stop("`suggested_input` must be a named list.")
  }
  structure(
    list(
      sheets = sheets,
      warnings = warnings,
      suggested_input = suggested_input,
      source = source
    ),
    class = "ImportReport"
  )
}

#' Test whether an object is an `ImportReport`
#'
#' @param x Object to test.
#'
#' @return Logical scalar.
#' @export
#' @family io
is_import_report <- function(x) {
  inherits(x, "ImportReport")
}

#' Print method for `ImportReport`
#'
#' @param x An `ImportReport`.
#' @param ... Unused.
#'
#' @return Invisibly returns `x`.
#' @export
print.ImportReport <- function(x, ...) {
  cat("<ImportReport>\n")
  if (!is.na(x$source)) {
    cat("  source   :", x$source, "\n")
  }
  cat("  sheets   :", nrow(x$sheets), "\n")
  if (nrow(x$sheets) > 0L) {
    for (i in seq_len(nrow(x$sheets))) {
      r <- x$sheets[i, , drop = FALSE]
      cat(sprintf("    - %-20s role=%-13s %dx%d  conf=%.2f  orient=%s\n",
                  r$name, r$role, as.integer(r$n_rows %||% 0L),
                  as.integer(r$n_cols %||% 0L),
                  as.numeric(r$confidence %||% NA_real_),
                  r$orientation %||% "NA"))
    }
  }
  if (length(x$warnings) > 0L) {
    cat("  warnings :", length(x$warnings), "\n")
    for (w in x$warnings) {
      cat("    !", w, "\n")
    }
  }
  if (length(x$suggested_input) > 0L) {
    cat("  suggested_input:\n")
    for (nm in names(x$suggested_input)) {
      val <- x$suggested_input[[nm]]
      cat("    -", nm, ":", format_short(val), "\n")
    }
  }
  invisible(x)
}

#' Access the sheets table from an `ImportReport`
#'
#' @param report An `ImportReport`.
#'
#' @return A `data.frame` of detected sheets.
#' @export
#' @family io
import_report_sheets <- function(report) {
  if (!is_import_report(report)) arg_stop("report", "an `ImportReport`", report)
  report$sheets
}

#' Access the warnings vector from an `ImportReport`
#'
#' @param report An `ImportReport`.
#'
#' @return Character vector of warnings.
#' @export
#' @family io
import_report_warnings <- function(report) {
  if (!is_import_report(report)) arg_stop("report", "an `ImportReport`", report)
  report$warnings
}

#' Append a warning to an `ImportReport`
#'
#' @param report An `ImportReport`.
#' @param msg Character scalar warning message.
#'
#' @return The updated `ImportReport`.
#' @keywords internal
add_import_warning <- function(report, msg) {
  stopifnot(is_import_report(report))
  report$warnings <- c(report$warnings, as.character(msg))
  report
}

# ---- internal helpers --------------------------------------------------

# The full sheets-table schema. Kept here so both the constructor and the
# classifier produce identically-shaped data.frames.
import_report_sheets_template <- function() {
  data.frame(
    name = character(0),
    role = character(0),
    n_rows = integer(0),
    n_cols = integer(0),
    confidence = numeric(0),
    orientation = character(0),
    notes = character(0),
    stringsAsFactors = FALSE
  )
}

coerce_sheet_table <- function(sheets) {
  tmpl <- import_report_sheets_template()
  if (is.null(sheets)) {
    return(tmpl)
  }
  if (!is.data.frame(sheets)) {
    stop("`sheets` must be a data.frame.")
  }
  for (col in setdiff(colnames(tmpl), colnames(sheets))) {
    sheets[[col]] <- switch(
      col,
      n_rows = integer(nrow(sheets)),
      n_cols = integer(nrow(sheets)),
      confidence = rep(NA_real_, nrow(sheets)),
      rep(NA_character_, nrow(sheets))
    )
  }
  sheets <- sheets[, colnames(tmpl), drop = FALSE]
  sheets$name <- as.character(sheets$name)
  sheets$role <- as.character(sheets$role)
  sheets$n_rows <- as.integer(sheets$n_rows)
  sheets$n_cols <- as.integer(sheets$n_cols)
  sheets$confidence <- as.numeric(sheets$confidence)
  sheets$orientation <- as.character(sheets$orientation)
  sheets$notes <- as.character(sheets$notes)
  rownames(sheets) <- NULL
  sheets
}

format_short <- function(x) {
  if (is.null(x)) return("NULL")
  if (is.character(x) && length(x) == 1L) return(x)
  if (is.atomic(x) && length(x) <= 4L) {
    return(paste(format(x), collapse = ", "))
  }
  paste0("<", class(x)[1L], "[", length(x), "]>")
}
