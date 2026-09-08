# One workbook, one sheet per analysis.
#
# What a reader of this project does next is open the tables in Excel:
# sort the differential by p, filter the pathways, paste a column into a
# figure. A rendered report answered a different question -- "what did
# this run do" -- and answered it in a form nobody could carry into the
# next step. This produces the tables.

WORKBOOK_SHEET_LIMIT <- 31L   # Excel's, not ours

#' Excel's rules for a sheet name, applied
#'
#' At most 31 characters, and none of `[ ] : * ? / \`. openxlsx enforces
#' this by stopping, which would fail a whole export over a project name.
#'
#' @keywords internal
#' @noRd
safe_sheet_name <- function(name, taken = character(0)) {
  out <- gsub("[\\[\\]:*?/\\\\]", "-", name, perl = TRUE)
  out <- trimws(out)
  if (!nzchar(out)) out <- "Sheet"
  out <- substr(out, 1L, WORKBOOK_SHEET_LIMIT)
  # Excel compares sheet names case-insensitively.
  if (tolower(out) %in% tolower(taken)) {
    for (i in 2:99) {
      suffix <- paste0(" (", i, ")")
      cand <- paste0(substr(out, 1L, WORKBOOK_SHEET_LIMIT - nchar(suffix)),
                     suffix)
      if (!tolower(cand) %in% tolower(taken)) return(cand)
    }
  }
  out
}

#' The expression matrix as a sheet, with its feature ids in a column
#'
#' A matrix written straight out loses its row names -- the reader gets
#' 129 columns of numbers and no way to tell which gene is which.
#'
#' @keywords internal
#' @noRd
matrix_sheet <- function(mat, id_name = "feature_id") {
  df <- as.data.frame(mat, stringsAsFactors = FALSE, check.names = FALSE)
  ids <- rownames(mat)
  if (is.null(ids)) ids <- as.character(seq_len(nrow(mat)))
  cbind(stats::setNames(data.frame(ids, stringsAsFactors = FALSE), id_name),
        df)
}

#' The parameters each analysis ran with, as one long table
#'
#' The workbook replaces a report, and the report carried this. Without it
#' a sheet of differential results does not say which contrast it is, which
#' method produced it, or what was filtered upstream -- and that is exactly
#' what someone asks three months later.
#'
#' @keywords internal
#' @noRd
params_sheet <- function(project) {
  rows <- list()
  add <- function(analysis, key, value) {
    rows[[length(rows) + 1L]] <<- data.frame(
      analysis = analysis, parameter = key,
      value = describe_param_value(value),
      stringsAsFactors = FALSE)
  }
  add("project", "name", project$name %||% NA_character_)
  add("project", "exported", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  pp <- project$preprocessing %||% list()
  for (tag in names(project$experiments %||% list())) {
    e <- project$experiments[[tag]]
    lab <- paste0("layer: ", tag)
    add(lab, "omics_type", e$omics_type %||% NA_character_)
    add(lab, "assay_type", e$assay_type %||% NA_character_)
    add(lab, "features (as imported)", nrow(e$expr_mat))
    add(lab, "samples (as imported)", ncol(e$expr_mat))
    # As imported and as analysed are different numbers once QC has been
    # applied, and the sheet used to give only the first -- so a reader
    # comparing it against a differential of 26,492 features had no way to
    # account for the other 36,749.
    prep <- pp[[tag]]
    if (!is.null(prep$cleaned)) {
      add(lab, "features (analysed)", nrow(prep$cleaned$expr_mat))
      add(lab, "samples (analysed)", ncol(prep$cleaned$expr_mat))
      for (k in names(prep$params %||% list())) {
        add(paste0("QC: ", tag), k, prep$params[[k]])
      }
    }
  }
  for (nm in names(project$bundles %||% list())) {
    b <- project$bundles[[nm]]
    # Which layer and which run, not only the settings. A workbook whose
    # Differential sheet does not say what it was computed on is a table
    # of numbers, and a project with two layers of one omics type has
    # nothing else to tell them apart.
    add(nm, "layer_tag",
        paste(b$input_info$layer_tag %||% NA_character_, collapse = ", "))
    add(nm, "analysis_id", b$provenance_binding$analysis_id %||% NA_character_)
    for (k in names(b$params %||% list())) add(nm, k, b$params[[k]])
  }
  # The results a displayed one depends on but which no slot shows. The
  # integration's second differential is registered and never displayed,
  # so a two-layer workbook described one side of a comparison.
  shown <- vapply(project$bundles %||% list(),
                  function(b) b$provenance_binding$analysis_id %||% "",
                  character(1))
  for (id in names(project$analyses %||% list())) {
    if (id %in% shown) next
    b <- project$analyses[[id]]
    lab <- paste0("upstream: ", b$analysis_name %||% "analysis")
    add(lab, "analysis_id", id)
    add(lab, "layer_tag",
        paste(b$input_info$layer_tag %||% NA_character_, collapse = ", "))
    for (k in names(b$params %||% list())) add(lab, k, b$params[[k]])
  }
  if (length(rows) == 0L) return(NULL)
  do.call(rbind, rows)
}

# A parameter as one cell. Lists and vectors reach here from `covariates`
# and from the imputation backends' extra arguments.
describe_param_value <- function(x) {
  if (is.null(x)) return(NA_character_)
  if (is.list(x)) x <- unlist(x, use.names = FALSE)
  if (length(x) == 0L) return(NA_character_)
  paste(format(x, trim = TRUE), collapse = ", ")
}

#' The sheets a project has to offer, in the order they were produced
#'
#' Returned as a list rather than written, so the caller can say what a
#' workbook would contain without building one -- the Report view shows
#' that before the download, and computing it twice by two routes is how
#' the two come to disagree.
#'
#' @param project An [`omics_project`][is_omics_project()].
#' @return A named list of data frames, possibly empty.
#' @export
#' @family persistence
workbook_sheets <- function(project) {
  if (!is_omics_project(project)) {
    stop("`project` must be an `omics_project`.")
  }
  bundles <- project$bundles %||% list()
  sheets <- list()

  # QC, and only when it changed the values. An "imputed" sheet identical
  # to the import is a sheet that invites someone to believe imputation
  # ran. Missingness and outlier tables go out whenever QC ran, because
  # those are the QC result whether or not anything was imputed.
  qc <- bundles$qc
  if (!is.null(qc)) {
    miss <- qc$results$qc_summary$missingness
    if (!is.null(miss$sample_metrics)) {
      sheets[["QC samples"]] <- miss$sample_metrics
    }
    if (!is.null(miss$feature_metrics)) {
      sheets[["QC features"]] <- miss$feature_metrics
    }
    cleaned <- qc$results$cleaned_input
    imputed <- !identical(qc$params$impute_method %||% "none", "none")
    if (imputed && !is.null(cleaned$expr_mat)) {
      sheets[["Imputed data"]] <- matrix_sheet(cleaned$expr_mat)
    }
  }

  if (!is.null(bundles$diff$results$diff_result_df)) {
    sheets[["Differential"]] <- bundles$diff$results$diff_result_df
  }

  # The differentials no slot displays. run_integration() works from two,
  # and only one of them is ever shown -- the second is computed by the
  # Integration view and registered, never displayed. A two-layer
  # workbook therefore described one side of a comparison and offered the
  # concordance as though both were in front of the reader.
  #
  # Named by layer, because "Differential (2)" says nothing.
  shown_ids <- vapply(bundles, function(b) {
    b$provenance_binding$analysis_id %||% ""
  }, character(1))
  for (id in names(project$analyses %||% list())) {
    if (id %in% shown_ids) next
    b <- project$analyses[[id]]
    df <- b$results$diff_result_df
    if (is.null(df) || nrow(df) == 0L) next
    tag <- (b$input_info$layer_tag %||% NA_character_)[[1L]]
    nm <- if (is.na(tag)) "Differential (upstream)" else
      paste0("Differential (", tag, ")")
    if (nm %in% names(sheets)) next
    sheets[[nm]] <- df
  }
  if (!is.null(bundles$enrich$results$enrich_result_df)) {
    sheets[["Enrichment"]] <- bundles$enrich$results$enrich_result_df
  }

  # The full grid, split by test rather than by database. ORA and GSEA
  # answer different questions and are read separately; databases are a
  # column within each, which is how the rows already carry them -- one
  # sheet per database would be fourteen sheets and no more information.
  #
  # The log goes out too, because "GO BP is not in this workbook" and "GO
  # BP came back with nothing" look identical from the sheets alone, and
  # only one of them means something went wrong.
  grid <- bundles$enrich_all
  if (!is.null(grid$results$enrich_result_df)) {
    all_rows <- grid$results$enrich_result_df
    if ("result_type" %in% names(all_rows)) {
      for (rt in unique(all_rows$result_type)) {
        part <- all_rows[all_rows$result_type == rt, , drop = FALSE]
        if (nrow(part) == 0L) next
        sheets[[paste("Enrichment", toupper(rt))]] <- part
      }
    } else if (nrow(all_rows) > 0L) {
      sheets[["Enrichment (all)"]] <- all_rows
    }
  }
  if (!is.null(grid$results$grid_log)) {
    sheets[["Enrichment log"]] <- grid$results$grid_log
  }
  if (!is.null(bundles$integration$results$integration_df)) {
    sheets[["Integration"]] <- bundles$integration$results$integration_df
  }

  # The parameters go in front, and only when there is something in front
  # of. On its own it is not a workbook -- it is a description of a project
  # nobody has analysed yet -- and offering that as a download makes the
  # button look like it worked. Empty here is what disables it.
  if (length(sheets) == 0L) return(list())
  p <- params_sheet(project)
  if (is.null(p)) return(sheets)
  c(list(`Analysis parameters` = p), sheets)
}

#' Write a project's tables to one Excel workbook
#'
#' One sheet per analysis: the parameters everything ran with, QC's
#' per-sample and per-feature metrics, the imputed matrix when imputation
#' actually ran, and the differential, enrichment and integration tables
#' that exist. Analyses that were not run contribute no sheet rather than
#' an empty one.
#'
#' @param project An [`omics_project`][is_omics_project()].
#' @param path Output path.
#' @param overwrite If `FALSE` (default) and `path` exists, raise an error.
#'
#' @return Invisibly, `path`.
#' @export
#' @family persistence
#' @examples
#' \dontrun{
#'   export_workbook(project, "results.xlsx")
#' }
export_workbook <- function(project, path, overwrite = FALSE) {
  assert_flag(overwrite, "overwrite")
  if (!is_omics_project(project)) {
    stop("`project` must be an `omics_project`.")
  }
  if (!is.character(path) || length(path) != 1L || !nzchar(path)) {
    stop("`path` must be a non-empty single string.")
  }
  if (file.exists(path) && !isTRUE(overwrite)) {
    stop("Path already exists: ", path,
         " (pass `overwrite = TRUE` to replace).")
  }
  if (!is_installed("openxlsx")) {
    stop("Package 'openxlsx' is required for export_workbook(). ",
         "Install with: install.packages('openxlsx').", call. = FALSE)
  }

  sheets <- workbook_sheets(project)
  if (length(sheets) == 0L) {
    stop("This project has nothing to export yet. Run QC, a differential, ",
         "an enrichment or an integration first.", call. = FALSE)
  }

  wb <- openxlsx::createWorkbook()
  used <- character(0)
  for (nm in names(sheets)) {
    sheet <- safe_sheet_name(nm, used)
    used <- c(used, sheet)
    openxlsx::addWorksheet(wb, sheet)
    openxlsx::writeData(wb, sheet, sheets[[nm]])
    openxlsx::freezePane(wb, sheet, firstRow = TRUE)
  }

  # Written beside the target and moved into place, so a failure part way
  # through leaves the previous file rather than a truncated one. The
  # download handler serves whatever is at `path`.
  tmp <- paste0(path, ".part")
  on.exit(unlink(tmp), add = TRUE)
  openxlsx::saveWorkbook(wb, tmp, overwrite = TRUE)
  if (!file.copy(tmp, path, overwrite = TRUE)) {
    stop("Could not write the workbook to ", path, call. = FALSE)
  }

  invisible(path)
}
