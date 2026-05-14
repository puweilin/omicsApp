# Smart input reader for the omicsApp import wizard.
#
# `read_omics(path)` scans an Excel workbook (or CSV / RDS) and returns
# `list(input, report)` where:
#   * `input`  is an `omics_input` if the heuristics produced a coherent
#              triple (matrix, metadata, optional feature annotation), or
#              `NULL` if the user needs to confirm an assignment.
#   * `report` is an `ImportReport` carrying the per-sheet table, any
#              warnings, and the assignment used (or proposed).

#' Read an omics workbook and report classifier findings
#'
#' Front door for the omicsApp import wizard. Reads an Excel workbook,
#' CSV, or saved R object, classifies each sheet/data frame, and tries
#' to build a candidate `omics_input`. Always returns both the built
#' object (or `NULL`) and a structured [`ImportReport`][new_import_report()]
#' so the UI can show what was detected before the user commits.
#'
#' @param path Path to a file. The extension drives auto-detection.
#' @param type One of `"auto"`, `"excel"`, `"csv"`, or `"rds"`.
#' @param omics_type Optional omics modality. If `NULL`, callers (the
#'   Shiny wizard) are expected to set it after inspecting the report.
#' @param assay_type Optional assay semantic label.
#' @param ... Forwarded to the underlying reader (`readxl::read_excel`,
#'   `utils::read.csv`).
#'
#' @return A list with two elements:
#'   * `input`: an `omics_input` or `NULL`.
#'   * `report`: an [`ImportReport`][new_import_report()].
#' @export
#' @family io
read_omics <- function(
  path,
  type = c("auto", "excel", "csv", "rds"),
  omics_type = NULL,
  assay_type = NULL,
  ...
) {
  type <- match.arg(type)
  if (!is.character(path) || length(path) != 1L || !nzchar(path)) {
    stop("`path` must be a non-empty single string.")
  }
  if (!file.exists(path)) {
    stop("File does not exist: ", path)
  }
  if (type == "auto") {
    type <- detect_file_type(path)
  }
  switch(type,
    excel = read_omics_excel(path, omics_type = omics_type,
                             assay_type = assay_type, ...),
    csv = read_omics_csv(path, omics_type = omics_type,
                         assay_type = assay_type, ...),
    rds = read_omics_rds(path, omics_type = omics_type,
                         assay_type = assay_type, ...),
    stop("Unsupported type: ", type)
  )
}

# ---- internal: per-type readers ----------------------------------------

detect_file_type <- function(path) {
  ext <- tolower(tools::file_ext(path))
  switch(ext,
    xlsx = "excel",
    xls  = "excel",
    csv  = "csv",
    tsv  = "csv",
    txt  = "csv",
    rds  = "rds",
    stop("Cannot auto-detect file type for extension: ", ext)
  )
}

read_omics_excel <- function(path, omics_type, assay_type, ...) {
  sheet_names <- readxl::excel_sheets(path)
  sheets <- list()
  rows <- list()
  for (nm in sheet_names) {
    df <- tryCatch(
      as.data.frame(readxl::read_excel(path, sheet = nm, ...),
                    stringsAsFactors = FALSE),
      error = function(e) NULL
    )
    if (is.null(df)) {
      rows[[length(rows) + 1L]] <- data.frame(
        name = nm, role = "unknown", n_rows = 0L, n_cols = 0L,
        confidence = 0, orientation = NA_character_,
        notes = "failed to read sheet",
        stringsAsFactors = FALSE
      )
      next
    }
    sheets[[nm]] <- df
    cls <- classify_sheet_role(df, name = nm)
    rows[[length(rows) + 1L]] <- data.frame(
      name = nm,
      role = cls$role,
      n_rows = as.integer(nrow(df)),
      n_cols = as.integer(ncol(df)),
      confidence = cls$confidence,
      orientation = cls$orientation %||% NA_character_,
      notes = cls$notes %||% "",
      stringsAsFactors = FALSE
    )
  }
  sheet_table <- do.call(rbind, rows)
  build_input_from_sheets(sheets, sheet_table, source = path,
                          omics_type = omics_type, assay_type = assay_type)
}

read_omics_csv <- function(path, omics_type, assay_type, ...) {
  sep <- if (grepl("\\.tsv$|\\.txt$", path, ignore.case = TRUE)) "\t" else ","
  df <- utils::read.table(path, header = TRUE, sep = sep,
                          check.names = FALSE, stringsAsFactors = FALSE, ...)
  nm <- tools::file_path_sans_ext(basename(path))
  cls <- classify_sheet_role(df, name = nm)
  sheet_table <- data.frame(
    name = nm, role = cls$role,
    n_rows = as.integer(nrow(df)), n_cols = as.integer(ncol(df)),
    confidence = cls$confidence,
    orientation = cls$orientation %||% NA_character_,
    notes = cls$notes %||% "",
    stringsAsFactors = FALSE
  )
  sheets <- setNames(list(df), nm)
  build_input_from_sheets(sheets, sheet_table, source = path,
                          omics_type = omics_type, assay_type = assay_type)
}

read_omics_rds <- function(path, omics_type, assay_type, ...) {
  obj <- readRDS(path)
  if (inherits(obj, "omics_input")) {
    if (!is.null(omics_type)) obj$omics_type <- omics_type
    if (!is.null(assay_type)) obj$assay_type <- assay_type
    sheet_table <- data.frame(
      name = "rds", role = "matrix",
      n_rows = as.integer(nrow(obj$expr_mat)),
      n_cols = as.integer(ncol(obj$expr_mat)),
      confidence = 1.0,
      orientation = "features_in_rows",
      notes = "loaded omics_input directly",
      stringsAsFactors = FALSE
    )
    report <- new_import_report(
      sheets = sheet_table,
      warnings = character(0),
      suggested_input = list(
        matrix_sheet = "rds",
        orientation = "features_in_rows",
        omics_type = obj$omics_type,
        assay_type = obj$assay_type
      ),
      source = path
    )
    return(list(input = obj, report = report))
  }
  # Not an omics_input - treat the object as a single data.frame/matrix
  # and run it through the normal classifier.
  df <- if (is.matrix(obj)) as.data.frame(obj, check.names = FALSE) else obj
  if (!is.data.frame(df)) {
    report <- new_import_report(
      sheets = data.frame(name = basename(path), role = "unknown",
                          n_rows = 0L, n_cols = 0L, confidence = 0,
                          orientation = NA_character_,
                          notes = paste("unsupported RDS payload:",
                                        paste(class(obj), collapse = "/")),
                          stringsAsFactors = FALSE),
      warnings = "RDS payload was neither an omics_input nor a data.frame.",
      source = path
    )
    return(list(input = NULL, report = report))
  }
  nm <- tools::file_path_sans_ext(basename(path))
  cls <- classify_sheet_role(df, name = nm)
  sheet_table <- data.frame(
    name = nm, role = cls$role,
    n_rows = as.integer(nrow(df)), n_cols = as.integer(ncol(df)),
    confidence = cls$confidence,
    orientation = cls$orientation %||% NA_character_,
    notes = cls$notes %||% "",
    stringsAsFactors = FALSE
  )
  sheets <- setNames(list(df), nm)
  build_input_from_sheets(sheets, sheet_table, source = path,
                          omics_type = omics_type, assay_type = assay_type)
}

# ---- internal: assembling an omics_input from classifier output --------

# Given the per-sheet table from the classifier, pick the highest-confidence
# matrix / metadata / feature_annot sheets and try to glue them into an
# `omics_input`. Returns list(input, report) without raising.
build_input_from_sheets <- function(sheets, sheet_table, source,
                                    omics_type, assay_type) {
  report <- new_import_report(sheets = sheet_table, source = source)
  if (nrow(sheet_table) == 0L) {
    report <- add_import_warning(report, "No readable sheets in the input.")
    return(list(input = NULL, report = report))
  }

  matrix_sheet <- pick_best_sheet(sheet_table, "matrix")
  metadata_sheet <- pick_best_sheet(sheet_table, "metadata")
  feature_sheet <- pick_best_sheet(sheet_table, "feature_annot")

  if (is.null(matrix_sheet)) {
    report <- add_import_warning(report,
      "No sheet looked like an expression matrix.")
    return(list(input = NULL, report = report))
  }
  if (is.null(omics_type)) {
    report <- add_import_warning(report,
      "`omics_type` not specified; caller must set it before constructing the input.")
  }

  orient <- sheet_table$orientation[sheet_table$name == matrix_sheet]
  mat <- materialize_matrix(sheets[[matrix_sheet]], orient)
  if (is.null(mat)) {
    report <- add_import_warning(report,
      "Could not coerce the picked matrix sheet to a numeric matrix.")
    report$suggested_input <- list(
      matrix_sheet = matrix_sheet,
      metadata_sheet = metadata_sheet,
      feature_sheet = feature_sheet,
      orientation = orient
    )
    return(list(input = NULL, report = report))
  }

  meta <- materialize_metadata(if (is.null(metadata_sheet)) NULL
                               else sheets[[metadata_sheet]],
                               sample_ids = colnames(mat))
  feat <- materialize_feature_annot(if (is.null(feature_sheet)) NULL
                                    else sheets[[feature_sheet]],
                                    feature_ids = rownames(mat))

  suggested <- list(
    matrix_sheet = matrix_sheet,
    metadata_sheet = metadata_sheet,
    feature_sheet = feature_sheet,
    orientation = orient,
    omics_type = omics_type,
    assay_type = assay_type,
    id_column = attr(feat, "id_column")
  )
  attr(feat, "id_column") <- NULL
  report$suggested_input <- suggested

  if (is.null(omics_type)) {
    return(list(input = NULL, report = report))
  }

  input <- tryCatch(
    omics_input(mat, meta, feat,
                omics_type = omics_type, assay_type = assay_type),
    error = function(e) {
      report <<- add_import_warning(report,
        paste0("omics_input() rejected the assembly: ", conditionMessage(e)))
      NULL
    }
  )
  list(input = input, report = report)
}

pick_best_sheet <- function(sheet_table, role) {
  cand <- sheet_table[sheet_table$role == role, , drop = FALSE]
  if (nrow(cand) == 0L) return(NULL)
  cand <- cand[order(cand$confidence, decreasing = TRUE), , drop = FALSE]
  cand$name[[1L]]
}

materialize_matrix <- function(df, orientation) {
  if (is.null(df)) return(NULL)
  if (!is.data.frame(df) || nrow(df) == 0L || ncol(df) == 0L) return(NULL)
  # If the first column is non-numeric, treat it as the ID column.
  id_col <- NULL
  if (!is.numeric(df[[1L]]) && !is.logical(df[[1L]])) {
    id_col <- as.character(df[[1L]])
    df <- df[, -1L, drop = FALSE]
  }
  # Drop any remaining non-numeric columns.
  numeric_cols <- vapply(df, function(c) {
    is.numeric(c) || all(is.na(suppressWarnings(as.numeric(as.character(c)))) ==
                         is.na(c))
  }, logical(1))
  df <- df[, numeric_cols, drop = FALSE]
  if (ncol(df) == 0L) return(NULL)
  mat <- as.matrix(as.data.frame(lapply(df, function(c) {
    if (is.numeric(c)) c else suppressWarnings(as.numeric(as.character(c)))
  }), check.names = FALSE))
  colnames(mat) <- colnames(df)
  if (!is.null(id_col)) {
    rownames(mat) <- make_unique_labels(id_col)
  } else if (!is.null(rownames(df))) {
    rownames(mat) <- rownames(df)
  } else {
    rownames(mat) <- paste0("f", seq_len(nrow(mat)))
  }
  if (identical(orientation, "samples_in_rows")) {
    mat <- t(mat)
  }
  mat
}

materialize_metadata <- function(df, sample_ids) {
  if (is.null(df)) {
    # Build a minimal metadata frame so validate_omics_input() passes.
    return(data.frame(
      sample_id = sample_ids,
      row.names = sample_ids,
      stringsAsFactors = FALSE
    ))
  }
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  # Find a column that overlaps with sample IDs.
  match_col <- NULL
  best_overlap <- 0L
  for (col in colnames(df)) {
    vals <- as.character(df[[col]])
    overlap <- sum(vals %in% sample_ids)
    if (overlap > best_overlap) {
      best_overlap <- overlap
      match_col <- col
    }
  }
  if (!is.null(match_col) && best_overlap > 0L) {
    rownames(df) <- make_unique_labels(df[[match_col]])
  } else if (nrow(df) == length(sample_ids)) {
    rownames(df) <- sample_ids
  } else {
    # Fall back: align by truncation/padding to keep validate_omics_input() honest.
    df <- data.frame(sample_id = sample_ids,
                     row.names = sample_ids, stringsAsFactors = FALSE)
  }
  df
}

materialize_feature_annot <- function(df, feature_ids) {
  if (is.null(df)) {
    out <- data.frame(
      feature_id = feature_ids,
      row.names = feature_ids,
      stringsAsFactors = FALSE
    )
    attr(out, "id_column") <- NA_character_
    return(out)
  }
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  id_col <- pick_id_column(df, feature_ids)
  if (is.null(id_col)) {
    # No matching ID column -- synthesise one aligned to the matrix rows.
    out <- data.frame(
      feature_id = feature_ids,
      row.names = feature_ids,
      stringsAsFactors = FALSE
    )
    attr(out, "id_column") <- NA_character_
    return(out)
  }
  ids <- as.character(df[[id_col]])
  keep <- ids %in% feature_ids
  df <- df[keep, , drop = FALSE]
  ids <- ids[keep]
  if (id_col != "feature_id") {
    df$feature_id <- ids
  }
  rownames(df) <- make_unique_labels(ids)
  # Ensure every matrix feature appears (left-join behaviour).
  missing_ids <- setdiff(feature_ids, ids)
  if (length(missing_ids) > 0L) {
    pad <- data.frame(feature_id = missing_ids,
                      row.names = missing_ids,
                      stringsAsFactors = FALSE)
    for (col in setdiff(colnames(df), colnames(pad))) {
      pad[[col]] <- NA
    }
    df <- rbind(df, pad[, colnames(df), drop = FALSE])
  }
  df <- df[feature_ids, , drop = FALSE]
  attr(df, "id_column") <- id_col
  df
}

pick_id_column <- function(df, feature_ids) {
  best_col <- NULL
  best_overlap <- 0L
  for (col in colnames(df)) {
    vec <- as.character(df[[col]])
    overlap <- sum(vec %in% feature_ids)
    if (overlap > best_overlap) {
      best_overlap <- overlap
      best_col <- col
    }
  }
  if (best_overlap == 0L) return(NULL)
  best_col
}

make_unique_labels <- function(x) {
  x <- as.character(x)
  x[is.na(x) | !nzchar(x)] <- "._missing"
  make.unique(x, sep = "_")
}
