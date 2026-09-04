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
#' @param sheet_roles Optional named character vector overriding what the
#'   classifier decided, as `c("<sheet name>" = "<role>")` with roles drawn
#'   from [IMPORT_REPORT_ROLES]. Sheets left out keep their inferred role.
#'
#'   The classifier is a heuristic, and a confident wrong answer is worse than
#'   an unconfident one: a metadata sheet read as the matrix produces an
#'   `omics_input` that analyses cleanly and means nothing. This is how a
#'   caller lets the user correct it without re-implementing the reader.
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
  sheet_roles = NULL,
  ...
) {
  type <- match.arg(type)
  sheet_roles <- validate_sheet_roles(sheet_roles)
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
                             assay_type = assay_type,
                             sheet_roles = sheet_roles, ...),
    csv = read_omics_csv(path, omics_type = omics_type,
                         assay_type = assay_type,
                         sheet_roles = sheet_roles, ...),
    rds = read_omics_rds(path, omics_type = omics_type,
                         assay_type = assay_type,
                         sheet_roles = sheet_roles, ...),
    stop("Unsupported type: ", type)
  )
}

# Roles the caller may assign, checked before any file is opened so a typo
# surfaces as a message about the argument rather than as a mysteriously
# unchanged import.
validate_sheet_roles <- function(sheet_roles) {
  if (is.null(sheet_roles)) return(NULL)
  if (!is.character(sheet_roles) || is.null(names(sheet_roles)) ||
      any(!nzchar(names(sheet_roles)))) {
    stop("`sheet_roles` must be a named character vector, ",
         "e.g. c(Sheet1 = \"matrix\").", call. = FALSE)
  }
  bad <- setdiff(sheet_roles, IMPORT_REPORT_ROLES)
  if (length(bad) > 0L) {
    stop("Unknown role(s) in `sheet_roles`: ",
         paste(sprintf("'%s'", bad), collapse = ", "),
         ". Valid roles are ",
         paste(sprintf("'%s'", IMPORT_REPORT_ROLES), collapse = ", "), ".",
         call. = FALSE)
  }
  sheet_roles
}

# Apply the caller's assignment on top of the classifier's. Overridden sheets
# get confidence 1 so pick_best_sheet() prefers them over anything the
# classifier was merely sure about, and the note records that a human chose.
apply_sheet_roles <- function(sheet_table, sheet_roles) {
  if (is.null(sheet_roles) || nrow(sheet_table) == 0L) return(sheet_table)
  for (nm in names(sheet_roles)) {
    idx <- which(sheet_table$name == nm)
    if (length(idx) == 0L) next
    sheet_table$role[idx] <- unname(sheet_roles[[nm]])
    sheet_table$confidence[idx] <- 1
    sheet_table$notes[idx] <- "role set by user"
  }
  sheet_table
}

# ---- internal: per-type readers ----------------------------------------

# Both Excel formats announce themselves in their first bytes: .xlsx is
# a zip, .xls is an OLE2 compound file. Anything else with those
# extensions is text pretending.
is_excel_file <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  magic <- readBin(con, "raw", n = 8L)
  if (length(magic) < 8L) return(FALSE)
  identical(magic[1:4], as.raw(c(0x50, 0x4b, 0x03, 0x04))) ||
    identical(magic, as.raw(c(0xd0, 0xcf, 0x11, 0xe0, 0xa1, 0xb1, 0x1a, 0xe1)))
}

detect_file_type <- function(path) {
  ext <- tolower(tools::file_ext(path))

  # Trust the bytes over the name. RNA-seq pipelines routinely write
  # tab-separated text and call it .xls -- a 170 MB expression table
  # with that extension is normal output, not a mistake -- and handing
  # it to readxl fails with a message about the workbook rather than
  # about the format, which sends you looking at the wrong thing.
  if (ext %in% c("xls", "xlsx") && file.exists(path) && !is_excel_file(path)) {
    return("csv")
  }

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

read_omics_excel <- function(path, omics_type, assay_type,
                             sheet_roles = NULL, ...) {
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
  sheet_table <- apply_sheet_roles(do.call(rbind, rows), sheet_roles)
  build_input_from_sheets(sheets, sheet_table, source = path,
                          omics_type = omics_type, assay_type = assay_type)
}

# From the header line, not the extension: a .xls that turned out to be
# text is tab-separated more often than not, and guessing comma there
# yields one column holding the whole row -- which classify_sheet_role()
# then reports as an unrecognisable sheet rather than as a parse
# failure. Counted on the header alone, where no field is quoted.
detect_delimiter <- function(path) {
  line <- tryCatch(readLines(path, n = 1L, warn = FALSE),
                   error = function(e) character(0))
  if (length(line) == 0L || !nzchar(line[1])) return(",")
  candidates <- c("\t", ",", ";")
  counts <- vapply(candidates,
                   function(d) lengths(regmatches(line[1], gregexpr(d, line[1], fixed = TRUE))),
                   integer(1L))
  if (max(counts) == 0L) return(",")
  candidates[which.max(counts)]
}

read_omics_csv <- function(path, omics_type, assay_type,
                           sheet_roles = NULL, ...) {
  sep <- detect_delimiter(path)

  args <- list(path, header = TRUE, sep = sep,
               check.names = FALSE, stringsAsFactors = FALSE, ...)

  # Quoting rules differ by ecosystem, so they follow the delimiter.
  # Tab-separated files from bioinformatics pipelines do not quote
  # anything, and their annotation columns are full of characters that
  # read.table treats as quotes -- a gene called "5'-nucleotidase" opens
  # a quote that stays open until the next apostrophe, swallowing whole
  # rows. The failure is `line N did not have K elements`, thousands of
  # lines after the one that caused it.
  #
  # Comma-separated files are usually written by a spreadsheet, where a
  # field containing a comma is quoted and the quotes must be honoured.
  if (is.null(args$quote)) {
    args$quote <- if (identical(sep, "\t")) "" else "\""
  }
  # Never: `#` is a legal character in a gene description and truncating
  # the line at one loses columns silently.
  if (is.null(args$comment.char)) args$comment.char <- ""

  df <- do.call(utils::read.table, args)
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
  sheet_table <- apply_sheet_roles(sheet_table, sheet_roles)
  build_input_from_sheets(sheets, sheet_table, source = path,
                          omics_type = omics_type, assay_type = assay_type)
}

read_omics_rds <- function(path, omics_type, assay_type,
                           sheet_roles = NULL, ...) {
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
  sheet_table <- apply_sheet_roles(sheet_table, sheet_roles)
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
  sym_col <- pick_symbol_column(df, exclude = id_col)
  if (!is.null(sym_col)) {
    df$feature_symbol <- first_gene_symbol(df[[sym_col]])
  }
  attr(df, "symbol_column") <- sym_col %||% NA_character_
  df
}

# Column names the instrument vendors use for the gene symbol.
# Matched after lowercasing and stripping everything that is not a
# letter, so "Gene names", "Gene.Name", "PG.Genes" and "gene_symbol" all
# collapse onto an entry here.
SYMBOL_COLUMN_NAMES <- c(
  "genesymbol", "genesymbols", "symbol", "symbols",
  "gene", "genes", "genename", "genenames",
  "pggenes",        # Spectronaut
  "hgncsymbol", "hgnc"
)

# Which column holds the gene symbol, or NULL.
#
# Matched by name rather than by content, unlike pick_id_column(), which
# can check its guess against the matrix row names. There is nothing to
# check a symbol against: any short uppercase string looks like a gene,
# and guessing wrong here would relabel every feature in every plot and
# every enrichment with something that is not its name. A vendor's column
# heading is a stated fact; the shape of the values is not.
pick_symbol_column <- function(df, exclude = NULL) {
  if (is.null(df) || !ncol(df)) return(NULL)
  if ("feature_symbol" %in% colnames(df)) return(NULL)  # already named
  norm <- tolower(gsub("[^A-Za-z]", "", colnames(df)))
  hit <- which(norm %in% SYMBOL_COLUMN_NAMES &
                 !colnames(df) %in% c(exclude, "feature_id"))
  if (!length(hit)) return(NULL)
  # Earliest match wins; the vendor list is ordered most specific first
  # but a workbook can carry more than one, and the leftmost is the one
  # a reader would have taken.
  colnames(df)[hit[1L]]
}

# MaxQuant and friends pack several genes into one cell, separated by
# ";" or ",". The first is the leading identification; keeping the whole
# string would give an enrichment query that matches nothing.
first_gene_symbol <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("[;,].*$", "", x)
  x[!nzchar(x)] <- NA_character_
  x
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
