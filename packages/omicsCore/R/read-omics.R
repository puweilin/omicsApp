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
detect_delimiter <- function(path, line = NULL) {
  if (is.null(line)) {
    line <- tryCatch(readLines(path, n = 1L, warn = FALSE),
                     error = function(e) character(0))
  }
  if (length(line) == 0L || !nzchar(line[1])) return(",")
  candidates <- c("\t", ",", ";")
  counts <- vapply(candidates,
                   function(d) lengths(regmatches(line[1], gregexpr(d, line[1], fixed = TRUE))),
                   integer(1L))
  if (max(counts) == 0L) return(",")
  candidates[which.max(counts)]
}

# How many leading lines to skip before the header. featureCounts writes
# its command line first, as "# Program:featureCounts ...", and
# read.table takes whatever comes first as the header -- one column,
# and then "more columns than column names" on the second line. Only
# leading lines are considered, and only ones that begin with "#" and
# carry fewer fields than the line after them: "#" inside data is an
# ordinary character (see comment.char below), and a header whose first
# cell happens to start with "#" has as many fields as its rows.
count_leading_comment_lines <- function(path) {
  head_lines <- tryCatch(readLines(path, n = 25L, warn = FALSE),
                         error = function(e) character(0))
  skip <- 0L
  while (skip + 1L < length(head_lines) && startsWith(head_lines[skip + 1L], "#")) {
    nxt <- head_lines[skip + 2L]
    sep <- detect_delimiter(path, line = nxt)
    n_fields <- function(x) lengths(regmatches(x, gregexpr(sep, x, fixed = TRUE)))
    if (n_fields(head_lines[skip + 1L]) >= n_fields(nxt)) break
    skip <- skip + 1L
  }
  skip
}

read_omics_csv <- function(path, omics_type, assay_type,
                           sheet_roles = NULL, ...) {
  skip <- count_leading_comment_lines(path)
  header_line <- tryCatch(readLines(path, n = skip + 1L, warn = FALSE)[skip + 1L],
                          error = function(e) NA_character_)
  sep <- detect_delimiter(path, line = header_line)

  args <- list(path, header = TRUE, sep = sep, skip = skip,
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
  # Vendor tables carry QC counts and annotation beside the samples; the
  # numeric ones would otherwise become samples. Done on the table, where
  # the column names still say what each column is.
  cols <- select_sample_columns(sheets[[matrix_sheet]])
  for (note in cols$notes) report <- add_import_warning(report, note)
  mat <- materialize_matrix(cols$df, orient)
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

  # Before metadata is matched against the column names, because the
  # names are what changes here.
  picked <- select_measurement_columns(mat)
  mat <- picked$mat
  # Surfaced rather than done quietly: dropping half the columns is a
  # decision made on the user's behalf, and the count they see in the
  # schema review should be one they can reconcile with their file.
  if (!is.null(picked$note)) {
    report <- add_import_warning(report, picked$note)
  }

  meta <- materialize_metadata(if (is.null(metadata_sheet)) NULL
                               else sheets[[metadata_sheet]],
                               sample_ids = colnames(mat))
  # A separate annotation sheet wins; failing that, the annotation
  # columns the matrix sheet itself carried.
  feat_source <- if (!is.null(feature_sheet)) sheets[[feature_sheet]] else cols$annotation
  feat <- materialize_feature_annot(feat_source, feature_ids = rownames(mat))

  # An Ensembl-keyed matrix needs symbols before any pathway database can
  # be matched against it. Done here rather than left to the user,
  # because the alternative -- an empty enrichment result -- is
  # indistinguishable from a real one that found nothing.
  sym <- attach_hgnc_symbols(feat, rownames(mat))
  feat <- sym$feature_df
  if (!is.null(sym$note)) report <- add_import_warning(report, sym$note)

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

# Suffixes a sequencing vendor uses to say what a column measures. Only
# these are recognised: an unknown trailing word is far more likely to be
# part of the sample name than a unit.
MEASURE_SUFFIXES <- c(
  count = "count", counts = "count", readcount = "count",
  fpkm = "fpkm", tpm = "tpm", cpm = "cpm", rpkm = "rpkm"
)

# The last underscore-delimited word only. `_[A-Za-z_]+$` would be
# greedy across separators and read `RD001_Folli_FPKM` as `Folli_FPKM`,
# which matches nothing -- the table then looks unit-free and both units
# survive into one matrix.
MEASURE_SUFFIX_RE <- "_[A-Za-z]+$"

measure_suffix <- function(x) {
  m <- regmatches(x, regexpr(MEASURE_SUFFIX_RE, x))
  if (length(m) == 0L) return(NA_character_)
  unname(MEASURE_SUFFIXES[tolower(sub("^_", "", m))])
}

#' Keep one measurement per sample
#'
#' Vendor reports often put every sample twice, once per unit --
#' `RD001_Folli_FPKM` beside `RD001_Folli_count`. Read whole, that is a
#' matrix where half the columns are on a different scale from the other
#' half, and nothing downstream can tell: PCA, clustering and
#' differential testing all run and all produce nonsense.
#'
#' Counts win when present, because they are what the differential
#' methods need; FPKM cannot be given to DESeq2 at all. The suffix is
#' then dropped, which is also what makes the two report layouts agree --
#' the counts-only file names the same sample `RD001_Folli`, so metadata
#' written for one matches the other.
#'
#' @param mat Numeric matrix, samples in columns.
#' @return `list(mat, note)`; `note` is `NULL` when nothing was changed.
#' @keywords internal
#' @noRd
select_measurement_columns <- function(mat) {
  unchanged <- list(mat = mat, note = NULL)
  nms <- colnames(mat)
  if (is.null(nms) || ncol(mat) == 0L) return(unchanged)

  suffix <- vapply(nms, measure_suffix, character(1L), USE.NAMES = FALSE)
  # Any column without a recognised unit means this is not a
  # units-in-the-header table, and guessing on a partial match would
  # silently drop real samples.
  if (anyNA(suffix)) return(unchanged)

  groups <- unique(suffix)
  keep_unit <- if ("count" %in% groups) "count" else groups[[1L]]
  keep <- suffix == keep_unit

  stripped <- sub(MEASURE_SUFFIX_RE, "", nms[keep])
  # Only rename when it stays unambiguous.
  if (anyDuplicated(stripped) == 0L) {
    out <- mat[, keep, drop = FALSE]
    colnames(out) <- stripped
  } else {
    out <- mat[, keep, drop = FALSE]
  }

  if (length(groups) == 1L && identical(colnames(out), nms)) return(unchanged)

  note <- if (length(groups) > 1L) {
    sprintf(
      "Columns carried units (%s); kept the %d %s column(s) and dropped %d.",
      paste(sort(groups), collapse = ", "), sum(keep), keep_unit, sum(!keep))
  } else {
    sprintf("Dropped the '%s' suffix from %d sample name(s).",
            keep_unit, ncol(out))
  }
  list(mat = out, note = note)
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
  # Text columns are read as numbers where every cell is a number or a
  # way of writing "missing". Drop the rest.
  df[] <- lapply(df, function(c) if (is.numeric(c)) c else clean_numeric_text(c))
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

# The ways a cell says "no value". Spectronaut writes "Filtered", Excel
# exports "#N/A", people write "n.d."; read as text, any one of them
# turned a whole sample column non-numeric, and the column was dropped
# without a word -- a missing sample rather than a missing value.
MISSING_TOKENS <- c("", "na", "n/a", "#n/a", "#na", "nan", "null", "none",
                    "filtered", "n.d.", "nd", "-", "--", "?")

# Numbers written for people: thousands separators, and the tokens above.
# Only a column where every non-missing cell reads as a number after the
# cleanup is treated as numeric, so a genuinely textual column is left
# alone.
clean_numeric_text <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (!is.character(x)) return(x)
  y <- trimws(x)
  y[tolower(y) %in% MISSING_TOKENS] <- NA_character_
  present <- y[!is.na(y)]
  if (length(present) == 0L) return(y)
  grouped <- "^[-+]?[0-9]{1,3}(,[0-9]{3})+(\\.[0-9]+)?$"
  plain <- "^[-+]?([0-9]+\\.?[0-9]*|\\.[0-9]+)([eE][-+]?[0-9]+)?$"
  if (any(grepl(grouped, present)) && all(grepl(grouped, present) | grepl(plain, present))) {
    y <- gsub(",", "", y, fixed = TRUE)
  }
  y
}

# Columns a quantification tool writes beside the samples, by the family
# they belong to. MaxQuant names every per-sample column "<measure>
# <sample>" and adds totals and QC counts under bare names; the sample
# block is the one family, and the rest is not data. Listed in order of
# preference: LFQ is what MaxQuant users analyse, raw Intensity is what
# they fall back to.
SAMPLE_COLUMN_FAMILIES <- c("LFQ intensity ", "iBAQ ", "Intensity ",
                            "Reporter intensity corrected ", "Reporter intensity ")

# featureCounts' fixed header. Start/End/Strand are ";"-joined text and
# fall away on their own; Length is numeric and became a seventh sample.
FEATURECOUNTS_ANNOTATION <- c("Chr", "Start", "End", "Strand", "Length")

# Vendor suffixes on sample columns: Spectronaut's "[1] S01.raw.PG.Quantity",
# DIA-NN's raw-file paths. Stripped so the sample is called what the
# metadata sheet calls it. Only these; a prefix shared by every sample
# ("Cheek_01", "Cheek_02") is part of the name and stays.
strip_vendor_decoration <- function(x) {
  x <- sub("^\\[[0-9]+\\]\\s*", "", x)
  x <- sub("^.*[\\\\/]", "", x)
  x <- sub("\\.(PG|PEP|EG|FG)\\.[A-Za-z]+$", "", x)
  x <- sub("\\.(raw|d|wiff|mzML|mzXML|bam|sam|cram)$", "", x, ignore.case = TRUE)
  x
}

#' Keep the columns that are samples
#'
#' @param df The matrix sheet as read, first column the feature ids.
#' @return `list(df, notes)`; `notes` is `character(0)` when nothing
#'   changed.
#' @keywords internal
#' @noRd
select_sample_columns <- function(df) {
  unchanged <- list(df = df, notes = character(0), annotation = NULL)
  if (!is.data.frame(df) || ncol(df) < 3L) return(unchanged)
  notes <- character(0)
  nms <- colnames(df)

  # The text columns beside the samples are the feature annotation a
  # single-table export carries -- "Gene names" in MaxQuant, "PG.Genes"
  # in Spectronaut, "Genes" in DIA-NN. Kept aside for the case where the
  # file has no separate annotation sheet, so the symbol column is
  # picked up by the same rules that read one.
  is_text <- vapply(df, function(c) !is.numeric(c) && !is.logical(c), logical(1))
  is_text[1L] <- TRUE
  annotation <- if (sum(is_text) > 1L) df[, is_text, drop = FALSE] else NULL

  if (all(c("Geneid", FEATURECOUNTS_ANNOTATION) %in% nms)) {
    df <- df[, setdiff(nms, FEATURECOUNTS_ANNOTATION), drop = FALSE]
    nms <- colnames(df)
    notes <- c(notes, sprintf(
      "Dropped featureCounts' annotation columns (%s); they are not samples.",
      paste(FEATURECOUNTS_ANNOTATION, collapse = ", ")))
  }

  id_col <- nms[[1L]]
  for (family in SAMPLE_COLUMN_FAMILIES) {
    members <- nms[startsWith(nms, family)]
    if (length(members) < 2L) next
    others <- setdiff(nms, c(id_col, members))
    keep <- df[, c(id_col, members), drop = FALSE]
    colnames(keep) <- c(id_col, sub(family, "", members, fixed = TRUE))
    notes <- c(notes, sprintf(
      "Kept the %d '%s' columns as samples and dropped %d other column(s) (%s).",
      length(members), trimws(family), length(others),
      paste(utils::head(others, 6L), collapse = ", ")))
    df <- keep
    nms <- colnames(df)
    break
  }

  decorated <- nms[-1L]
  plain <- strip_vendor_decoration(decorated)
  if (!identical(plain, decorated) && !anyDuplicated(plain) && all(nzchar(plain))) {
    colnames(df) <- c(nms[[1L]], plain)
    notes <- c(notes, sprintf(
      "Shortened %d sample name(s) by their vendor decoration, e.g. '%s' to '%s'.",
      sum(plain != decorated), decorated[plain != decorated][[1L]],
      plain[plain != decorated][[1L]]))
  }

  list(df = df, notes = notes, annotation = annotation)
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
