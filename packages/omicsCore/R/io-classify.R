# Heuristics used by `read_omics()` to label each tab/file in an Excel
# workbook (or a single data frame from a CSV/RDS) as one of:
#   * "matrix"        - numeric expression matrix, samples in one axis
#   * "metadata"      - sample metadata table (groups, batches, age, ...)
#   * "feature_annot" - feature metadata (gene symbols, descriptions, ...)
#   * "unknown"       - could not place this with reasonable confidence
#
# Each classifier returns a list(role, confidence, orientation, notes)
# so the caller can build the `ImportReport` sheets row directly.

# ID-column regular expressions used to detect feature-annotation tables
# and to choose an ID column inside one. The patterns are intentionally
# loose; the classifier only needs to know that *some* column matches a
# biological identifier convention.
ID_COLUMN_PATTERNS <- list(
  uniprot = "^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$",
  ensembl_gene = "^ENS[A-Z]*G[0-9]{6,}",
  ensembl_tx   = "^ENS[A-Z]*T[0-9]{6,}",
  refseq       = "^N[MR]_[0-9]+",
  hgnc_symbol  = "^[A-Z][A-Z0-9-]{1,14}$"
)

# Column-name hints that drive the metadata vs feature_annot decision.
METADATA_NAME_HINTS <- c(
  "sample", "sample_id", "subject", "donor", "donor_id", "patient",
  "group", "condition", "batch", "site", "tissue", "age", "sex",
  "gender", "time", "timepoint", "replicate"
)

FEATURE_NAME_HINTS <- c(
  "feature", "feature_id", "gene", "gene_id", "gene_symbol",
  "symbol", "uniprot", "protein", "protein_id", "ensembl", "ensg",
  "transcript", "description", "annotation"
)

#' Detect biological ID columns in a data.frame
#'
#' Scans column names and the first non-missing value in each column;
#' tags any column that matches a UniProt / Ensembl / RefSeq / HGNC
#' regex.
#'
#' @param df A `data.frame`.
#' @param max_check Number of rows to sample per column when probing
#'   values. Defaults to 200 to keep large workbooks responsive.
#'
#' @return A `data.frame` with columns `column`, `pattern`, `match_rate`.
#'   Empty data frame if no matches were found.
#' @export
#' @family io
detect_id_columns <- function(df, max_check = 200L) {
  out <- data.frame(
    column = character(0),
    pattern = character(0),
    match_rate = numeric(0),
    stringsAsFactors = FALSE
  )
  if (!is.data.frame(df) || ncol(df) == 0L) return(out)
  for (col in colnames(df)) {
    vec <- df[[col]]
    if (is.factor(vec)) vec <- as.character(vec)
    if (!is.character(vec)) next
    vec <- vec[!is.na(vec) & nzchar(vec)]
    if (length(vec) == 0L) next
    probe <- utils::head(vec, max_check)
    for (pat_name in names(ID_COLUMN_PATTERNS)) {
      rate <- mean(grepl(ID_COLUMN_PATTERNS[[pat_name]], probe, perl = TRUE))
      if (rate >= 0.6) {
        out <- rbind(out, data.frame(
          column = col,
          pattern = pat_name,
          match_rate = rate,
          stringsAsFactors = FALSE
        ))
        break
      }
    }
  }
  out
}

#' Detect the orientation of a candidate expression matrix
#'
#' Two-way guess: are samples in columns (rows = features) or in rows
#' (rows = samples)? The heuristic favours the axis whose labels look
#' more like biological identifiers (long, alphanumeric, no spaces) and
#' tiebreaks on the typical "many features, few samples" shape.
#'
#' @param df A `data.frame` whose first column may be an ID column.
#'
#' @return A list with `orientation` (one of `"features_in_rows"`,
#'   `"samples_in_rows"`, `"ambiguous"`), `confidence`, and `notes`.
#' @export
#' @family io
detect_orientation <- function(df) {
  if (!is.data.frame(df)) {
    return(list(orientation = "ambiguous", confidence = 0,
                notes = "input is not a data.frame"))
  }
  # Strip a leading ID column if present (by content or, for numeric
  # identifiers such as Entrez ids, by name -- see first_column_is_id()).
  body <- df
  first_col_is_id <- ncol(df) > 1L && first_column_is_id(df)
  if (first_col_is_id) {
    row_labels <- as.character(df[[1L]])
    body <- df[, -1L, drop = FALSE]
  } else {
    row_labels <- rownames(df)
  }
  # Trimmed before judging. A cell " ENSG00000000001 " is not an
  # Ensembl id to the regex, and with no feature-like rows the sample
  # names "S01".."S06" -- which pass as gene symbols -- decided the
  # orientation, and the matrix came back transposed.
  row_labels <- trimws(row_labels)
  col_labels <- trimws(colnames(body))
  rows_look_like_features <- looks_like_feature_labels(row_labels)
  cols_look_like_features <- looks_like_feature_labels(col_labels)

  shape_hint <- if (nrow(body) > ncol(body) * 2L) "features_in_rows"
                else if (ncol(body) > nrow(body) * 2L) "samples_in_rows"
                else NA_character_

  if (rows_look_like_features && !cols_look_like_features) {
    return(list(orientation = "features_in_rows", confidence = 0.9,
                notes = "row labels look like feature IDs"))
  }
  # Column labels alone are weaker evidence than row labels: "S01" passes
  # for a gene symbol, so a matrix with no id column and a handful of
  # samples used to come back transposed. Unless the columns carry a
  # strong identifier pattern (UniProt, Ensembl, RefSeq), the shape has
  # to agree -- a transposed matrix has far more columns than rows.
  if (cols_look_like_features && !rows_look_like_features &&
      (strongly_feature_like(col_labels) || ncol(body) >= nrow(body))) {
    return(list(orientation = "samples_in_rows", confidence = 0.9,
                notes = "column labels look like feature IDs"))
  }
  if (!is.na(shape_hint)) {
    return(list(orientation = shape_hint, confidence = 0.55,
                notes = "decided by row vs column count"))
  }
  list(orientation = "features_in_rows", confidence = 0.4,
       notes = "fallback default (features in rows)")
}

#' Classify a sheet/data.frame's role in an omics import
#'
#' Walks a small decision tree to label one sheet as `"matrix"`,
#' `"metadata"`, `"feature_annot"`, or `"unknown"`. Returns the
#' confidence in `[0, 1]` so callers can present alternatives.
#'
#' @param df A `data.frame`.
#' @param name Optional sheet name (used as a tiebreaker hint).
#'
#' @return A list with `role`, `confidence`, `orientation`, `notes`.
#' @export
#' @family io
classify_sheet_role <- function(df, name = NA_character_) {
  if (!is.data.frame(df) || nrow(df) == 0L || ncol(df) == 0L) {
    return(list(role = "unknown", confidence = 0,
                orientation = NA_character_,
                notes = "empty or non-tabular"))
  }
  numeric_fraction <- compute_numeric_fraction(df)
  id_hits <- detect_id_columns(df)
  name_lc <- tolower(name %||% "")
  has_meta_hint <- any(tolower(colnames(df)) %in% METADATA_NAME_HINTS)
  has_feat_hint <- any(tolower(colnames(df)) %in% FEATURE_NAME_HINTS)
  name_meta <- grepl("meta|sample|pheno|design|cohort", name_lc)
  name_feat <- grepl("feat|gene|protein|annot|symbol", name_lc)
  name_matrix <- grepl("expr|intens|count|matrix|abund|fpkm|tpm", name_lc)

  # Matrix: a large fraction of cells is numeric, plus optional name hint.
  if (numeric_fraction >= 0.8) {
    orient <- detect_orientation(df)
    conf <- 0.7 + 0.2 * (numeric_fraction - 0.8) / 0.2
    if (name_matrix) conf <- min(conf + 0.1, 0.99)
    return(list(role = "matrix",
                confidence = round(conf, 3),
                orientation = orient$orientation,
                notes = sprintf("%.0f%% numeric cells", 100 * numeric_fraction)))
  }

  # A vendor report: a block of numeric sample columns beside a handful
  # of annotation columns. MaxQuant's proteinGroups.txt carries five text
  # columns and a dozen QC counts before the first LFQ intensity, a
  # Spectronaut report four text columns before its quantities -- and at
  # six samples neither comes near 80% numeric, so both were "no sheet
  # looked like an expression matrix" and the import stopped there. The
  # block is what makes it a matrix: several numeric columns, at least
  # half the sheet, behind a first column of identifiers, on a table
  # with more rows than a sample sheet has, and none of the column names
  # that mark a sample sheet. A metadata sheet with three numeric
  # columns and no recognisable name is the case this could get wrong,
  # which is why the confidence is low enough for the review step to
  # show it as a guess.
  n_numeric <- round(numeric_fraction * ncol(df))
  first_col_ids <- first_column_is_id(df) &&
    (is.numeric(df[[1L]]) || looks_like_feature_labels(df[[1L]]))
  if (n_numeric >= 3L && numeric_fraction >= 0.5 && first_col_ids &&
      !has_meta_hint && nrow(df) >= 20L) {
    orient <- detect_orientation(df)
    conf <- if (name_matrix) 0.7 else 0.6
    return(list(role = "matrix",
                confidence = conf,
                orientation = orient$orientation,
                notes = sprintf("%d numeric column(s) beside %d annotation column(s)",
                                n_numeric, ncol(df) - n_numeric)))
  }

  # Metadata is checked before feature_annot so that a column named
  # `sample_id` (whose values may also match the HGNC regex) doesn't get
  # mistaken for a gene-symbol table.
  meta_hint_strength <- sum(tolower(colnames(df)) %in% METADATA_NAME_HINTS)
  feat_hint_strength <- sum(tolower(colnames(df)) %in% FEATURE_NAME_HINTS)
  if (meta_hint_strength > 0L || name_meta) {
    if (meta_hint_strength >= feat_hint_strength || name_meta) {
      conf <- 0.55 + 0.07 * meta_hint_strength
      if (numeric_fraction < 0.4) conf <- conf + 0.1
      if (name_meta) conf <- conf + 0.1
      return(list(role = "metadata",
                  confidence = round(min(conf, 0.95), 3),
                  orientation = NA_character_,
                  notes = paste0("metadata column hints (",
                                 meta_hint_strength, ")")))
    }
  }

  # Feature annotation: looks like a wide-ish table whose first column or
  # one of its named columns matches an ID regex, and that has very few
  # numeric columns.
  if (nrow(id_hits) > 0L && numeric_fraction < 0.4) {
    conf <- 0.6 + 0.05 * nrow(id_hits)
    if (has_feat_hint || name_feat) conf <- conf + 0.15
    return(list(role = "feature_annot",
                confidence = round(min(max(conf, 0), 0.99), 3),
                orientation = NA_character_,
                notes = sprintf("ID column(s): %s",
                                paste(id_hits$column, collapse = ", "))))
  }

  if (has_feat_hint || name_feat) {
    return(list(role = "feature_annot",
                confidence = 0.55,
                orientation = NA_character_,
                notes = "feature column hints in header"))
  }

  list(role = "unknown",
       confidence = round(0.3 + 0.2 * numeric_fraction, 3),
       orientation = NA_character_,
       notes = "no strong signal")
}

# ---- internal helpers --------------------------------------------------

compute_numeric_fraction <- function(df) {
  if (ncol(df) == 0L) return(0)
  # Probe up to 1000 cells per column without materialising the whole table.
  sample_rows <- min(nrow(df), 500L)
  if (sample_rows == 0L) return(0)
  rows_idx <- seq_len(sample_rows)
  numeric_cols <- vapply(df, function(col) {
    if (is.numeric(col)) return(TRUE)
    if (is.logical(col)) return(FALSE)
    sub <- col[rows_idx]
    # The same cleanup the reader applies, so a sample column that says
    # "Filtered" in a few cells counts as the numbers it mostly is.
    sub <- clean_numeric_text(as.character(sub))
    sub <- sub[!is.na(sub) & nzchar(sub)]
    if (length(sub) == 0L) return(FALSE)
    coerced <- suppressWarnings(as.numeric(sub))
    mean(!is.na(coerced)) >= 0.9
  }, logical(1))
  mean(numeric_cols)
}

# UniProt, Ensembl or RefSeq accessions, which nothing else looks like.
strongly_feature_like <- function(labels) {
  labels <- trimws(as.character(labels))
  labels <- labels[!is.na(labels) & nzchar(labels)]
  if (length(labels) < 4L) return(FALSE)
  strong <- ID_COLUMN_PATTERNS[c("uniprot", "ensembl_gene", "ensembl_tx", "refseq")]
  probe <- utils::head(labels, 200L)
  hits <- vapply(probe, function(s) {
    any(vapply(strong, function(pat) grepl(pat, s, perl = TRUE), logical(1)))
  }, logical(1))
  mean(hits) >= 0.5
}

looks_like_feature_labels <- function(labels) {
  labels <- trimws(as.character(labels))
  labels <- labels[!is.na(labels) & nzchar(labels)]
  if (length(labels) < 4L) return(FALSE)
  probe <- utils::head(labels, 200L)
  hits <- vapply(probe, function(s) {
    any(vapply(ID_COLUMN_PATTERNS, function(pat) {
      grepl(pat, s, perl = TRUE)
    }, logical(1)))
  }, logical(1))
  # Long alphanumeric tokens without spaces are weakly feature-like, too.
  longish <- nchar(probe) >= 4L & !grepl("[[:space:]]", probe)
  mean(hits) >= 0.5 || mean(hits | longish) >= 0.7
}
