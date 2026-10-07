# SummarizedExperiment and DESeqDataSet objects, from an .rds upload or
# handed over in R.
#
# What is taken: one assay as the matrix (the read counts when there are
# any, otherwise the first assay, and the report says which), colData as
# the sample information, rowData as the feature annotation. Nothing
# else -- not the design formula, the fitted dispersions, the metadata
# list -- because the analyses here start from the matrix and re-estimate
# the rest.
#
# HOW AN UPLOADED OBJECT IS READ SAFELY. An .rds from anywhere is checked
# by check_project_structure() (persistence.R): data only, of classes the
# app writes. A SummarizedExperiment fails that check by design -- it is
# an S4 object, of classes the app never writes -- and a DESeqDataSet
# carries code outright: a `dispersionFunction` slot holding a function,
# and a `design` formula holding an environment. Loosening the check to
# admit them would admit those too. Instead the object is never used as
# what it claims to be:
#   * only the three exact classes in SE_CLASSES are recognised, by the
#     class attribute, without asking the methods package;
#   * the slots read are fetched with attr(), which runs no method and no
#     code the file supplied, and each is checked to be the plain shape it
#     should be (a list, a character vector, a matrix of numbers) before
#     anything else is done with it;
#   * of the sample and feature tables, only columns that are plain
#     vectors or factors are copied out; anything else (a nested table, a
#     list column, an S4 vector) is left behind and named in the report;
#   * what was copied out is base R -- a matrix and two data frames -- and
#     goes through check_project_structure() before it is used, the same
#     gate an uploaded table goes through.
# The slots not read (design, dispersionFunction, metadata, the rest of
# rowRanges) are never touched, so whatever they hold is never run; it
# is dropped with the object. Deserialising itself is readRDS()'s, and
# was already done for any .rds before this point.

SE_CLASSES <- c("SummarizedExperiment", "RangedSummarizedExperiment", "DESeqDataSet")

# Classes a column copied out of a DataFrame may carry.
SE_COLUMN_CLASSES <- c("factor", "ordered", "Date", "POSIXct", "POSIXt", "difftime", "AsIs")

# The class an object says it has, read without dispatch.
s4_class_name <- function(x) {
  cls <- attr(x, "class", exact = TRUE)
  if (is.character(cls) && length(cls) >= 1L) as.character(cls[[1L]]) else NA_character_
}

# A slot, read as the attribute it is stored in. An empty slot is stored
# as the symbol `\001NULL\001`.
se_slot <- function(obj, name) {
  v <- attr(obj, name, exact = TRUE)
  if (is.name(v) && identical(as.character(v), "\001NULL\001")) NULL else v
}

se_refuse <- function(message) {
  stop(structure(class = c("omics_se_error", "error", "condition"),
                 list(message = message, call = NULL)))
}

is_plain_list <- function(x) {
  typeof(x) == "list" && !isS4(x) && is.null(attr(x, "class", exact = TRUE))
}

#' Is this a SummarizedExperiment the reader handles?
#'
#' @param x Any object.
#' @param trusted `TRUE` for an object handed over in R, where a subclass
#'   (`DESeqTransform`, `SingleCellExperiment`) is recognised through the
#'   methods package; `FALSE` for one read from an uploaded file, where
#'   only the exact classes in `SE_CLASSES` are.
#' @keywords internal
#' @noRd
is_se_object <- function(x, trusted = FALSE) {
  if (!isS4(x)) return(FALSE)
  if (s4_class_name(x) %in% SE_CLASSES) return(TRUE)
  isTRUE(trusted) && is_installed("SummarizedExperiment") &&
    methods::is(x, "SummarizedExperiment")
}

# ---- copying the parts out --------------------------------------------------

# A DataFrame (S4Vectors) as a base data frame of its plain columns.
# Returns list(df, skipped), `df` NULL when the slot is empty.
se_data_frame <- function(df, what) {
  if (is.null(df)) return(list(df = NULL, skipped = character(0)))
  cls <- s4_class_name(df)
  if (!isS4(df) || !cls %in% c("DFrame", "DataFrame")) {
    se_refuse(sprintf("Its %s is stored as '%s', not as the table this reader knows.",
                      what, cls))
  }
  n <- se_slot(df, "nrows")
  rn <- se_slot(df, "rownames")
  ld <- se_slot(df, "listData")
  if (!(typeof(n) %in% c("integer", "double")) || length(n) != 1L || is.na(n) || n < 0) {
    se_refuse(sprintf("Its %s does not say how many rows it has.", what))
  }
  n <- as.integer(n)
  if (!is.null(rn) && (typeof(rn) != "character" || isS4(rn) || length(rn) != n)) {
    se_refuse(sprintf("The row names of its %s are not a list of names.", what))
  }
  if (!is.null(ld) && !is_plain_list(ld)) {
    se_refuse(sprintf("The columns of its %s are not stored as a list.", what))
  }
  nms <- attr(ld, "names", exact = TRUE)
  if (length(ld) && (typeof(nms) != "character" || length(nms) != length(ld))) {
    se_refuse(sprintf("The columns of its %s have no names.", what))
  }
  cols <- list()
  skipped <- character(0)
  for (i in seq_along(ld)) {
    col <- ld[[i]]
    cls <- attr(col, "class", exact = TRUE)
    plain <- !isS4(col) &&
      typeof(col) %in% c("logical", "integer", "double", "character", "complex") &&
      is.null(attr(col, "dim", exact = TRUE)) &&
      (is.null(cls) || (is.character(cls) && all(cls %in% SE_COLUMN_CLASSES))) &&
      length(col) == n
    if (plain) cols[[nms[[i]]]] <- col else skipped <- c(skipped, nms[[i]])
  }
  # The kept columns are walked before anything is built from them: a
  # factor's levels, a date's class, any attribute could hold code, and
  # as.data.frame() would be the first to touch it.
  check_project_structure(cols)
  out <- if (length(cols)) {
    as.data.frame(cols, stringsAsFactors = FALSE, optional = TRUE)
  } else {
    as.data.frame(matrix(nrow = n, ncol = 0L))
  }
  if (length(cols)) names(out) <- make.unique(names(cols))
  rownames(out) <- if (!is.null(rn)) make_unique_labels(rn) else NULL
  attr(out, "source_rownames") <- rn
  list(df = out, skipped = skipped)
}

# One assay as a base numeric matrix: an ordinary matrix as it is, a
# sparse dgCMatrix (the Matrix package's) spelled out from its slots.
se_assay_matrix <- function(m, name,
                            max_mb = getOption("omicsCore.max_unpacked_mb", 2048)) {
  limit <- max_mb * 1024^2
  dense_ok <- function(d) {
    if (prod(as.numeric(d)) * 8 > limit) {
      se_refuse(sprintf(paste(
        "Its assay '%s' holds %s values, more than the %s this server reads.",
        "Save fewer features or samples and upload that."),
        name, format(prod(as.numeric(d)), big.mark = ","), format_mb(limit)))
    }
  }
  good_dim <- function(d) {
    typeof(d) %in% c("integer", "double") && length(d) == 2L && !anyNA(d) && all(d >= 0)
  }
  good_dimnames <- function(dn) {
    is.null(dn) || (is_plain_list(dn) && length(dn) == 2L &&
                      all(vapply(dn, function(v) is.null(v) || (typeof(v) == "character" && !isS4(v)),
                                 logical(1))))
  }
  cls <- s4_class_name(m)
  if (!isS4(m) && typeof(m) %in% c("double", "integer", "logical") &&
      (is.null(attr(m, "class", exact = TRUE)) || all(attr(m, "class", exact = TRUE) %in% c("matrix", "array")))) {
    d <- attr(m, "dim", exact = TRUE)
    dn <- attr(m, "dimnames", exact = TRUE)
    if (!good_dim(d) || !good_dimnames(dn)) {
      se_refuse(sprintf("Its assay '%s' is not a two-dimensional table.", name))
    }
    out <- matrix(as.numeric(m), nrow = d[[1L]], ncol = d[[2L]])
    if (!is.null(dn)) dimnames(out) <- dn
    return(out)
  }
  if (isS4(m) && identical(cls, "dgCMatrix")) {
    d <- se_slot(m, "Dim")
    i <- se_slot(m, "i")
    p <- se_slot(m, "p")
    x <- se_slot(m, "x")
    dn <- se_slot(m, "Dimnames")
    ok <- good_dim(d) && typeof(i) == "integer" && typeof(p) == "integer" &&
      typeof(x) == "double" && length(i) == length(x) && length(p) == d[[2L]] + 1L &&
      !anyNA(p) && !anyNA(i) && p[[1L]] == 0L && !is.unsorted(p) &&
      p[[length(p)]] == length(x) && all(i >= 0L & i < d[[1L]]) &&
      good_dimnames(if (is.list(dn) && all(vapply(dn, is.null, logical(1)))) NULL else dn)
    if (!ok) se_refuse(sprintf("Its sparse assay '%s' is not laid out as a sparse matrix is.", name))
    dense_ok(d)
    out <- matrix(0, nrow = d[[1L]], ncol = d[[2L]])
    if (length(x)) out[cbind(i + 1L, rep.int(seq_len(d[[2L]]), diff(p)))] <- x
    if (is.list(dn) && length(dn) == 2L) dimnames(out) <- dn
    return(out)
  }
  se_refuse(sprintf(paste(
    "Its assay '%s' is stored as '%s' (on disk, or in a form this reader does",
    "not know). Turn it into an ordinary matrix in R first, e.g. with",
    "assay(x, \"%s\") <- as.matrix(assay(x, \"%s\")), and save the object again."),
    name, if (is.na(cls)) typeof(m) else cls, name, name))
}

# The assays, unconverted, by name.
se_assay_list <- function(x) {
  a <- se_slot(x, "assays")
  cls <- s4_class_name(a)
  if (!isS4(a) || !identical(cls, "SimpleAssays")) {
    se_refuse(sprintf(paste(
      "Its assays are stored as '%s', the layout of an older SummarizedExperiment.",
      "Update the object in R with updateObject() and save it again."),
      if (is.na(cls)) typeof(a) else cls))
  }
  d <- se_slot(a, "data")
  if (!isS4(d) || !identical(s4_class_name(d), "SimpleList")) {
    se_refuse("Its assays are not stored as a list.")
  }
  ld <- se_slot(d, "listData")
  if (!is_plain_list(ld) || !length(ld)) se_refuse("It holds no assay.")
  nms <- attr(ld, "names", exact = TRUE)
  if (typeof(nms) != "character" || length(nms) != length(ld)) {
    nms <- rep("", length(ld))
  }
  nms[is.na(nms) | !nzchar(nms)] <- paste0("assay", which(is.na(nms) | !nzchar(nms)))
  attr(ld, "names") <- nms
  ld
}

# Feature names: the NAMES slot, or the names of rowRanges.
se_row_names <- function(x) {
  nm <- se_slot(x, "NAMES")
  if (typeof(nm) == "character" && !isS4(nm)) return(as.character(nm))
  rr <- se_slot(x, "rowRanges")
  holder <- switch(s4_class_name(rr) %|NA|% "",
                   GRanges = se_slot(rr, "ranges"),
                   CompressedGRangesList = se_slot(rr, "partitioning"),
                   NULL)
  if (isS4(holder)) {
    nm <- se_slot(holder, "NAMES")
    if (typeof(nm) == "character" && !isS4(nm)) return(as.character(nm))
  }
  NULL
}

# The feature annotation: rowRanges' metadata columns for a ranged
# object, the elementMetadata slot otherwise. A rowRanges of another
# shape is left alone (and said), not refused: the matrix does not need
# it.
se_row_table <- function(x) {
  rr <- se_slot(x, "rowRanges")
  if (is.null(rr)) return(list(df = se_slot(x, "elementMetadata"), note = NULL))
  if (isS4(rr) && s4_class_name(rr) %in% c("GRanges", "CompressedGRangesList")) {
    return(list(df = se_slot(rr, "elementMetadata"), note = NULL))
  }
  list(df = NULL, note = sprintf(
    "The feature annotation (rowRanges, stored as '%s') was not read.",
    s4_class_name(rr) %|NA|% typeof(rr)))
}

# DESeq2 records which of rowData's columns it computed itself (base
# means, dispersions, test results) in the columns' own metadata, as
# type "intermediate" or "results". Those are a fit, not annotation.
# The metadata table has one row per column, in the columns' order.
deseq_computed_columns <- function(row_df) {
  meta <- tryCatch(se_data_frame(se_slot(row_df, "elementMetadata"), "row annotation")$df,
                   omics_se_error = function(e) NULL)
  if (is.null(meta) || !"type" %in% names(meta)) return(character(0))
  cols <- attr(se_slot(row_df, "listData"), "names", exact = TRUE)
  if (typeof(cols) != "character" || length(cols) != nrow(meta)) return(character(0))
  cols[as.character(meta$type) %in% c("intermediate", "results")]
}

#' The parts of a SummarizedExperiment, read slot by slot
#'
#' @return `list(class, assays, row_names, col_df, row_df, notes)`; the
#'   assays still unconverted, by name.
#' @keywords internal
#' @noRd
se_parts_from_slots <- function(x) {
  notes <- character(0)
  assays <- se_assay_list(x)
  col <- se_data_frame(se_slot(x, "colData"), "sample information (colData)")
  rt <- se_row_table(x)
  if (!is.null(rt$note)) notes <- c(notes, rt$note)
  drop_cols <- character(0)
  if (identical(s4_class_name(x), "DESeqDataSet") && isS4(rt$df)) {
    drop_cols <- deseq_computed_columns(rt$df)
  }
  row <- se_data_frame(rt$df, "feature annotation (rowData)")
  if (length(drop_cols) && !is.null(row$df)) {
    src <- attr(row$df, "source_rownames")
    row$df <- row$df[, setdiff(names(row$df), drop_cols), drop = FALSE]
    attr(row$df, "source_rownames") <- src
    row$skipped <- setdiff(row$skipped, drop_cols)
    notes <- c(notes, sprintf(
      "Left out %d feature column(s) DESeq2 had computed (base means, dispersions, test results); they are recomputed when DESeq2 runs.",
      length(drop_cols)))
  }
  for (part in list(list(col, "sample information"), list(row, "feature annotation"))) {
    if (length(part[[1L]]$skipped)) {
      notes <- c(notes, sprintf(
        "Left out %d column(s) of the %s that are not plain values: %s.",
        length(part[[1L]]$skipped), part[[2L]],
        paste(utils::head(part[[1L]]$skipped, 6L), collapse = ", ")))
    }
  }
  list(class = s4_class_name(x), assays = assays, row_names = se_row_names(x),
       col_df = col$df, row_df = row$df, notes = notes)
}

# The same parts through SummarizedExperiment's own accessors, for an
# object handed over in R whose slots are laid out some other way (an
# assay kept on disk, an object saved by an old release). Never for an
# uploaded file: these run the object's methods.
se_parts_from_accessors <- function(x) {
  an <- SummarizedExperiment::assayNames(x)
  n <- length(SummarizedExperiment::assays(x, withDimnames = FALSE))
  if (!n) se_refuse("It holds no assay.")
  if (is.null(an) || length(an) != n) an <- paste0("assay", seq_len(n))
  assays <- stats::setNames(lapply(seq_len(n), function(i) {
    as.matrix(SummarizedExperiment::assay(x, i, withDimnames = TRUE))
  }), an)
  col <- se_data_frame(SummarizedExperiment::colData(x), "sample information (colData)")
  row <- se_data_frame(SummarizedExperiment::rowData(x), "feature annotation (rowData)")
  notes <- character(0)
  if (length(col$skipped) || length(row$skipped)) {
    notes <- sprintf("Left out %d column(s) of the sample and feature tables that are not plain values.",
                     length(col$skipped) + length(row$skipped))
  }
  list(class = class(x)[[1L]], assays = assays, row_names = rownames(x),
       col_df = col$df, row_df = row$df, notes = notes)
}

# Which assay is the matrix.
choose_se_assay <- function(names, assay = NULL) {
  if (!is.null(assay)) {
    i <- if (is.numeric(assay)) {
      if (assay >= 1 && assay <= length(names) && assay == round(assay)) as.integer(assay) else NA
    } else {
      match(as.character(assay), names)
    }
    if (is.na(i)) {
      stop(sprintf("`assay` must name one of the object's assays (%s), not %s.",
                   paste(sprintf("'%s'", names), collapse = ", "), describe_value(assay)),
           call. = FALSE)
    }
    return(list(i = i, why = "as asked"))
  }
  if ("counts" %in% names) return(list(i = match("counts", names), why = "the read counts"))
  list(i = 1L, why = "the first")
}

# ---- building the omics_input ----------------------------------------------

#' A SummarizedExperiment or DESeqDataSet as an omics_input
#'
#' @param x The object.
#' @param trusted Whether `x` was handed over in R (`TRUE`) or read from
#'   an uploaded file (`FALSE`; see the top of this file).
#' @return `list(input, report)`, as [read_omics()].
#' @keywords internal
#' @noRd
read_se_object <- function(x, omics_type = NULL, assay_type = NULL, assay = NULL,
                           source = NA_character_, trusted = FALSE) {
  cls <- s4_class_name(x)
  fail <- function(msg) {
    list(input = NULL, report = new_import_report(
      sheets = data.frame(name = cls, role = "unknown", n_rows = 0L, n_cols = 0L,
                          confidence = 0, orientation = NA_character_,
                          notes = "could not be read", stringsAsFactors = FALSE),
      warnings = sprintf("This %s could not be read. %s", cls, msg),
      source = source))
  }
  unsafe_message <- function(e) {
    sub("which a project never holds", "which a data table never holds",
        conditionMessage(e), fixed = TRUE)
  }
  parts <- tryCatch(se_parts_from_slots(x), omics_se_error = function(e) e,
                    omp_unsafe_error = function(e) {
                      structure(class = c("omp_unsafe_error", "error", "condition"),
                                list(message = unsafe_message(e), call = NULL))
                    })
  # Code where data should be is final; a layout this reader does not
  # know is not, for an object handed over in R.
  if (inherits(parts, "omp_unsafe_error")) return(fail(conditionMessage(parts)))
  if (inherits(parts, "omics_se_error") && isTRUE(trusted) &&
      is_installed("SummarizedExperiment")) {
    parts <- tryCatch(se_parts_from_accessors(x), error = function(e) {
      structure(class = c("omics_se_error", "error", "condition"),
                list(message = unsafe_message(e), call = NULL))
    })
  }
  if (inherits(parts, "omics_se_error")) return(fail(conditionMessage(parts)))
  # One assay as a matrix; for an object handed over in R, an assay in
  # a form se_assay_matrix() does not know (a DelayedMatrix, kept on
  # disk) is asked for through the object's own accessor.
  assay_matrix <- function(i) {
    name <- names(parts$assays)[[i]]
    out <- tryCatch(se_assay_matrix(parts$assays[[i]], name), omics_se_error = function(e) e)
    if (inherits(out, "omics_se_error") && isTRUE(trusted) &&
        is_installed("SummarizedExperiment")) {
      dense <- tryCatch(as.matrix(SummarizedExperiment::assay(x, i, withDimnames = FALSE)),
                        error = function(e) NULL)
      if (!is.null(dense)) {
        out <- tryCatch(se_assay_matrix(dense, name), omics_se_error = function(e) e)
      }
    }
    out
  }

  notes <- character(0)
  pick <- choose_se_assay(names(parts$assays), assay)
  assay_name <- names(parts$assays)[[pick$i]]
  mat <- assay_matrix(pick$i)
  if (inherits(mat, "omics_se_error")) return(fail(conditionMessage(mat)))
  if (length(parts$assays) > 1L) {
    notes <- c(notes, sprintf(
      "The object holds %d assays (%s); '%s' was read as the matrix (%s).",
      length(parts$assays), paste(names(parts$assays), collapse = ", "),
      assay_name, pick$why))
  }

  # Names: the object's, then the assay's own, then made up.
  rn <- parts$row_names %||% rownames(mat) %||%
    attr(parts$row_df, "source_rownames")
  if (length(rn) != nrow(mat)) rn <- NULL
  cn <- attr(parts$col_df, "source_rownames") %||% colnames(mat)
  if (length(cn) != ncol(mat)) cn <- NULL
  if (is.null(rn)) {
    rn <- paste0("feature_", seq_len(nrow(mat)))
    notes <- c(notes, "The features had no names, so they were numbered feature_1 onwards.")
  }
  if (is.null(cn)) {
    cn <- paste0("sample_", seq_len(ncol(mat)))
    notes <- c(notes, "The samples had no names, so they were numbered sample_1 onwards.")
  }
  if (anyDuplicated(rn) || anyNA(rn) || !all(nzchar(rn))) {
    notes <- c(notes, sprintf("%d feature name(s) were repeated or empty and were made unique.",
                              sum(duplicated(rn) | is.na(rn) | !nzchar(rn))))
  }
  if (anyDuplicated(cn) || anyNA(cn) || !all(nzchar(cn))) {
    notes <- c(notes, "Some sample names were repeated or empty and were made unique.")
  }
  rn <- make_unique_labels(rn)
  cn <- make_unique_labels(cn)
  dimnames(mat) <- list(rn, cn)

  # Read counts are RNA-seq, whatever the caller expected.
  counts <- identical(assay_name, "counts")
  asked <- omics_type
  if (counts) {
    omics_type <- "rnaseq"
    if (!is.null(asked) && !identical(asked, "rnaseq")) {
      notes <- c(notes, sprintf(
        "The '%s' assay holds RNA-seq read counts, so this was read as RNA-seq, not %s.",
        assay_name, asked))
      assay_type <- NULL
    }
  }

  # Sample information, one row per sample in the matrix's order.
  meta <- parts$col_df
  if (is.null(meta) || nrow(meta) != ncol(mat)) {
    if (!is.null(meta)) notes <- c(notes, "The sample information did not have one row per sample and was not used.")
    meta <- data.frame(row.names = seq_len(ncol(mat)))
  }
  attr(meta, "source_rownames") <- NULL
  if (identical(cls, "DESeqDataSet") && "sizeFactor" %in% names(meta)) {
    # DESeq2's own estimate, re-made when DESeq2 runs; as a column it
    # would be offered as a covariate.
    meta$sizeFactor <- NULL
  }
  if (!"sample_id" %in% names(meta)) meta <- cbind(data.frame(sample_id = cn, stringsAsFactors = FALSE), meta)
  rownames(meta) <- cn

  # GENCODE's empty chromosome-Y copies (see feature-symbols.R).
  par <- par_y_rows(mat)
  notes <- c(notes, par$notes)

  # Feature annotation, one row per feature.
  feat <- data.frame(feature_id = rn, stringsAsFactors = FALSE)
  row_df <- parts$row_df
  if (!is.null(row_df) && nrow(row_df) == nrow(mat) && ncol(row_df) > 0L) {
    attr(row_df, "source_rownames") <- NULL
    names(row_df)[names(row_df) == "feature_id"] <- "source_feature_id"
    rownames(row_df) <- NULL
    feat <- cbind(feat, row_df)
    sym_col <- pick_symbol_column(feat, exclude = "feature_id")
    if (!is.null(sym_col)) feat$feature_symbol <- first_gene_symbol(feat[[sym_col]])
  }
  rownames(feat) <- rn
  if (!all(par$keep)) {
    mat <- mat[par$keep, , drop = FALSE]
    feat <- feat[par$keep, , drop = FALSE]
  }
  sym <- attach_gene_symbols(feat, rownames(mat))
  feat <- sym$feature_df
  if (!is.null(sym$note)) notes <- c(notes, sym$note)

  # Effective lengths, when the object carries them beside its counts:
  # tximeta's "length" assay, DESeqDataSetFromTximport's "avgTxLength".
  len <- NULL
  if (counts) {
    len_name <- intersect(c("avgTxLength", "length"), names(parts$assays))[1L]
    if (!is.na(len_name)) {
      len <- assay_matrix(match(len_name, names(parts$assays)))
      if (is.matrix(len) && identical(dim(len), c(length(rn), length(cn)))) {
        dimnames(len) <- list(rn, cn)
        len <- len[par$keep, , drop = FALSE]
        len[!is.finite(len) | len <= 0] <- 1
        notes <- c(notes, sprintf(
          "The '%s' assay was kept as the features' effective lengths, which DESeq2 and edgeR correct for.",
          len_name))
      } else {
        len <- NULL
      }
    }
  }

  # What was copied out is plain data; checked as an uploaded table is.
  unsafe <- tryCatch({
    check_project_structure(list(mat, meta, feat, len))
    NULL
  }, omp_unsafe_error = function(e) conditionMessage(e))
  if (!is.null(unsafe)) {
    return(fail(sub("which a project never holds", "which a data table never holds",
                    unsafe, fixed = TRUE)))
  }

  if (is.null(assay_type) && !is.null(omics_type)) {
    assay_type <- infer_assay_type(mat, omics_type)
    if (is.na(assay_type)) assay_type <- if (counts) "raw_count" else NULL
  }

  sheet_table <- data.frame(
    name = c(names(parts$assays), "colData", "rowData"),
    role = c(ifelse(seq_along(parts$assays) == pick$i, "matrix", "unknown"),
             "metadata", "feature_annot"),
    n_rows = c(rep(nrow(mat), length(parts$assays)), ncol(mat), nrow(mat)),
    n_cols = c(rep(ncol(mat), length(parts$assays)), ncol(meta), ncol(feat)),
    confidence = 1,
    orientation = c(ifelse(seq_along(parts$assays) == pick$i, "features_in_rows", NA_character_),
                    NA_character_, NA_character_),
    notes = c(ifelse(seq_along(parts$assays) == pick$i, sprintf("assay of a %s", cls),
                     "another assay, not read"),
              "sample information", "feature annotation"),
    stringsAsFactors = FALSE)
  report <- new_import_report(
    sheets = sheet_table,
    warnings = c(sprintf("Read a %s: the '%s' assay of %s features and %d samples.",
                         cls, assay_name, format(nrow(mat), big.mark = ","), ncol(mat)),
                 parts$notes, notes),
    suggested_input = list(matrix_sheet = assay_name, metadata_sheet = "colData",
                           feature_sheet = "rowData", orientation = "features_in_rows",
                           orientation_confidence = 1, orientation_source = "detected",
                           omics_type = omics_type, assay_type = assay_type,
                           se_class = cls, se_assay = assay_name),
    source = source)

  if (is.null(omics_type)) {
    report <- add_import_warning(report,
      "`omics_type` not specified; caller must set it before constructing the input.")
    return(list(input = NULL, report = report))
  }
  input <- tryCatch(
    omics_input(mat, meta, feat, omics_type = omics_type, assay_type = assay_type),
    error = function(e) {
      report <<- add_import_warning(report,
        paste0("omics_input() rejected the assembly: ", conditionMessage(e)))
      NULL
    })
  if (!is.null(input) && !is.null(len)) {
    input$misc <- list(tximport = list(length = len, counts_from_abundance = "no",
                                       source = cls))
  }
  list(input = input, report = report)
}

#' Read a SummarizedExperiment or DESeqDataSet as an omics layer
#'
#' Converts a Bioconductor `SummarizedExperiment` (or
#' `RangedSummarizedExperiment`, or a `DESeqDataSet`) into an
#' [omics_input()], with an [`ImportReport`][new_import_report()] saying
#' what was taken:
#'
#' * one assay as the matrix -- the `counts` assay when there is one (a
#'   `DESeqDataSet` always has), otherwise the first, or the one named by
#'   `assay`; the report names the assays left out;
#' * `colData` as the sample information, and `rowData` as the feature
#'   annotation, with any column of gene symbols picked up as for a
#'   workbook and Ensembl ids mapped to symbols as for any import;
#' * for counts, the effective lengths beside them when the object carries
#'   them (tximeta's `length` assay, `DESeqDataSetFromTximport()`'s
#'   `avgTxLength`), kept for DESeq2 and edgeR.
#'
#' Read counts are RNA-seq, so a `counts` assay sets `omics_type` to
#' `"rnaseq"`. The assay type is inferred from the values with
#' [infer_assay_type()] unless given. Only columns of plain values are
#' copied out of `colData` and `rowData`; DESeq2's size factors and the
#' columns it computed (dispersions, test results) are left out, since
#' every analysis re-estimates them.
#'
#' The same conversion runs when [read_omics()] is given an `.rds` file
#' holding one of these objects; there the object is read slot by slot and
#' its code (a `DESeqDataSet`'s design formula and dispersion function) is
#' never touched.
#'
#' @param x A `SummarizedExperiment`, `RangedSummarizedExperiment` or
#'   `DESeqDataSet`.
#' @param omics_type Optional omics modality, one of
#'   [SUPPORTED_OMICS_TYPES]. Ignored for a `counts` assay (RNA-seq).
#'   With neither, `input` is `NULL` and the report says so, as for
#'   [read_omics()].
#' @param assay_type Optional assay semantic label; inferred from the
#'   values when `NULL`.
#' @param assay Optional name or position of the assay to read.
#'
#' @return A list with two elements, like [read_omics()]:
#'   * `input`: an `omics_input`, or `NULL`;
#'   * `report`: an [`ImportReport`][new_import_report()].
#' @export
#' @family io
#' @examples
#' \dontrun{
#' dds <- DESeq2::DESeqDataSetFromMatrix(counts, coldata, design = ~ group)
#' res <- read_summarized_experiment(dds)
#' res$report
#' run_diff(res$input, method = "deseq2", group_col = "group",
#'          control_group = "ctrl", case_group = "ko")
#' }
read_summarized_experiment <- function(x, omics_type = NULL, assay_type = NULL,
                                       assay = NULL) {
  if (!is_se_object(x, trusted = TRUE)) {
    arg_stop("x", "a SummarizedExperiment or DESeqDataSet", x)
  }
  assert_choice(omics_type, "omics_type", SUPPORTED_OMICS_TYPES, allow_null = TRUE)
  assert_string(assay_type, "assay_type", allow_null = TRUE)
  assert_label(assay, "assay", allow_null = TRUE)
  read_se_object(x, omics_type = omics_type, assay_type = assay_type, assay = assay,
                 trusted = TRUE)
}
