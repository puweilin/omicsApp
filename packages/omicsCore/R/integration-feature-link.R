# Working out which feature of one layer is the same gene as which
# feature of another.
#
# Until this existed the two layers were joined on the gene symbol each
# carried, one row per symbol: where several proteins named one gene
# (isoforms such as P04637 and P04637-2, or several protein groups) all
# but the most abundant were dropped before anything was compared. Two
# things change here.
#
#   1. Every feature is kept. A gene measured by two proteins gives two
#      protein-gene pairs, and each method says what it does with them:
#      concordance and correlation report each pair as a row of its own;
#      ActivePathways, which scores genes, takes each gene's most
#      significant feature in a layer and corrects its p-value for how
#      many features were on offer (see run_integration_active_pathways()).
#
#   2. The join can be stated rather than inferred: a `feature_link`
#      table, one column per layer, each row saying "this feature of one
#      layer is this feature of the other" (a UniProt accession and a gene
#      symbol, say). When one is given it replaces symbol matching
#      outright -- a feature the table does not mention is not paired --
#      so what was integrated is exactly what the table says.
#
# A value in the table names a feature of its layer when it is, in order
# of preference: the feature's id (or one member of a `;`-separated
# protein group id); the canonical accession of an isoform id
# (P04637-2 -> P04637), so that a table written at accession level covers
# every isoform of the protein; or the feature's symbol. Each feature
# takes the most specific of these that the table offers, so an entry for
# P04637-2 itself wins over the one for P04637.

# A UniProt accession with an isoform suffix, captured without it. The
# accession part is UniProt's own format, so a gene symbol that happens
# to end in -<digit> (NKX2-1, HLA-DRB1-3) is not mistaken for one.
UNIPROT_ISOFORM_PATTERN <-
  "^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})-[0-9]+$"

# The canonical accession of each isoform accession; NA for anything else.
uniprot_canonical <- function(x) {
  x <- toupper(trimws(as.character(x)))
  out <- rep(NA_character_, length(x))
  hit <- !is.na(x) & grepl(UNIPROT_ISOFORM_PATTERN, x)
  out[hit] <- sub("-[0-9]+$", "", x[hit])
  out
}

# Every key a feature answers to, with its tier (1 = its id or a member
# of it, 2 = the canonical accession of an isoform id, 3 = its symbol).
feature_link_keys <- function(ids, symbols) {
  n <- length(ids)
  ids <- as.character(ids)
  pieces <- strsplit(ifelse(is.na(ids), "", ids), ";", fixed = TRUE)
  idx <- rep(seq_len(n), lengths(pieces))
  acc <- trimws(unlist(pieces, use.names = FALSE))
  # FASTA-style headers (sp|P04637|P53_HUMAN) carry the accession inside.
  acc <- sub("^(sp|tr)\\|([^|]+)\\|.*$", "\\2", acc)
  canon <- uniprot_canonical(acc)
  out <- rbind(
    data.frame(i = c(seq_len(n), idx), key = integration_join_key(c(ids, acc)),
               tier = 1L, stringsAsFactors = FALSE),
    data.frame(i = idx, key = integration_join_key(canon), tier = 2L,
               stringsAsFactors = FALSE),
    data.frame(i = seq_len(n), key = integration_join_key(symbols), tier = 3L,
               stringsAsFactors = FALSE)
  )
  out <- out[!is.na(out$key), , drop = FALSE]
  out[!duplicated(out[c("i", "key")]), , drop = FALSE]
}

# Which features of one layer each link row names: data.frame(i, row,
# tier), `tier` saying how (see feature_link_keys()).
match_link_side <- function(ids, symbols, values) {
  pieces <- strsplit(ifelse(is.na(values), "", as.character(values)), ";", fixed = TRUE)
  lk <- data.frame(row = rep(seq_along(values), lengths(pieces)),
                   key = integration_join_key(unlist(pieces, use.names = FALSE)),
                   stringsAsFactors = FALSE)
  lk <- unique(lk[!is.na(lk$key), , drop = FALSE])
  m <- merge(feature_link_keys(ids, symbols), lk, by = "key")
  if (!nrow(m)) return(data.frame(i = integer(0), row = integer(0), tier = integer(0)))
  best <- stats::ave(m$tier, m$i, FUN = min)
  unique(m[m$tier == best, c("i", "row", "tier"), drop = FALSE])
}

#' Pair the features of two layers
#'
#' @param ids_a,ids_b Feature ids of each layer.
#' @param sym_a,sym_b The join column of each layer (symbols), aligned
#'   with the ids; `NULL` when a layer has none.
#' @param link `NULL` to match symbols, or a two-column data frame whose
#'   first column names features of layer A and second of layer B.
#' @return A data frame with one row per pair: `i_a`, `i_b` (row indices
#'   into the inputs), `feature_a`, `feature_b`, `symbol_a`, `symbol_b`,
#'   `label` (the gene the pair stands for) and `feature_id` (the label,
#'   made unique where one gene has several pairs). Attribute `info`
#'   carries the counts [feature_link_counts()] makes.
#' @keywords internal
#' @noRd
link_features <- function(ids_a, sym_a, ids_b, sym_b, link = NULL) {
  ids_a <- as.character(ids_a)
  ids_b <- as.character(ids_b)
  sym_a <- if (is.null(sym_a)) rep(NA_character_, length(ids_a)) else as.character(sym_a)
  sym_b <- if (is.null(sym_b)) rep(NA_character_, length(ids_b)) else as.character(sym_b)

  if (is.null(link)) {
    a <- data.frame(i_a = seq_along(ids_a), key = integration_join_key(sym_a))
    b <- data.frame(i_b = seq_along(ids_b), key = integration_join_key(sym_b))
    pairs <- merge(a[!is.na(a$key), , drop = FALSE],
                   b[!is.na(b$key), , drop = FALSE], by = "key")
    pairs <- pairs[c("i_a", "i_b")]
    b_names_gene <- rep(FALSE, nrow(pairs))
  } else {
    ma <- match_link_side(ids_a, sym_a, link[[1L]])
    mb <- match_link_side(ids_b, sym_b, link[[2L]])
    names(ma) <- c("i_a", "row", "tier_a")
    names(mb) <- c("i_b", "row", "tier_b")
    pairs <- merge(ma, mb, by = "row")
    pairs <- pairs[!duplicated(pairs[c("i_a", "i_b")]), , drop = FALSE]
    # The gene a pair stands for is the symbol the link names it by: when
    # only layer B's side of the link is written in symbols (accessions
    # against gene names), B's symbol, whatever A's own symbol column
    # says. Otherwise A's, as with symbol matching.
    b_names_gene <- pairs$tier_b == 3L & pairs$tier_a != 3L
  }

  blank <- function(x) is.na(x) | !nzchar(trimws(x))
  label <- ifelse(b_names_gene, sym_b[pairs$i_b], sym_a[pairs$i_a])
  label[blank(label)] <- sym_b[pairs$i_b][blank(label)]
  label[blank(label)] <- sym_a[pairs$i_a][blank(label)]
  label[blank(label)] <- ids_a[pairs$i_a][blank(label)]
  out <- data.frame(
    i_a = pairs$i_a, i_b = pairs$i_b,
    feature_a = ids_a[pairs$i_a], feature_b = ids_b[pairs$i_b],
    symbol_a = sym_a[pairs$i_a], symbol_b = sym_b[pairs$i_b],
    label = label,
    stringsAsFactors = FALSE
  )
  out <- out[order(integration_join_key(out$label), out$i_a, out$i_b), , drop = FALSE]
  rownames(out) <- NULL
  out$feature_id <- unique_pair_ids(out)
  attr(out, "info") <- feature_link_counts(out, ids_a, ids_b)
  out
}

# One id per pair. The gene's own name where it has one pair -- which is
# what every result held before isoforms were kept, so those results are
# unchanged -- and otherwise the gene followed by whichever side's
# features tell its pairs apart: "TP53 (P04637-2)".
unique_pair_ids <- function(pairs) {
  ids <- pairs$label
  if (!nrow(pairs)) return(character(0))
  key <- integration_join_key(ids)
  key[is.na(key)] <- ""
  dup <- key %in% key[duplicated(key)]
  for (k in unique(key[dup])) {
    rows <- which(key == k)
    va <- length(unique(pairs$feature_a[rows])) > 1L
    vb <- length(unique(pairs$feature_b[rows])) > 1L
    parts <- if (va && vb) paste(pairs$feature_a[rows], pairs$feature_b[rows], sep = " / ")
             else if (va) pairs$feature_a[rows]
             else pairs$feature_b[rows]
    ids[rows] <- sprintf("%s (%s)", pairs$label[rows], parts)
  }
  make.unique(ids, sep = " #")
}

# How the pairing went, in counts: features of each layer, how many were
# paired, how many pairs, and how many features share their partner with
# another feature of their layer (several proteins of one gene).
feature_link_counts <- function(pairs, ids_a, ids_b) {
  per_b <- table(pairs$i_b)
  per_a <- table(pairs$i_a)
  multi_b <- as.integer(names(per_b)[per_b > 1L])
  multi_a <- as.integer(names(per_a)[per_a > 1L])
  list(
    n_features_a = length(ids_a),
    n_features_b = length(ids_b),
    n_pairs = nrow(pairs),
    n_paired_a = length(unique(pairs$i_a)),
    n_paired_b = length(unique(pairs$i_b)),
    # Features of A sharing their B partner with another feature of A,
    # and the B features so shared.
    n_a_sharing = length(unique(pairs$i_a[pairs$i_b %in% multi_b])),
    n_b_shared = length(multi_b),
    n_b_sharing = length(unique(pairs$i_b[pairs$i_a %in% multi_a])),
    n_a_shared = length(multi_a)
  )
}

# ---- the link table ------------------------------------------------------

# Internal: check a feature_link against the project's layer names. Mirrors
# validate_sample_link(): only the shape is checked when the project has
# no layers yet.
validate_feature_link <- function(feature_link, tags) {
  if (!is.data.frame(feature_link)) {
    stop("`feature_link` must be a data.frame or NULL.")
  }
  if (ncol(feature_link) < 2L) {
    stop("`feature_link` needs a column for each of two layers, named after ",
         "the layers (e.g. data.frame(prot = \"P04637\", rna = \"TP53\")).",
         call. = FALSE)
  }
  if (length(tags) > 0L && sum(names(feature_link) %in% tags) < 2L) {
    stop("`feature_link` columns must be named after the layers they map; ",
         "the project's layers are: ", paste(tags, collapse = ", "), ".",
         call. = FALSE)
  }
  invisible(TRUE)
}

# The two columns of a link for one pair of layers, as clean character
# vectors: blank cells and incomplete or repeated rows dropped.
link_columns <- function(feature_link, tag_a, tag_b) {
  clean <- function(x) {
    x <- trimws(as.character(x))
    x[!is.na(x) & !nzchar(x)] <- NA_character_
    x
  }
  out <- data.frame(a = clean(feature_link[[tag_a]]), b = clean(feature_link[[tag_b]]),
                    stringsAsFactors = FALSE)
  out <- out[!is.na(out$a) & !is.na(out$b), , drop = FALSE]
  out <- unique(out)
  names(out) <- c(tag_a, tag_b)
  rownames(out) <- NULL
  out
}

# The link a run uses: the one passed to run_integration(), else the
# project's when it names both layers, else none (match symbols).
# Returns list(link, source, file): `source` is "supplied", "project" or
# "symbol"; `file` is where the table was read from, when known.
resolve_feature_link <- function(project, tag_a, tag_b, feature_link = NULL) {
  if (!is.null(feature_link)) {
    assert_data_frame(feature_link, "feature_link")
    missing <- setdiff(c(tag_a, tag_b), names(feature_link))
    if (length(missing)) {
      stop("`feature_link` must have a column named after each layer being ",
           "integrated; missing: ", paste(missing, collapse = ", "), ".",
           call. = FALSE)
    }
    source <- "supplied"
  } else {
    feature_link <- project$feature_link
    if (is.null(feature_link) || !all(c(tag_a, tag_b) %in% names(feature_link))) {
      return(list(link = NULL, source = "symbol", file = NULL))
    }
    source <- "project"
  }
  link <- link_columns(feature_link, tag_a, tag_b)
  if (!nrow(link)) {
    stop("The feature link has no row naming a feature of both '", tag_a,
         "' and '", tag_b, "'.", call. = FALSE)
  }
  list(link = link, source = source, file = attr(feature_link, "source"))
}

#' Read a feature link table
#'
#' Reads a table that says which feature of one layer is which feature of
#' another -- typically UniProt accessions against gene symbols or
#' Ensembl gene ids -- for [run_integration()]'s `feature_link`.
#' Delimited text (comma, tab, semicolon) and Excel files are read; an
#' Excel file is read from its first sheet.
#'
#' The result has one column per layer, named after the layer, holding
#' the two chosen columns of the file as text. Blank cells, rows missing
#' either side, and repeated rows are dropped. Where the file came from
#' is kept (attribute `"source"`), so [export_script()] can read the same
#' file again.
#'
#' A row may name several features of a layer separated by `;`. A
#' UniProt accession without an isoform suffix (`P04637`) covers every
#' isoform of the protein (`P04637-2`, ...) unless the table also names
#' the isoform itself.
#'
#' @param path Path to the file.
#' @param columns Named vector: names are the two layers, values the
#'   columns of the file holding each layer's identifiers, by name or by
#'   position, e.g. `c(prot = "UniProt", rna = "Gene")`. `NULL` returns
#'   every column of the file as text, unchanged, to choose from.
#' @return A data frame with one column per layer (with `columns = NULL`,
#'   the file's own columns).
#' @export
#' @family integration
#' @examples
#' f <- tempfile(fileext = ".csv")
#' write.csv(data.frame(UniProt = c("P04637", "P38398"),
#'                      Gene = c("TP53", "BRCA1")), f, row.names = FALSE)
#' read_feature_link(f, columns = c(prot = "UniProt", rna = "Gene"))
read_feature_link <- function(path, columns) {
  assert_string(path, "path")
  if (!file.exists(path)) stop("File not found: ", path, call. = FALSE)
  if (!is.null(columns) && (length(columns) != 2L || is.null(names(columns)) ||
      any(!nzchar(names(columns))) || anyDuplicated(names(columns)))) {
    stop("`columns` must name two layers, e.g. c(prot = \"UniProt\", rna = \"Gene\").",
         call. = FALSE)
  }
  guard_archive(path)
  tab <- if (is_excel_file(path)) {
    as.data.frame(read_excel_sheet(path, readxl::excel_sheets(path)[[1L]],
                                   col_types = "text"),
                  stringsAsFactors = FALSE)
  } else {
    read_sample_table(path)
  }
  attr(tab, "encoding") <- NULL
  if (is.null(columns)) {
    tab[] <- lapply(tab, function(x) trimws(as.character(x)))
    return(tab)
  }
  pick <- function(col) {
    if (is.numeric(col)) {
      if (col < 1L || col > ncol(tab)) {
        stop("The file has ", ncol(tab), " column(s); there is no column ", col, ".",
             call. = FALSE)
      }
      return(tab[[col]])
    }
    if (!col %in% names(tab)) {
      stop("The file has no column '", col, "'; its columns are: ",
           paste(names(tab), collapse = ", "), ".", call. = FALSE)
    }
    tab[[col]]
  }
  raw <- stats::setNames(
    data.frame(pick(columns[[1L]]), pick(columns[[2L]]), stringsAsFactors = FALSE),
    names(columns))
  out <- link_columns(raw, names(columns)[[1L]], names(columns)[[2L]])
  if (!nrow(out)) stop("The feature link file has no complete rows.", call. = FALSE)
  attr(out, "source") <- list(path = path, columns = columns)
  out
}

#' Preview how two layers' features will be paired
#'
#' The feature-level counterpart of [sample_pairing_preview()]: how
#' [run_integration()] will match the features of two layers, and what
#' came of it, without running anything.
#'
#' @param project An `omics_project`.
#' @param tag_a,tag_b The two layers.
#' @param feature_link Optional link table (see [run_integration()]); when
#'   `NULL`, the project's own `feature_link` is used if it names both
#'   layers, and symbols are matched otherwise.
#' @param by Column of each layer's `feature_df` matched when there is no
#'   link (default `"feature_symbol"`).
#' @return A list: `source` (`"supplied"`, `"project"` or `"symbol"`),
#'   `pairs` (data frame with `feature_a`, `feature_b`, `feature`), and the
#'   counts `n_features_a`, `n_features_b`, `n_pairs`, `n_paired_a`,
#'   `n_paired_b`, `n_a_sharing` (features of `tag_a` that share their
#'   partner with another feature of `tag_a`, as isoforms of one gene do),
#'   `n_b_shared` (the partners so shared), and the same two the other way
#'   round (`n_b_sharing`, `n_a_shared`).
#' @export
#' @family integration
feature_pairing_preview <- function(project, tag_a, tag_b, feature_link = NULL,
                                    by = "feature_symbol") {
  resolve_experiment_pair(project, c(tag_a, tag_b))
  res <- resolve_feature_link(project, tag_a, tag_b, feature_link)
  fa <- project$experiments[[tag_a]]$feature_df
  fb <- project$experiments[[tag_b]]$feature_df
  pairs <- link_features(fa$feature_id, fa[[by]], fb$feature_id, fb[[by]],
                         link = res$link)
  c(list(source = res$source,
         pairs = data.frame(feature_a = pairs$feature_a, feature_b = pairs$feature_b,
                            feature = pairs$feature_id, stringsAsFactors = FALSE)),
    attr(pairs, "info"))
}
