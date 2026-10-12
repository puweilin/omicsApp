# Per-sample quantification files: Salmon's quant.sf, RSEM's
# *.genes.results / *.isoforms.results, kallisto's abundance.tsv.
#
# Each holds one sample: a feature id, its lengths, and two or three
# measures of it (estimated reads, TPM, FPKM). Read as a matrix, the
# numeric columns became four "samples" called Length, EffectiveLength,
# TPM and NumReads, and the import analysed them. They are recognised by
# their column signature instead, and several of them are merged into one
# counts layer the way tximport does it -- estimated reads as the counts,
# effective lengths kept beside them so DESeq2 and edgeR can correct for
# the length each transcript or gene had in each sample. tximport itself
# is not used: it is a Bioconductor install for twenty lines of
# arithmetic, and it is deliberately not in the production image.

# The columns that identify each format, and which of them hold the
# estimated reads, the effective length and the abundance. RSEM names its
# id column gene_id or transcript_id depending on the level; the id is
# whichever comes first.
QUANT_FORMATS <- list(
  salmon = list(label = "Salmon",
                columns = c("Name", "Length", "EffectiveLength", "TPM", "NumReads"),
                ids = "Name", counts = "NumReads", length = "EffectiveLength",
                abundance = "TPM"),
  rsem = list(label = "RSEM",
              columns = c("length", "effective_length", "expected_count", "TPM"),
              ids = c("gene_id", "transcript_id"), counts = "expected_count",
              length = "effective_length", abundance = "TPM"),
  kallisto = list(label = "kallisto",
                  columns = c("target_id", "length", "eff_length", "est_counts", "tpm"),
                  ids = "target_id", counts = "est_counts", length = "eff_length",
                  abundance = "tpm")
)

#' Which quantification format a table is in
#'
#' @param x A data frame, or its column names.
#' @return `"salmon"`, `"rsem"`, `"kallisto"`, or `NA_character_`.
#' @keywords internal
#' @noRd
detect_quant_format <- function(x) {
  nms <- if (is.data.frame(x)) colnames(x) else as.character(x)
  if (!length(nms)) return(NA_character_)
  for (fmt in names(QUANT_FORMATS)) {
    spec <- QUANT_FORMATS[[fmt]]
    if (all(spec$columns %in% nms) && nms[[1L]] %in% spec$ids) return(fmt)
  }
  NA_character_
}

# File names the quantifiers write, taken off to leave the sample name:
# "ctrl_1.quant.sf", "ctrl_1.genes.results", "ctrl_1_abundance.tsv".
QUANT_FILE_SUFFIX_RE <- paste0(
  "([._-]?(quant([._]genes)?[.]sf|abundance[.](tsv|txt)|",
  "(genes|isoforms)[.]results|quant|abundance))?",
  "([.](sf|tsv|txt|csv|results))?([.]gz)?$")

#' Sample names from quantification file paths
#'
#' The file name without what the quantifier calls every file; for
#' Salmon's and kallisto's one-directory-per-sample layout
#' (`ctrl_1/quant.sf`), where nothing is left, the directory's name.
#'
#' @param paths File paths, or bare file names as a browser sends them.
#' @return Character vector, one name per path, possibly empty strings
#'   where neither the file nor its directory names the sample.
#' @keywords internal
#' @noRd
quant_sample_names <- function(paths) {
  paths <- gsub("\\\\", "/", as.character(paths))
  file <- sub("^.*/", "", paths)
  stem <- sub(QUANT_FILE_SUFFIX_RE, "", file, ignore.case = TRUE)
  dir <- ifelse(grepl("/", paths), sub("^.*/", "", sub("/[^/]*$", "", paths)), "")
  ifelse(nzchar(stem), stem, dir)
}

# One file, as the counts, effective lengths and abundance it holds,
# keyed by feature id.
read_quant_table <- function(path, shown = basename(path)) {
  read <- read_delimited_table(path)
  df <- read$df
  fmt <- detect_quant_format(df)
  if (is.na(fmt)) {
    stop(sprintf(paste(
      "'%s' is not a Salmon, RSEM or kallisto quantification file;",
      "its columns are %s."),
      shown, paste(utils::head(colnames(df), 8L), collapse = ", ")),
      call. = FALSE)
  }
  spec <- QUANT_FORMATS[[fmt]]
  ids <- trimws(as.character(df[[1L]]))
  num <- function(col) suppressWarnings(as.numeric(df[[col]]))
  list(format = fmt, ids = ids, counts = num(spec$counts),
       length = num(spec$length), abundance = num(spec$abundance),
       gene_id = if (identical(colnames(df)[[1L]], "transcript_id") &&
                     "gene_id" %in% colnames(df)) as.character(df$gene_id),
       encoding = read$encoding)
}

#' Read Salmon, RSEM or kallisto quantification files as one counts layer
#'
#' Merges one quantification file per sample -- Salmon's `quant.sf`,
#' RSEM's `*.genes.results` or `*.isoforms.results`, kallisto's
#' `abundance.tsv` -- into a single RNA-seq counts layer, the way
#' `tximport` does with `countsFromAbundance = "no"`:
#'
#' * the counts are the estimated reads (`NumReads`, `expected_count`,
#'   `est_counts`), which are fractional;
#' * the effective lengths are kept beside them, in the `tximport` slot
#'   that `run_diff(method = "deseq2")` and `run_diff(method = "edger")`
#'   read, so both correct for the length each feature had in each
#'   sample (a gene whose dominant isoform changes between groups
#'   otherwise looks differentially expressed);
#' * with `tx2gene`, transcripts are summed to genes: counts and TPM are
#'   summed, and a gene's length in a sample is the mean of its
#'   transcripts' effective lengths weighted by their TPM there.
#'
#' The `tximport` package itself is not needed.
#'
#' @param paths Paths to the quantification files, one per sample, all of
#'   one format and quantified against the same reference.
#' @param sample_names Optional sample names, one per file. By default
#'   each file's name without what the quantifier calls every file
#'   (`ctrl_1.genes.results` gives `ctrl_1`), or, for Salmon's and
#'   kallisto's one-directory-per-sample layout (`ctrl_1/quant.sf`), the
#'   directory's name.
#' @param tx2gene Optional data frame mapping transcripts to genes: the
#'   first column the transcript id, the second the gene id (the
#'   `tximport` convention). Ids are matched with and without their
#'   Ensembl version suffix; transcripts it does not list are left out,
#'   and the report says how many.
#' @param sample_sheet Optional path to a sample sheet, as for
#'   [read_omics()]: its rows are matched to the samples by name.
#' @param file_names Optional original names of the files, one per path,
#'   for when `paths` are temporary copies (a browser upload): sample
#'   names are then taken from these rather than from `paths`.
#'
#' @return A list with two elements, like [read_omics()]:
#'   * `input`: an `omics_input` with `omics_type = "rnaseq"` and
#'     `assay_type = "raw_count"`, or `NULL` when the files could not be
#'     merged;
#'   * `report`: an [`ImportReport`][new_import_report()] whose warnings
#'     say what was done.
#' @export
#' @family io
#' @examples
#' \dontrun{
#' files <- file.path("salmon", c("ctrl_1", "ctrl_2", "ko_1", "ko_2"), "quant.sf")
#' res <- read_quant_files(files, tx2gene = read.csv("tx2gene.csv"))
#' run_diff(res$input, method = "edger", group_col = "group",
#'          control_group = "ctrl", case_group = "ko")
#' }
read_quant_files <- function(paths, sample_names = NULL, tx2gene = NULL,
                             sample_sheet = NULL, file_names = NULL) {
  assert_names(paths, "paths")
  assert_character(sample_names, "sample_names", allow_null = TRUE)
  assert_character(file_names, "file_names", allow_null = TRUE)
  assert_string(sample_sheet, "sample_sheet", allow_null = TRUE)
  if (!is.null(tx2gene) && (!is.data.frame(tx2gene) || ncol(tx2gene) < 2L)) {
    arg_stop("tx2gene", "a data frame of transcript and gene ids", tx2gene)
  }
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop("File does not exist: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (!is.null(sample_sheet) && !file.exists(sample_sheet)) {
    stop("Sample sheet does not exist: ", sample_sheet, call. = FALSE)
  }
  for (arg in c("sample_names", "file_names")) {
    val <- get(arg)
    if (!is.null(val) && length(val) != length(paths)) {
      stop(sprintf("`%s` must name each of the %d file(s), not %d.",
                   arg, length(paths), length(val)), call. = FALSE)
    }
  }
  notes <- character(0)

  # Errors name the files as the user knows them, not as temporary copies.
  shown <- sub("^.*[\\\\/]", "", file_names %||% paths)
  tables <- unname(Map(read_quant_table, paths, shown))
  formats <- unique(vapply(tables, `[[`, character(1), "format"))
  if (length(formats) > 1L) {
    stop(sprintf(paste(
      "The files come from different quantifiers (%s); a layer is merged",
      "from files of one."), paste(vapply(formats, function(f) QUANT_FORMATS[[f]]$label,
                                          character(1)), collapse = ", ")),
      call. = FALSE)
  }
  fmt <- formats[[1L]]
  label <- QUANT_FORMATS[[fmt]]$label
  for (enc in unique(vapply(tables, `[[`, character(1), "encoding"))) {
    notes <- c(notes, encoding_note(enc, "The quantification file"))
  }

  # GENCODE transcript ids from Salmon run without --gencode carry the
  # whole FASTA header, "ENST...|ENSG...|...|DDX11L1-202|..."; the first
  # field is the id (tximport's ignoreAfterBar).
  ids <- tables[[1L]]$ids
  if (any(grepl("|", ids, fixed = TRUE))) {
    for (i in seq_along(tables)) tables[[i]]$ids <- sub("[|].*$", "", tables[[i]]$ids)
    ids <- tables[[1L]]$ids
    notes <- c(notes, "Feature ids carried '|'-separated annotation; the first field was kept as the id.")
  }
  # One reference for all: the rows are aligned by id, and files
  # quantified against different annotations cannot be merged.
  for (i in seq_along(tables)[-1L]) {
    if (!setequal(tables[[i]]$ids, ids) || length(tables[[i]]$ids) != length(ids)) {
      stop(sprintf(paste(
        "'%s' and '%s' list different features (%d vs %d); they were",
        "quantified against different references and cannot be merged."),
        shown[[1L]], shown[[i]],
        length(ids), length(tables[[i]]$ids)), call. = FALSE)
    }
  }
  if (anyDuplicated(ids)) {
    stop(sprintf("'%s' lists %d feature id(s) more than once, e.g. '%s'.",
                 shown[[1L]], sum(duplicated(ids)), ids[duplicated(ids)][[1L]]),
         call. = FALSE)
  }
  as_matrix <- function(what) {
    m <- vapply(tables, function(t) t[[what]][match(ids, t$ids)], numeric(length(ids)))
    matrix(m, nrow = length(ids), dimnames = list(ids, NULL))
  }
  counts <- as_matrix("counts")
  len <- as_matrix("length")
  abundance <- as_matrix("abundance")

  # Sample names.
  given <- !is.null(sample_names)
  samples <- if (given) trimws(sample_names) else quant_sample_names(file_names %||% paths)
  if (!all(nzchar(samples)) || anyDuplicated(samples)) {
    if (given) {
      stop("`sample_names` must be unique and non-empty.", call. = FALSE)
    }
    samples <- sprintf("sample_%d", seq_along(paths))
    notes <- c(notes, sprintf(paste(
      "The file names do not tell the samples apart, so they were named",
      "sample_1 to sample_%d in the order given. Rename the files to carry",
      "the sample names (e.g. ctrl_1.quant.sf) for a sample sheet to match them."),
      length(paths)))
  }
  colnames(counts) <- colnames(len) <- colnames(abundance) <- samples

  # Transcripts to genes.
  gene_of <- tables[[1L]]$gene_id
  if (!is.null(tx2gene)) {
    summ <- summarise_to_genes(counts, len, abundance, tx2gene)
    counts <- summ$counts
    len <- summ$length
    notes <- c(notes, summ$notes)
    gene_of <- NULL
  }
  if (anyNA(counts)) {
    stop("Some estimated read counts could not be read as numbers.", call. = FALSE)
  }
  # GENCODE's chromosome-Y copies of the pseudoautosomal genes, which
  # hold no reads (see feature-symbols.R). After summing to genes, so a
  # tx2gene that maps the copies to their own _PAR_Y genes drops those.
  par <- par_y_rows(counts)
  if (!all(par$keep)) {
    counts <- counts[par$keep, , drop = FALSE]
    len <- len[par$keep, , drop = FALSE]
    if (!is.null(gene_of)) gene_of <- gene_of[par$keep]
  }
  notes <- c(notes, par$notes)

  # A length of zero is what RSEM and kallisto write for a feature
  # shorter than the fragments; it has no reads, and a zero would make
  # the length offsets infinite. tximport's importers set it to one too.
  len[!is.finite(len) | len <= 0] <- 1

  feat <- data.frame(feature_id = rownames(counts), row.names = rownames(counts),
                     stringsAsFactors = FALSE)
  if (!is.null(gene_of)) feat$gene_id <- gene_of
  sym <- attach_gene_symbols(feat, rownames(counts))
  feat <- sym$feature_df
  if (!is.null(sym$note)) notes <- c(notes, sym$note)
  meta <- data.frame(sample_id = samples, row.names = samples, stringsAsFactors = FALSE)

  notes <- c(sprintf(paste(
    "Read %d %s quantification file(s): the estimated read counts (%s) of",
    "%s %s, with their effective lengths kept for DESeq2 and edgeR."),
    length(paths), label, QUANT_FORMATS[[fmt]]$counts,
    format(nrow(counts), big.mark = ","),
    if (!is.null(tx2gene)) "genes" else "features"), notes)
  if (length(paths) == 1L) {
    notes <- c(notes, paste(
      "This is one sample, and one sample cannot be compared with anything.",
      "Select the quantification files of all samples together to import",
      "them as one layer."))
  }

  sheet_table <- data.frame(
    name = samples, role = "matrix", n_rows = nrow(counts), n_cols = 1L,
    confidence = 1, orientation = "features_in_rows",
    notes = sprintf("%s quantification", label), stringsAsFactors = FALSE)
  report <- new_import_report(
    sheets = sheet_table, warnings = notes,
    suggested_input = list(matrix_sheet = samples[[1L]], orientation = "features_in_rows",
                           orientation_confidence = 1, orientation_source = "detected",
                           omics_type = "rnaseq", assay_type = "raw_count",
                           quant_format = fmt, quant_files = length(paths)),
    source = paste(shown, collapse = ", "))

  input <- tryCatch(
    suppressWarnings(omics_input(counts, meta, feat, omics_type = "rnaseq",
                                 assay_type = "raw_count")),
    error = function(e) {
      report <<- add_import_warning(report,
        paste0("omics_input() rejected the assembly: ", conditionMessage(e)))
      NULL
    })
  if (!is.null(input)) {
    input$misc <- list(tximport = list(length = len, counts_from_abundance = "no",
                                       source = fmt))
  }
  out <- list(input = input, report = report)
  if (!is.null(sample_sheet) && !is.null(input)) out <- attach_sample_sheet(out, sample_sheet)
  out
}

#' Transcript-level quantification summed to genes, as tximport does
#'
#' Counts and abundance are summed over a gene's transcripts. Its length
#' in a sample is its transcripts' effective lengths averaged with their
#' abundance there as weights -- the length a read from that gene came
#' from, on average. Where a gene has no abundance in a sample the
#' weights are all zero; tximport then takes the geometric mean of the
#' gene's length in the other samples, or, with no abundance anywhere,
#' the plain mean of its transcripts' lengths, and so does this.
#'
#' @return `list(counts, length, abundance, notes)`.
#' @keywords internal
#' @noRd
summarise_to_genes <- function(counts, len, abundance, tx2gene) {
  tx <- trimws(as.character(tx2gene[[1L]]))
  gene <- trimws(as.character(tx2gene[[2L]]))
  ids <- rownames(counts)
  notes <- character(0)
  idx <- match(ids, tx)
  # Versions dropped on both sides when that matches more: a tx2gene
  # built from one GENCODE release against an index from another, or
  # one side written without versions.
  idx_unv <- match(strip_ensembl_version(ids), strip_ensembl_version(tx))
  if (sum(!is.na(idx_unv)) > sum(!is.na(idx))) {
    idx <- idx_unv
    notes <- c(notes, "Transcripts were matched to tx2gene ignoring their version suffixes.")
  }
  if (all(is.na(idx))) {
    stop(sprintf(paste(
      "None of the transcript ids (e.g. '%s') is in the first column of",
      "`tx2gene` (e.g. '%s')."), ids[[1L]], tx[[1L]]), call. = FALSE)
  }
  if (anyNA(idx)) {
    notes <- c(notes, sprintf(
      "%s of %s transcripts are not in tx2gene and were left out (e.g. %s).",
      format(sum(is.na(idx)), big.mark = ","), format(length(ids), big.mark = ","),
      paste(utils::head(ids[is.na(idx)], 3L), collapse = ", ")))
  }
  keep <- !is.na(idx)
  g <- gene[idx[keep]]
  counts <- counts[keep, , drop = FALSE]
  len <- len[keep, , drop = FALSE]
  abundance <- abundance[keep, , drop = FALSE]

  # rowsum() orders the groups; that order is kept throughout.
  counts_g <- rowsum(counts, g, reorder = TRUE)
  abundance_g <- rowsum(abundance, g, reorder = TRUE)
  len_g <- rowsum(abundance * len, g, reorder = TRUE) / abundance_g
  # The plain mean of each transcript's length over samples, averaged
  # over the gene's transcripts: tximport's fallback when a gene has no
  # abundance in any sample.
  ave_len <- tapply(rowMeans(len), g, mean)[rownames(len_g)]
  for (i in which(apply(len_g, 1L, function(r) any(!is.finite(r))))) {
    bad <- !is.finite(len_g[i, ])
    len_g[i, bad] <- if (all(bad)) ave_len[[i]] else exp(mean(log(len_g[i, !bad])))
  }
  notes <- c(notes, sprintf(
    "Summed %s transcripts to %s genes with tx2gene; a gene's length in each sample is its transcripts' effective lengths weighted by their TPM.",
    format(nrow(counts), big.mark = ","), format(nrow(counts_g), big.mark = ",")))
  list(counts = counts_g, length = len_g, abundance = abundance_g, notes = notes)
}
