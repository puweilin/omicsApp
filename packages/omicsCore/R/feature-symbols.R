# Ensembl gene id -> HGNC symbol (human) or MGI symbol (mouse; see
# "mouse" below).
#
# An RNA-seq counts matrix is keyed on Ensembl ids, and every pathway
# database is keyed on symbols. Without a mapping, enrichment matches
# nothing and returns an empty result -- which is indistinguishable from
# "no pathway was enriched", so the failure arrives disguised as an
# answer.
#
# Vendors ship a `gene_name` column, but it is their annotation frozen on
# the day they ran the pipeline: symbols are renamed, merged and retired
# continuously, so a two-year-old export names genes that no database
# calls that any more. This maps from the id, which does not change.
#
# The table is bundled rather than fetched. The deployment's build host
# cannot reach most of the internet, and a runtime lookup would make
# every import depend on a network that is not there. 198 KB is a small
# price for an import that behaves the same on any machine.
#
# UNMAPPED IDS GET NA, DELIBERATELY, NOT THE ID ITSELF. About a sixth of
# the analysable genes in a human counts matrix have no HGNC symbol --
# unnamed lncRNAs, pseudogenes -- and no pathway database contains them
# either, so they cannot contribute to enrichment whatever we call them.
# Writing the id in would put 21,000 strings that can never match into
# ORA's universe, and the universe is the denominator of the
# hypergeometric test: inflating it makes every p-value look better than
# it is. NA keeps them out, and both enrichment paths already drop NA
# (enrich-ora.R universe/features, diff-utils.R ranked list).
#
# Differential analysis is unaffected: it works on feature_id, so an
# unnamed lncRNA is still tested and still reported.

HGNC_MAP_FILE <- "hgnc_ensembl.rds"

.hgnc_cache <- new.env(parent = emptyenv())

#' The bundled Ensembl-to-HGNC table
#'
#' @return A data frame with `ensembl_gene_id`, `symbol`, `locus_group`,
#'   carrying `source`, `upstream_modified` and `retrieved` attributes.
#'   `NULL` if the file is missing, which is not fatal -- callers fall
#'   back to leaving symbols alone.
#' @keywords internal
#' @noRd
hgnc_ensembl_map <- function() {
  if (!is.null(.hgnc_cache$map)) return(.hgnc_cache$map)
  path <- system.file("extdata", HGNC_MAP_FILE, package = "omicsCore")
  if (!nzchar(path) || !file.exists(path)) return(NULL)
  .hgnc_cache$map <- readRDS(path)
  .hgnc_cache$map
}

#' Where the symbol table came from
#'
#' Recorded in the import report and the analysis report, because a
#' result that depends on gene names should be able to say which
#' vintage of gene names produced it.
#'
#' @return A single string, or `NA_character_` when no table is bundled.
#' @export
#' @family enrich
hgnc_map_provenance <- function() {
  m <- hgnc_ensembl_map()
  if (is.null(m)) return(NA_character_)
  sprintf("%s (upstream %s, retrieved %s, %s symbols)",
          attr(m, "source") %||% "HGNC",
          attr(m, "upstream_modified") %||% "unknown",
          attr(m, "retrieved") %||% "unknown",
          format(nrow(m), big.mark = ","))
}

# ---- mouse -----------------------------------------------------------------
#
# Mouse Ensembl ids (ENSMUSG...) map to MGI symbols through babelgene's
# ortholog table, the one msigdbr itself uses to write MSigDB's gene sets
# in mouse symbols (Trp53 for TP53). It covers the mouse genes that have
# a human counterpart, about 20,000, and nothing else -- which is exactly
# the set the mouse gene sets can contain, since they are the human sets
# carried over gene by gene. A mouse gene outside it (a Gm-numbered gene,
# an unnamed lncRNA) gets NA for the reason human unmapped ids do: no
# gene set could match it, and an id in its place would only inflate
# ORA's universe.
#
# babelgene comes with msigdbr (< 10) and is listed in the production
# image; elsewhere it is optional, and without it mouse ids keep no
# symbol and the import says why.

#' Where the mouse symbols came from
#' @return A single string, or `NA_character_` without babelgene.
#' @keywords internal
#' @noRd
mgi_map_provenance <- function() {
  if (!is_installed("babelgene")) return(NA_character_)
  sprintf("babelgene %s (MGI symbols of the mouse genes with a human ortholog, from HCOP)",
          as.character(utils::packageVersion("babelgene")))
}

#' Map mouse Ensembl gene ids (unversioned) to MGI symbols
#' @param ids Unversioned ENSMUSG ids.
#' @return Character vector the same length; `NA` where unmapped or when
#'   babelgene is not installed.
#' @keywords internal
#' @noRd
map_mouse_ensembl <- function(ids) {
  out <- rep(NA_character_, length(ids))
  uniq <- unique(ids[!is.na(ids)])
  if (!length(uniq) || !is_installed("babelgene")) return(out)
  tab <- tryCatch(
    suppressWarnings(babelgene::orthologs(genes = uniq, species = "mouse", human = FALSE)),
    error = function(e) NULL)
  if (!is.data.frame(tab) || !nrow(tab) || !all(c("ensembl", "symbol") %in% names(tab))) {
    return(out)
  }
  # One mouse gene can be the ortholog of several human ones, and is
  # listed once per pairing; its symbol is the same in each but for a
  # handful of ids, where the best-supported pairing is taken.
  if ("support_n" %in% names(tab)) tab <- tab[order(-tab$support_n), , drop = FALSE]
  tab <- tab[!is.na(tab$ensembl) & !is.na(tab$symbol) & nzchar(tab$symbol), , drop = FALSE]
  tab <- tab[!duplicated(tab$ensembl), , drop = FALSE]
  idx <- match(ids, tab$ensembl)
  out[!is.na(idx)] <- as.character(tab$symbol[idx[!is.na(idx)]])
  out
}

# ---- ids --------------------------------------------------------------------

# `ENSG00000141510.17` and `ENSG00000141510` are the same gene; the
# suffix is the annotation version. Vendors differ on whether they keep
# it, and a table keyed on the unversioned id matches neither if we do
# not strip it. GENCODE writes the chromosome-Y copy of a
# pseudoautosomal gene as `ENSG00000182378.15_PAR_Y`; the version goes
# and the `_PAR_Y` stays, so the copy is still told apart from its X
# twin where ids are compared (tx2gene matching, for one).
strip_ensembl_version <- function(x) sub("\\.\\d+(_PAR_Y)?$", "\\1", as.character(x))

# The gene an id names, for looking its symbol up: no version, and the
# PAR_Y copy reduced to the gene it copies.
ensembl_base_id <- function(x) sub("_PAR_Y$", "", strip_ensembl_version(x))

ENSEMBL_GENE_RE <- c(
  human = "^ENSG\\d{6,}(\\.\\d+)?(_PAR_Y)?$",
  mouse = "^ENSMUSG\\d{6,}(\\.\\d+)?(_PAR_Y)?$"
)

#' Which species' Ensembl gene ids these are
#'
#' @param ids Character vector.
#' @param min_fraction How many must match to call the whole set Ensembl.
#' @return `"human"`, `"mouse"`, or `NA_character_`.
#' @keywords internal
#' @noRd
ensembl_species <- function(ids, min_fraction = 0.5) {
  ids <- ids[!is.na(ids)]
  if (length(ids) == 0L) return(NA_character_)
  for (sp in names(ENSEMBL_GENE_RE)) {
    if (mean(grepl(ENSEMBL_GENE_RE[[sp]], ids)) >= min_fraction) return(sp)
  }
  NA_character_
}

#' Do these ids look like Ensembl gene ids?
#'
#' Asked before mapping so the step is skipped for anything else --
#' running it on symbols or UniProt accessions would be a no-op, but a
#' no-op that costs a table load and reports a mapping rate of zero.
#'
#' @param ids Character vector.
#' @param min_fraction How many must match to call the whole set Ensembl.
#' @keywords internal
#' @noRd
looks_like_ensembl <- function(ids, min_fraction = 0.5) {
  !is.na(ensembl_species(ids, min_fraction))
}

#' Map Ensembl gene ids to gene symbols
#'
#' Human ids (`ENSG...`) are mapped to HGNC symbols from a table bundled
#' with the package (see [hgnc_map_provenance()]). Mouse ids
#' (`ENSMUSG...`) are mapped to MGI symbols through the ortholog table of
#' the babelgene package -- the one msigdbr uses to write its mouse gene
#' sets -- when that is installed; it covers the mouse genes with a human
#' counterpart, which are the genes a mouse gene set can contain.
#'
#' The version suffix is ignored (`ENSG00000141510.17` is TP53), and so
#' is GENCODE's `_PAR_Y` suffix on the chromosome-Y copy of a
#' pseudoautosomal gene, which names the same gene as its X copy.
#'
#' @param ids Character vector of Ensembl gene ids, with or without the
#'   version suffix.
#' @return Character vector the same length as `ids`; `NA` where the id
#'   has no symbol (no approved HGNC symbol for a human id; no human
#'   ortholog, or babelgene not installed, for a mouse id).
#' @export
#' @family enrich
#' @examples
#' map_ensembl_symbols(c("ENSG00000141510", "ENSG00000012048"))
map_ensembl_symbols <- function(ids) {
  if (is.null(ids)) return(character(0))
  if (is.factor(ids) || is.numeric(ids)) ids <- as.character(ids)
  if (!is.character(ids)) arg_stop("ids", "a character vector", ids)
  out <- rep(NA_character_, length(ids))
  if (length(ids) == 0L) return(out)
  base <- ensembl_base_id(ids)
  human <- !is.na(base) & startsWith(base, "ENSG")
  mouse <- !is.na(base) & startsWith(base, "ENSMUSG")
  m <- hgnc_ensembl_map()
  if (any(human) && !is.null(m)) {
    out[human] <- m$symbol[match(base[human], m$ensembl_gene_id)]
  }
  if (any(mouse)) out[mouse] <- map_mouse_ensembl(base[mouse])
  out
}

#' Fill in `feature_symbol` for an Ensembl-keyed feature table
#'
#' Only fills gaps. A file that carried its own symbol column already had
#' it picked up by `materialize_feature_annot()`, and overwriting that
#' would silently replace what the user supplied.
#'
#' @param feature_df Feature annotation, one row per feature.
#' @param feature_ids Matrix row ids, in matrix order.
#' @return `list(feature_df, note)`; `note` is `NULL` when nothing was
#'   mapped, and otherwise says how many of how many, because that ratio
#'   is what predicts whether enrichment will return anything.
#' @keywords internal
#' @noRd
attach_gene_symbols <- function(feature_df, feature_ids) {
  unchanged <- list(feature_df = feature_df, note = NULL)
  species <- ensembl_species(feature_ids)
  if (is.na(species)) return(unchanged)

  existing <- feature_df$feature_symbol
  if (!is.null(existing) && !all(is.na(existing))) return(unchanged)

  mouse <- identical(species, "mouse")
  if (mouse && !is_installed("babelgene")) {
    return(list(feature_df = feature_df, note = paste0(
      "Feature ids look like mouse Ensembl gene ids, but the table that maps ",
      "them to gene symbols (the R package babelgene) is not installed, so no ",
      "symbols were added. Enrichment will find nothing, because pathway ",
      "databases are keyed on symbols.")))
  }

  mapped <- map_ensembl_symbols(feature_ids)
  # The PAR_Y copy of a gene whose X copy is also a row gets no symbol of
  # its own: one gene, one symbol, or enrichment would count it twice.
  is_par <- grepl("_PAR_Y$", feature_ids)
  if (any(is_par)) {
    twin <- is_par & ensembl_base_id(feature_ids) %in% ensembl_base_id(feature_ids[!is_par])
    mapped[twin] <- NA_character_
  }
  n_ok <- sum(!is.na(mapped))
  if (n_ok == 0L) {
    return(list(
      feature_df = feature_df,
      note = paste0(
        "Feature ids look like ", if (mouse) "mouse ", "Ensembl gene ids but ",
        "none matched an ", if (mouse) "MGI" else "HGNC", " ",
        "symbol. Enrichment will find nothing, because pathway databases ",
        "are keyed on symbols.")
    ))
  }

  # Aligned by id rather than by position: materialize_feature_annot()
  # may have reordered or subset the rows.
  feature_df$feature_symbol <-
    mapped[match(feature_df$feature_id, feature_ids)]

  n_of <- c(format(n_ok, big.mark = ","), format(length(feature_ids), big.mark = ","))
  note <- if (mouse) {
    sprintf(
      paste0("Gene symbol: mapped from mouse Ensembl id -- %s of %s features ",
             "matched an MGI symbol. The rest, mostly genes without a human ",
             "counterpart, keep their id and are excluded from enrichment: the ",
             "mouse gene sets are carried over from human ones and could not ",
             "have matched them anyway. Choose mouse as the species in ",
             "Enrichment. Source: %s"),
      n_of[[1L]], n_of[[2L]], mgi_map_provenance())
  } else {
    sprintf(
      paste0("Gene symbol: mapped from Ensembl id -- %s of %s features ",
             "matched an HGNC symbol. The rest keep their id and are ",
             "excluded from enrichment, which no pathway database would ",
             "have matched anyway. Source: %s"),
      n_of[[1L]], n_of[[2L]], hgnc_map_provenance())
  }
  list(feature_df = feature_df, note = note)
}

# ---- GENCODE's chromosome-Y PAR copies --------------------------------------
#
# GENCODE (to release 43) annotates each gene of the pseudoautosomal
# regions twice, once on X and once on Y, and gives the Y copy the X
# copy's id with `_PAR_Y` appended. The two sequences are identical, so
# an aligner or quantifier puts the reads on one of them -- in practice
# the X copy, the Y region being masked in the usual analysis set -- and
# the Y rows come out all zero. GENCODE's own advice is to drop them.
#
# Dropped, then, when they hold nothing. A PAR_Y row that does hold
# values is kept as a feature of its own rather than merged into the X
# row: summing is right for counts but wrong for intensities or
# log-scale values, which cannot be told apart here, and keeping it
# loses nothing -- its id stays distinct, so nothing downstream sees a
# duplicate, and attach_gene_symbols() leaves its symbol empty while the
# X copy carries it.

#' Which rows to keep of a matrix whose ids may carry `_PAR_Y`
#'
#' @param mat Numeric matrix, features in rows, ids as row names.
#' @return `list(keep, notes)`: a logical vector over the rows, and the
#'   notes for the report (`character(0)` when no id carries the suffix).
#' @keywords internal
#' @noRd
par_y_rows <- function(mat) {
  ids <- rownames(mat)
  keep <- rep(TRUE, nrow(mat))
  none <- list(keep = keep, notes = character(0))
  if (is.null(ids)) return(none)
  par <- grepl("_PAR_Y$", ids)
  if (!any(par)) return(none)
  vals <- mat[par, , drop = FALSE]
  empty <- rowSums(!is.na(vals) & vals != 0) == 0L
  drop <- which(par)[empty]
  # Never every row: a table of nothing else is not GENCODE's layout.
  if (length(drop) == nrow(mat)) drop <- integer(0)
  keep[drop] <- FALSE
  notes <- character(0)
  if (length(drop)) {
    notes <- c(notes, sprintf(paste(
      "Dropped %d row(s) whose ids end in '_PAR_Y' (e.g. %s): GENCODE's copies",
      "of the pseudoautosomal genes on chromosome Y, which held no values --",
      "their reads are counted on the X-chromosome copy."),
      length(drop), ids[drop[[1L]]]))
  }
  kept <- which(par & keep)
  if (length(kept)) {
    notes <- c(notes, sprintf(paste(
      "Kept %d row(s) whose ids end in '_PAR_Y' (e.g. %s) because they hold",
      "values: each is the chromosome-Y copy of a pseudoautosomal gene, kept",
      "as a separate feature from its X-chromosome copy."),
      length(kept), ids[kept[[1L]]]))
  }
  list(keep = keep, notes = notes)
}

