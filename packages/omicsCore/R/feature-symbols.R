# Ensembl gene id -> HGNC symbol.
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

# `ENSG00000141510.17` and `ENSG00000141510` are the same gene; the
# suffix is the annotation version. Vendors differ on whether they keep
# it, and a table keyed on the unversioned id matches neither if we do
# not strip it.
strip_ensembl_version <- function(x) sub("\\.\\d+$", "", as.character(x))

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
  ids <- ids[!is.na(ids)]
  if (length(ids) == 0L) return(FALSE)
  mean(grepl("^ENSG\\d{6,}(\\.\\d+)?$", ids)) >= min_fraction
}

#' Map Ensembl gene ids to HGNC symbols
#'
#' @param ids Character vector of Ensembl gene ids, with or without the
#'   version suffix.
#' @return Character vector the same length as `ids`; `NA` where the id
#'   has no approved HGNC symbol.
#' @export
#' @family enrich
#' @examples
#' map_ensembl_symbols(c("ENSG00000141510", "ENSG00000012048"))
map_ensembl_symbols <- function(ids) {
  out <- rep(NA_character_, length(ids))
  m <- hgnc_ensembl_map()
  if (is.null(m) || length(ids) == 0L) return(out)
  idx <- match(strip_ensembl_version(ids), m$ensembl_gene_id)
  out[!is.na(idx)] <- m$symbol[idx[!is.na(idx)]]
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
attach_hgnc_symbols <- function(feature_df, feature_ids) {
  unchanged <- list(feature_df = feature_df, note = NULL)
  if (!looks_like_ensembl(feature_ids)) return(unchanged)

  existing <- feature_df$feature_symbol
  if (!is.null(existing) && !all(is.na(existing))) return(unchanged)

  mapped <- map_ensembl_symbols(feature_ids)
  n_ok <- sum(!is.na(mapped))
  if (n_ok == 0L) {
    return(list(
      feature_df = feature_df,
      note = paste0(
        "Feature ids look like Ensembl gene ids but none matched an HGNC ",
        "symbol. Enrichment will find nothing, because pathway databases ",
        "are keyed on symbols.")
    ))
  }

  # Aligned by id rather than by position: materialize_feature_annot()
  # may have reordered or subset the rows.
  feature_df$feature_symbol <-
    mapped[match(feature_df$feature_id, feature_ids)]

  list(
    feature_df = feature_df,
    note = sprintf(
      paste0("Gene symbol: mapped from Ensembl id -- %s of %s features ",
             "matched an HGNC symbol. The rest keep their id and are ",
             "excluded from enrichment, which no pathway database would ",
             "have matched anyway. Source: %s"),
      format(n_ok, big.mark = ","),
      format(length(feature_ids), big.mark = ","),
      hgnc_map_provenance())
  )
}
