# Gene-set provider for the enrichment layer. All databases are fetched from
# msigdbr in **symbol space** so that ORA / GSEA can be run via the unified
# `clusterProfiler::enricher()` / `GSEA()` interface without any OrgDb
# dependency. Results are memoised in a private package env for the duration
# of the R session (no disk writes — see `omicsApp/docs/export-manifest.md`).

# Manifest-spec database keys → msigdbr (collection, subcollection).
# All keys are lowercase to match the export-manifest contract.
DB_MSIGDBR_MAP <- list(
  hallmark     = list(collection = "H",  subcollection = NA_character_),
  kegg         = list(collection = "C2", subcollection = "CP:KEGG_LEGACY"),
  reactome     = list(collection = "C2", subcollection = "CP:REACTOME"),
  wikipathways = list(collection = "C2", subcollection = "CP:WIKIPATHWAYS"),
  go_bp        = list(collection = "C5", subcollection = "GO:BP"),
  go_mf        = list(collection = "C5", subcollection = "GO:MF"),
  go_cc        = list(collection = "C5", subcollection = "GO:CC")
)

SUPPORTED_ENRICH_DATABASES <- names(DB_MSIGDBR_MAP)

# Common organism shorthands.
ORGANISM_ALIASES <- c(
  Hs                = "Homo sapiens",
  human             = "Homo sapiens",
  `Homo sapiens`    = "Homo sapiens",
  Mm                = "Mus musculus",
  mouse             = "Mus musculus",
  `Mus musculus`    = "Mus musculus"
)

normalize_organism <- function(organism) {
  if (length(organism) != 1L || is.na(organism)) {
    stop("`organism` must be a single non-missing string.")
  }
  key <- as.character(organism)
  if (!key %in% names(ORGANISM_ALIASES)) {
    stop("Unsupported organism: ", organism,
         ". Supported: ", paste(unique(ORGANISM_ALIASES), collapse = ", "), ".")
  }
  unname(ORGANISM_ALIASES[key])
}

normalize_enrich_database <- function(database) {
  if (length(database) != 1L || is.na(database)) {
    stop("`database` must be a single non-missing string.")
  }
  key <- tolower(as.character(database))
  alias_map <- c(
    h = "hallmark", hallmarks = "hallmark",
    gobp = "go_bp", gomf = "go_mf", gocc = "go_cc",
    wiki = "wikipathways"
  )
  if (key %in% names(alias_map)) key <- unname(alias_map[key])
  if (!key %in% SUPPORTED_ENRICH_DATABASES) {
    stop("Unsupported enrichment database: ", database,
         ". Supported: ", paste(SUPPORTED_ENRICH_DATABASES, collapse = ", "), ".")
  }
  key
}

ensure_enrichment_deps <- function(extras = character(0)) {
  needed <- c("clusterProfiler", "msigdbr", extras)
  missing <- needed[!vapply(needed, is_installed, logical(1))]
  if (length(missing) > 0L) {
    stop(
      "The enrichment layer needs the following packages: ",
      paste(missing, collapse = ", "),
      ". Install with: omicsCore::install_optional('enrichment').",
      call. = FALSE
    )
  }
}

# ---- session memoisation ------------------------------------------------

.omicsCore_enrich_cache <- new.env(parent = emptyenv())

cache_get <- function(key) {
  if (!exists(key, envir = .omicsCore_enrich_cache, inherits = FALSE)) {
    return(NULL)
  }
  get(key, envir = .omicsCore_enrich_cache, inherits = FALSE)
}

cache_set <- function(key, value) {
  assign(key, value, envir = .omicsCore_enrich_cache)
  invisible(value)
}

# ---- msigdbr fetch ------------------------------------------------------

fetch_msigdbr_table <- function(database, organism) {
  database <- normalize_enrich_database(database)
  organism <- normalize_organism(organism)
  cache_key <- paste0("msig::", database, "::", organism)
  hit <- cache_get(cache_key)
  if (!is.null(hit)) return(hit)

  cfg <- DB_MSIGDBR_MAP[[database]]
  msig_args <- names(formals(msigdbr::msigdbr))

  df <- if ("collection" %in% msig_args) {
    msigdbr::msigdbr(
      species = organism,
      collection = cfg$collection,
      subcollection = if (is.na(cfg$subcollection)) NULL else cfg$subcollection
    )
  } else {
    # Legacy msigdbr (<v10) API.
    raw <- msigdbr::msigdbr(species = organism, category = cfg$collection)
    if (!is.na(cfg$subcollection) && "gs_subcat" %in% colnames(raw)) {
      raw <- raw[raw$gs_subcat == cfg$subcollection, , drop = FALSE]
    }
    raw
  }

  cache_set(cache_key, df)
  df
}

# Build long-form (term, gene) and (term, name) tables in symbol space.
# Returns a list with `term2gene` and `term2name` data.frames, suitable for
# `clusterProfiler::enricher()` / `GSEA()`.
build_term_tables <- function(database, organism) {
  df <- fetch_msigdbr_table(database, organism)
  gene_col <- if ("gene_symbol" %in% colnames(df)) "gene_symbol" else "human_gene_symbol"
  name_col <- if ("gs_description" %in% colnames(df)) "gs_description" else NA_character_

  keep <- !is.na(df[[gene_col]]) & !is.na(df$gs_name)
  t2g <- data.frame(
    term = df$gs_name[keep],
    gene = df[[gene_col]][keep],
    stringsAsFactors = FALSE
  )
  t2g <- unique(t2g)

  if (!is.na(name_col)) {
    name_keep <- !is.na(df$gs_name)
    t2n <- data.frame(
      term = df$gs_name[name_keep],
      name = dplyr::coalesce(df[[name_col]][name_keep], df$gs_name[name_keep]),
      stringsAsFactors = FALSE
    )
    t2n <- unique(t2n)
  } else {
    t2n <- data.frame(
      term = unique(t2g$term),
      name = unique(t2g$term),
      stringsAsFactors = FALSE
    )
  }

  list(term2gene = t2g, term2name = t2n)
}

# Build a named list of character vectors: pathway_name -> gene symbols.
# Used by GSVA.
get_gene_set_list <- function(database, organism = "Hs",
                              min_size = 10L, max_size = 500L) {
  ensure_enrichment_deps()
  df <- fetch_msigdbr_table(database, organism)
  gene_col <- if ("gene_symbol" %in% colnames(df)) "gene_symbol" else "human_gene_symbol"
  keep <- !is.na(df[[gene_col]]) & !is.na(df$gs_name)
  gs_list <- split(df[[gene_col]][keep], df$gs_name[keep])
  gs_list <- lapply(gs_list, unique)
  lens <- lengths(gs_list)
  gs_list[lens >= min_size & lens <= max_size]
}

# ---- public ------------------------------------------------------------

#' List the gene sets available in a pathway database
#'
#' Returns a `tibble` with one row per gene-set member, sourced from MSigDB
#' (`msigdbr`). All gene identifiers are returned as HGNC / official gene
#' symbols. Results are cached for the duration of the R session so repeated
#' calls are cheap.
#'
#' Supported databases (lowercase): `hallmark`, `kegg`, `reactome`,
#' `wikipathways`, `go_bp`, `go_mf`, `go_cc`.
#'
#' @param database One of the supported database keys (see above).
#' @param organism Organism, accepting `"Hs"`, `"human"`, `"Homo sapiens"`,
#'   `"Mm"`, `"mouse"`, or `"Mus musculus"`.
#'
#' @return A `tibble` with columns `database`, `pathway_id`, `pathway_name`,
#'   `gene_symbol`.
#' @export
#' @family enrich
#' @examples
#' \dontrun{
#'   list_gene_sets("hallmark")
#'   list_gene_sets("go_bp", organism = "Hs")
#' }
list_gene_sets <- function(database, organism = "Hs") {
  ensure_enrichment_deps()
  database <- normalize_enrich_database(database)
  organism <- normalize_organism(organism)

  df <- fetch_msigdbr_table(database, organism)
  gene_col <- if ("gene_symbol" %in% colnames(df)) "gene_symbol" else "human_gene_symbol"
  id_col <- if ("gs_id" %in% colnames(df)) "gs_id" else "gs_name"
  desc_col <- if ("gs_description" %in% colnames(df)) "gs_description" else "gs_name"

  out <- tibble::tibble(
    database = database,
    pathway_id = as.character(df[[id_col]]),
    pathway_name = as.character(df[[desc_col]]),
    gene_symbol = as.character(df[[gene_col]])
  )
  out <- out[!is.na(out$gene_symbol) & !is.na(out$pathway_name), , drop = FALSE]
  out <- out[!duplicated(out[, c("pathway_id", "gene_symbol")]), , drop = FALSE]
  out
}
