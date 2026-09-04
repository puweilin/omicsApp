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

cache_drop <- function(key) {
  if (exists(key, envir = .omicsCore_enrich_cache, inherits = FALSE)) {
    rm(list = key, envir = .omicsCore_enrich_cache)
  }
  invisible(NULL)
}

# ---- on-disk cache ------------------------------------------------------
#
# The first `msigdbr()` call in a process pays ~10s to lazy-load the
# package's data, which every user of a fresh container would otherwise
# pay before their first enrichment. Pointing `OMICSCORE_GENESET_CACHE`
# at a directory of pre-built tables removes that: the deploy image
# bakes them in at build time (see deploy/docker/prewarm_genesets.R).
#
# The cached tables carry only the columns this file actually reads.
# Dropping the rest takes GO:BP from 106 MB to 23 MB in memory and
# 2.4 MB on disk, which matters when several collections sit resident
# in every one of a handful of per-user containers. Tables written by
# refresh_geneset_cache() additionally carry `gs_source` (provenance).
GENESET_CACHE_COLUMNS <- c("gs_name", "gene_symbol",
                           "gs_description", "gs_id")

# OMICSCORE_GENESET_CACHE (explicit, wins when set) or a per-user cache
# directory that refresh_geneset_cache() populates. Read mode returns NULL
# when the resolved directory does not exist -- an env var pointing at a
# missing directory is respected as "no cache", not silently redirected.
geneset_cache_dir <- function(write = FALSE) {
  dir <- Sys.getenv("OMICSCORE_GENESET_CACHE", "")
  if (!nzchar(dir)) dir <- tools::R_user_dir("omicsCore", which = "cache")
  if (write) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    return(dir)
  }
  if (!dir.exists(dir)) return(NULL)
  dir
}

#' Path of the cache file for one database / organism pair
#'
#' @param dir Cache directory.
#' @param database Normalised database key.
#' @param organism Normalised organism name.
#'
#' @return A file path.
#' @keywords internal
#' @noRd
geneset_cache_file <- function(dir, database, organism) {
  file.path(dir, sprintf("%s__%s.qs2", database,
                         gsub("[^A-Za-z0-9]+", "_", organism)))
}

# Fails soft in every direction: a missing directory, an absent file, no
# qs2, a truncated write, or a table from an older column contract all
# just mean "no cache", and the caller falls back to msigdbr. A stale
# cache costs 10 seconds, never a wrong answer.
read_geneset_cache <- function(database, organism) {
  dir <- geneset_cache_dir()
  if (is.null(dir)) return(NULL)
  if (!is_installed("qs2")) return(NULL)
  path <- geneset_cache_file(dir, database, organism)
  if (!file.exists(path)) return(NULL)
  df <- tryCatch(qs2::qs_read(path), error = function(e) NULL)
  if (!is.data.frame(df) || nrow(df) == 0L) return(NULL)
  has_gene <- any(c("gene_symbol", "human_gene_symbol") %in% colnames(df))
  if (!("gs_name" %in% colnames(df)) || !has_gene) return(NULL)
  df
}

# ---- msigdbr fetch ------------------------------------------------------

fetch_msigdbr_table <- function(database, organism) {
  database <- normalize_enrich_database(database)
  organism <- normalize_organism(organism)
  cache_key <- paste0("msig::", database, "::", organism)
  hit <- cache_get(cache_key)
  if (!is.null(hit)) return(hit)

  cached <- read_geneset_cache(database, organism)
  if (!is.null(cached)) {
    if (identical(database, "kegg")) {
      cached <- maybe_refresh_kegg_cache(cached, organism)
    }
    cache_set(cache_key, cached)
    return(cached)
  }

  df <- fetch_msigdbr_raw(database, organism)
  cache_set(cache_key, df)
  df
}

# Translate a DB_MSIGDBR_MAP subcollection to the name the pre-v10
# msigdbr data uses. Only KEGG differs: MSigDB 7.5 calls it "CP:KEGG",
# and renamed it "CP:KEGG_LEGACY" when KEGG_MEDICUS arrived. The map
# carries the modern name, so on old msigdbr the filter would match no
# rows and hand back an empty table instead of the KEGG gene sets --
# an empty result, not an error, which is the worse of the two.
legacy_subcollection <- function(subcollection) {
  if (identical(subcollection, "CP:KEGG_LEGACY")) "CP:KEGG" else subcollection
}

fetch_msigdbr_raw <- function(database, organism) {
  cfg <- DB_MSIGDBR_MAP[[database]]
  msig_args <- names(formals(msigdbr::msigdbr))

  if ("collection" %in% msig_args) {
    msigdbr::msigdbr(
      species = organism,
      collection = cfg$collection,
      subcollection = if (is.na(cfg$subcollection)) NULL else cfg$subcollection
    )
  } else {
    # Legacy msigdbr (<v10) API.
    sub <- legacy_subcollection(cfg$subcollection)
    raw <- msigdbr::msigdbr(species = organism, category = cfg$collection)
    if (!is.na(sub) && "gs_subcat" %in% colnames(raw)) {
      raw <- raw[raw$gs_subcat == sub, , drop = FALSE]
    }
    raw
  }
}

# Build long-form (term, gene) and (term, name) tables in symbol space.
# Returns a list with `term2gene` and `term2name` data.frames, suitable for
# `clusterProfiler::enricher()` / `GSEA()`.
build_term_tables <- function(database, organism) {
  df <- fetch_msigdbr_table(database, organism)
  gene_col <- if ("gene_symbol" %in% colnames(df)) "gene_symbol" else "human_gene_symbol"
  # The display name comes from gs_name, not gs_description.
  #
  # For Hallmark a description is a short sentence, which is why this
  # went unnoticed. For GO and Reactome it is the full definition with
  # its citations -- "The chemical reactions and pathways involving
  # 10-formyltetrahydrofolate, the formylated derivative of
  # tetrahydrofolate. [GOC:ai]" -- and that was arriving as the label on
  # a dotplot axis.
  name_col <- NA_character_

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
    terms <- unique(t2g$term)
    t2n <- data.frame(
      term = terms,
      name = prettify_gene_set_name(terms),
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

# MSigDB set names carry their collection as a prefix and join words
# with underscores: HALLMARK_TNFA_SIGNALING_VIA_NFKB. Neither belongs on
# a plot axis.
GENE_SET_NAME_PREFIXES <- c(
  "HALLMARK_", "GOBP_", "GOMF_", "GOCC_", "KEGG_MEDICUS_", "KEGG_LEGACY_",
  "KEGG_", "REACTOME_", "WP_", "BIOCARTA_", "PID_"
)

# Turn a set name into a label.
#
# Casing is left as MSigDB wrote it. Lowercasing reads better on a long
# GO term -- "10 formyltetrahydrofolate metabolic process" -- but every
# rule for deciding which tokens to spare turns some gene symbol into a
# word, and TNFA rendered as "Tnfa" is wrong in a way that ALL CAPS is
# only ugly. A label that is hard to read can still be looked up; one
# that is wrong cannot.
prettify_gene_set_name <- function(x) {
  x <- as.character(x)
  for (p in GENE_SET_NAME_PREFIXES) {
    hit <- startsWith(x, p) & !is.na(x)
    x[hit] <- substring(x[hit], nchar(p) + 1L)
  }
  trimws(gsub("_", " ", x, fixed = TRUE))
}
