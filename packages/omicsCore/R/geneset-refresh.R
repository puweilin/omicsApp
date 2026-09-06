# Live KEGG refresh for the on-disk gene-set cache (enrich-genesets.R).
#
# KEGG's license forbids redistributing pathway data in bulk, so no package
# can ship current KEGG gene sets -- the bundled MSigDB `CP:KEGG_LEGACY`
# collection is frozen at ~2011 (186 human pathways vs ~370 today). What the
# license does allow is per-query access to the public REST API. This file
# implements that route: three GET requests assemble a symbol-space gene-set
# table and write it to the same qs2 cache the prewarmed MSigDB tables use,
# so `database = "kegg"` serves current pathways wherever the cache has been
# refreshed and falls back to the legacy snapshot everywhere else.
#
# Design constraints, in order:
#   1. No surprise network calls. A live fetch happens only inside
#      refresh_geneset_cache(), or through the TTL auto-refresh -- which
#      arms itself only after a cache is already live-sourced (someone
#      opted in once) and can be disabled with OMICSCORE_GENESET_TTL_DAYS=0.
#   2. Never worse than the fallback. Every failure path keeps the previous
#      cache file (writes are staged to a temp file and renamed) or degrades
#      to the msigdbr snapshot.
#   3. No new dependencies. The endpoints are plain TSV read with base R;
#      no KEGGREST, no org.* annotation package.
#
# Refreshed tables carry a `gs_source` column ("KEGG REST <date>" or
# "MSigDB msigdbr <version>"); run_enrichment() records it in the bundle
# params so a result can always say which pathway definitions produced it.

KEGG_REST_BASE <- "https://rest.kegg.jp"

KEGG_ORG_CODES <- c("Homo sapiens" = "hsa", "Mus musculus" = "mmu")

kegg_org_code <- function(organism) {
  code <- KEGG_ORG_CODES[organism]
  if (is.na(code)) {
    stop("Live KEGG refresh supports: ",
         paste(names(KEGG_ORG_CODES), collapse = ", "), ".", call. = FALSE)
  }
  unname(code)
}

# Single seam for the network: every KEGG request goes through here so
# tests can mock one binding. `path` is the REST path, e.g. "/list/hsa".
kegg_rest_table <- function(path) {
  utils::read.delim(
    paste0(KEGG_REST_BASE, path),
    header = FALSE, sep = "\t", quote = "", comment.char = "",
    colClasses = "character", stringsAsFactors = FALSE
  )
}

# /list/<org> column 4 is "SYM1, SYM2; description" when the entry has gene
# symbols and a plain description (no semicolon) when it does not.
parse_kegg_symbol <- function(annotation) {
  has_symbols <- grepl(";", annotation, fixed = TRUE)
  first <- sub(";.*$", "", annotation)
  symbol <- trimws(sub(",.*$", "", first))
  ifelse(has_symbols & nzchar(symbol), symbol, NA_character_)
}

# Assemble a cache-contract table (gs_name / gene_symbol / gs_description /
# gs_id / gs_source) from the current KEGG release. Three requests:
# gene->pathway links, pathway names, and the organism gene list that maps
# KEGG gene entries to symbols. KEGG's own symbol annotation is used
# directly, which keeps org.* packages out of the picture; symbols may
# differ from the OrgDb mapping for a handful of aliased genes.
fetch_kegg_live <- function(organism) {
  org <- kegg_org_code(organism)

  link  <- kegg_rest_table(paste0("/link/pathway/", org))
  paths <- kegg_rest_table(paste0("/list/pathway/", org))
  genes <- kegg_rest_table(paste0("/list/", org))

  if (ncol(link) < 2L || ncol(paths) < 2L || ncol(genes) < 4L) {
    stop("Unexpected KEGG REST response shape.", call. = FALSE)
  }

  gene_id <- sub(paste0("^", org, ":"), "", link[[1L]])
  path_id <- sub("^path:", "", link[[2L]])

  path_name <- stats::setNames(paths[[2L]], sub("^path:", "", paths[[1L]]))
  # "Glycolysis / Gluconeogenesis - Homo sapiens (human)" -> the name.
  path_name <- trimws(sub(sprintf(" - %s \\([^)]*\\)$", organism), "", path_name))

  symbol <- stats::setNames(
    parse_kegg_symbol(genes[[4L]]),
    sub(paste0("^", org, ":"), "", genes[[1L]])
  )

  df <- data.frame(
    gs_name = unname(path_name[path_id]),
    gene_symbol = unname(symbol[gene_id]),
    gs_description = unname(path_name[path_id]),
    gs_id = path_id,
    stringsAsFactors = FALSE
  )
  df <- unique(df[!is.na(df$gs_name) & !is.na(df$gene_symbol), , drop = FALSE])
  if (nrow(df) == 0L) {
    stop("KEGG REST returned no usable pathway members.", call. = FALSE)
  }
  df$gs_source <- paste("KEGG REST", format(Sys.Date()))
  df
}

# msigdbr snapshot trimmed to the cache contract, for the non-KEGG
# databases (and as provenance-stamped storage of the legacy tables).
fetch_msigdbr_snapshot <- function(database, organism) {
  if (!is_installed("msigdbr")) {
    stop("Refreshing '", database, "' needs the msigdbr package.",
         call. = FALSE)
  }
  df <- fetch_msigdbr_raw(database, organism)
  keep <- intersect(c(GENESET_CACHE_COLUMNS, "human_gene_symbol"), colnames(df))
  df <- as.data.frame(df)[, keep, drop = FALSE]
  df$gs_source <- paste("MSigDB msigdbr",
                        as.character(utils::packageVersion("msigdbr")))
  df
}

geneset_cache_age_days <- function(path) {
  if (!file.exists(path)) return(NA_real_)
  as.numeric(difftime(Sys.time(), file.mtime(path), units = "days"))
}

# TTL for the auto-refresh of live-sourced KEGG caches. Unset -> 30 days;
# 0, negative, or non-numeric -> disabled.
geneset_ttl_days <- function() {
  raw <- Sys.getenv("OMICSCORE_GENESET_TTL_DAYS", "")
  if (!nzchar(raw)) return(30)
  val <- suppressWarnings(as.numeric(raw))
  if (is.na(val) || val <= 0) return(0)
  val
}

is_live_kegg_table <- function(df) {
  is.data.frame(df) && "gs_source" %in% colnames(df) &&
    any(startsWith(as.character(df$gs_source), "KEGG REST"))
}

# Called by fetch_msigdbr_table() when a kegg table comes off disk. Only a
# live-sourced cache past its TTL triggers a re-fetch; a prewarmed msigdbr
# table never touches the network. Failure returns the stale table -- stale
# costs freshness, never the run.
maybe_refresh_kegg_cache <- function(cached, organism) {
  if (!is_live_kegg_table(cached)) return(cached)
  ttl <- geneset_ttl_days()
  if (ttl <= 0) return(cached)
  dir <- geneset_cache_dir()
  if (is.null(dir)) return(cached)
  age <- geneset_cache_age_days(geneset_cache_file(dir, "kegg", organism))
  if (is.na(age) || age < ttl) return(cached)

  message("KEGG gene-set cache is ", round(age), " days old; refreshing from ",
          "KEGG REST (set OMICSCORE_GENESET_TTL_DAYS=0 to disable).")
  tryCatch(
    refresh_geneset_cache("kegg", organism, force = TRUE, quiet = TRUE),
    error = function(e) NULL
  )
  read_geneset_cache("kegg", organism) %||% cached
}

# Provenance of the table a run actually used, for the bundle params.
# The session memo is guaranteed populated by the run that just happened.
geneset_table_source <- function(database, organism) {
  df <- cache_get(paste0("msig::", database, "::", organism))
  if (is.data.frame(df) && "gs_source" %in% colnames(df)) {
    return(as.character(df$gs_source[[1L]]))
  }
  "MSigDB (msigdbr)"
}

#' Refresh the on-disk gene-set cache
#'
#' Rebuilds the qs2 gene-set cache that [run_enrichment()], [run_gsva()],
#' and [list_gene_sets()] read. For `"kegg"` the table is fetched live from
#' the KEGG REST API (current pathway definitions, in symbol space -- about
#' twice the pathways of the bundled 2011 `KEGG_LEGACY` snapshot); all other
#' databases are re-snapshotted from the installed `msigdbr`. Every table is
#' stamped with a `gs_source` so downstream results can name the pathway
#' definitions they used.
#'
#' The cache directory is `OMICSCORE_GENESET_CACHE` when set, otherwise a
#' per-user directory (`tools::R_user_dir("omicsCore", "cache")`). Once a
#' live KEGG table is in place it is kept current automatically: an
#' enrichment run that finds it older than `OMICSCORE_GENESET_TTL_DAYS`
#' (default 30) re-fetches it first, falling back to the stale table if the
#' network is unavailable. Set `OMICSCORE_GENESET_TTL_DAYS=0` to disable.
#'
#' KEGG-derived caches are for local use only; KEGG's license does not
#' permit redistributing them.
#'
#' @param databases Database keys to refresh. Defaults to all supported
#'   databases; only `"kegg"` involves the network.
#' @param organism Organism shorthand (e.g. `"Hs"`, `"Mm"`).
#' @param max_age_days Skip databases whose cache file is younger than this.
#' @param force If `TRUE`, refresh regardless of age.
#' @param quiet If `TRUE`, suppress progress messages.
#'
#' @return Invisibly, a `data.frame` with one row per database: `database`,
#'   `action` (`"refreshed"`, `"fresh"`, or `"failed"`), `n_sets`, `path`.
#' @export
#' @family enrich
#' @examples
#' \dontrun{
#'   refresh_geneset_cache("kegg")
#'   geneset_cache_status()
#' }
refresh_geneset_cache <- function(
  databases = SUPPORTED_ENRICH_DATABASES,
  organism = "Hs",
  max_age_days = 30,
  force = FALSE,
  quiet = FALSE
) {
  organism <- normalize_organism(organism)
  databases <- unique(vapply(databases, normalize_enrich_database, character(1L)))
  if (length(databases) == 0L) {
    stop("`databases` must name at least one database.", call. = FALSE)
  }
  dir <- geneset_cache_dir(write = TRUE)

  rows <- lapply(databases, function(db) {
    path <- geneset_cache_file(dir, db, organism)
    row <- function(action, n_sets = NA_integer_) {
      data.frame(database = db, action = action, n_sets = n_sets,
                 path = path, stringsAsFactors = FALSE)
    }

    age <- geneset_cache_age_days(path)
    if (!force && !is.na(age) && age < max_age_days) {
      if (!quiet) {
        message("'", db, "' cache is ", round(age, 1),
                " days old; keeping it (force = TRUE to refresh).")
      }
      return(row("fresh"))
    }

    df <- tryCatch(
      if (db == "kegg") fetch_kegg_live(organism)
      else fetch_msigdbr_snapshot(db, organism),
      error = function(e) {
        warning("Refresh of '", db, "' failed, keeping the previous cache: ",
                conditionMessage(e), call. = FALSE)
        NULL
      }
    )
    if (is.null(df)) return(row("failed"))

    # Stage-and-rename so an interrupted write can never clobber a good
    # cache with a truncated one.
    # Unique per call: the monthly cron and a manual refresh can
    # overlap, and a shared temp name lets one truncate the other.
    tmp <- tempfile(pattern = paste0(basename(path), "."),
                    tmpdir = dirname(path), fileext = ".tmp")
    qs2::qs_save(df, tmp)
    file.rename(tmp, path)
    cache_drop(paste0("msig::", db, "::", organism))

    n_sets <- length(unique(df$gs_name))
    if (!quiet) {
      message("Refreshed '", db, "': ", n_sets, " gene sets (",
              df$gs_source[[1L]], ").")
    }
    row("refreshed", n_sets)
  })

  invisible(do.call(rbind, rows))
}

#' Report the state of the on-disk gene-set cache
#'
#' One row per database: whether a cache file is present and readable, how
#' many gene sets it holds, where it came from (`gs_source`, e.g.
#' `"KEGG REST 2026-08-19"`), and how old it is. Useful as a pre-flight
#' check before deciding whether [refresh_geneset_cache()] is worth a
#' network round-trip.
#'
#' @param databases Database keys to report on. Defaults to all supported.
#' @param organism Organism shorthand (e.g. `"Hs"`, `"Mm"`).
#'
#' @return A `data.frame` with columns `database`, `cached`, `n_sets`,
#'   `source`, `cached_at`, `age_days`.
#' @export
#' @family enrich
geneset_cache_status <- function(
  databases = SUPPORTED_ENRICH_DATABASES,
  organism = "Hs"
) {
  assert_names(databases, "databases")
  organism <- normalize_organism(organism)
  databases <- unique(vapply(databases, normalize_enrich_database, character(1L)))
  dir <- geneset_cache_dir()

  rows <- lapply(databases, function(db) {
    df <- if (!is.null(dir)) read_geneset_cache(db, organism) else NULL
    path <- if (!is.null(dir)) geneset_cache_file(dir, db, organism) else ""
    on_disk <- nzchar(path) && file.exists(path)
    data.frame(
      database = db,
      cached = !is.null(df),
      n_sets = if (!is.null(df)) length(unique(df$gs_name)) else NA_integer_,
      source = if (is.null(df)) NA_character_
               else if ("gs_source" %in% colnames(df)) as.character(df$gs_source[[1L]])
               else "prewarm (msigdbr)",
      cached_at = if (on_disk) format(file.mtime(path)) else NA_character_,
      age_days = if (on_disk) round(geneset_cache_age_days(path), 1) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
