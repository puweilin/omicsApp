#!/usr/bin/env Rscript
#
# Bake the MSigDB gene-set tables into the image.
#
# The first msigdbr() call in an R process spends ~10s lazy-loading the
# package's data. In a per-user-container deployment every user pays
# that before their first enrichment, and pays it again after every
# container recycle. Pre-building the tables here turns it into ~0.3s of
# disk read, once per collection, for everyone.
#
# omicsCore reads this directory when OMICSCORE_GENESET_CACHE points at
# it (see packages/omicsCore/R/enrich-genesets.R). The contract is
# small and deliberately so: one qs2 file per database/organism named
# `<database>__<organism>.qs2`, holding a data.frame with the columns
# below. Anything omicsCore cannot read it simply ignores, falling back
# to a live msigdbr fetch -- a stale or missing cache costs 10 seconds,
# never a wrong answer.
#
# This script depends only on msigdbr and qs2 so it can run before
# omicsCore is installed, which keeps its Docker layer cached across
# changes to the packages. The collection map below therefore mirrors
# DB_MSIGDBR_MAP in enrich-genesets.R. If the two drift, the missing
# database is simply not pre-warmed.

suppressPackageStartupMessages({
  library(msigdbr)
  library(qs2)
})

OUT_DIR <- Sys.getenv("OMICSCORE_GENESET_CACHE", "/opt/genesets")
ORGANISMS <- c("Homo sapiens", "Mus musculus")

# Mirrors omicsCore:::DB_MSIGDBR_MAP.
COLLECTIONS <- list(
  hallmark     = list(collection = "H",  subcollection = NA_character_),
  kegg         = list(collection = "C2", subcollection = "CP:KEGG_LEGACY"),
  reactome     = list(collection = "C2", subcollection = "CP:REACTOME"),
  wikipathways = list(collection = "C2", subcollection = "CP:WIKIPATHWAYS"),
  go_bp        = list(collection = "C5", subcollection = "GO:BP"),
  go_mf        = list(collection = "C5", subcollection = "GO:MF"),
  go_cc        = list(collection = "C5", subcollection = "GO:CC")
)

# Only the columns omicsCore reads. Keeping the rest would put ~4x the
# bytes resident in every container for no benefit.
KEEP <- c("gs_name", "gene_symbol", "gs_description", "gs_id")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cache_file <- function(database, organism) {
  file.path(OUT_DIR, sprintf("%s__%s.qs2", database,
                             gsub("[^A-Za-z0-9]+", "_", organism)))
}

fetch <- function(cfg, organism) {
  args <- names(formals(msigdbr::msigdbr))
  if ("collection" %in% args) {
    msigdbr::msigdbr(
      species = organism,
      collection = cfg$collection,
      subcollection = if (is.na(cfg$subcollection)) NULL else cfg$subcollection
    )
  } else {
    # msigdbr < 10 carries MSigDB 7.5, whose KEGG subcategory is
    # "CP:KEGG". MSigDB renamed it "CP:KEGG_LEGACY" on adding
    # KEGG_MEDICUS, and COLLECTIONS above uses the modern name -- so
    # without this the filter matches nothing and the kegg table is
    # silently empty rather than absent.
    sub <- cfg$subcollection
    if (identical(sub, "CP:KEGG_LEGACY")) sub <- "CP:KEGG"
    raw <- msigdbr::msigdbr(species = organism, category = cfg$collection)
    if (!is.na(sub) && "gs_subcat" %in% colnames(raw)) {
      raw <- raw[raw$gs_subcat == sub, , drop = FALSE]
    }
    raw
  }
}

total <- 0
failed <- character(0)

for (organism in ORGANISMS) {
  for (database in names(COLLECTIONS)) {
    path <- cache_file(database, organism)
    ok <- tryCatch({
      df <- fetch(COLLECTIONS[[database]], organism)
      cols <- intersect(KEEP, colnames(df))
      if (!"gs_name" %in% cols) stop("no gs_name column")
      qs2::qs_save(as.data.frame(df[, cols, drop = FALSE]), path)
      TRUE
    }, error = function(e) {
      message(sprintf("  ! %s / %s: %s", database, organism,
                      conditionMessage(e)))
      FALSE
    })
    if (!isTRUE(ok)) {
      failed <- c(failed, sprintf("%s/%s", database, organism))
      next
    }
    size <- file.size(path)
    total <- total + size
    message(sprintf("  %-14s %-14s %6.2f MB", database, organism,
                    size / 1024^2))
  }
}

message(sprintf("\nGene-set cache: %d files, %.1f MB total in %s",
                length(list.files(OUT_DIR)), total / 1024^2, OUT_DIR))

# A collection that cannot be built is not fatal -- omicsCore falls back
# to fetching it live. Fail the build only if nothing at all was cached,
# which means the prewarm is broken rather than one collection renamed.
if (length(list.files(OUT_DIR)) == 0L) {
  stop("Gene-set prewarm produced no files; the cache would be useless.")
}
if (length(failed) > 0L) {
  message("Not pre-warmed (will be fetched live on first use): ",
          paste(failed, collapse = ", "))
}
