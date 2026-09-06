# GSEA backend. Uses `clusterProfiler::GSEA()` against an MSigDB-derived
# TERM2GENE table in symbol space (same path as ORA — see enrich-ora.R).
# The public dispatcher is `run_enrichment()` in run-enrichment.R.

run_gsea_database <- function(
  ranked_features,
  database,
  organism = "Hs",
  p_cutoff = 0.05,
  output_p_cutoff = NULL,
  p_adjust_method = "BH",
  min_size = 10L,
  max_size = 500L,
  seed = TRUE
) {
  ensure_enrichment_deps()
  database <- normalize_enrich_database(database)
  organism <- normalize_organism(organism)

  ranked_features <- ranked_features[!is.na(ranked_features)]
  if (length(ranked_features) == 0L) return(NULL)
  ranked_features <- sort(ranked_features, decreasing = TRUE)
  ranked_features <- ranked_features[!duplicated(names(ranked_features))]

  terms <- build_term_tables(database = database, organism = organism)

  # The permutations draw from the random stream. Pinned so that two
  # runs on one ranking agree, and put back so that the caller's own
  # stream is where they left it -- it used to be advanced by every
  # GSEA run, and in a fresh session a global seed was planted.
  with_fixed_seed(if (isTRUE(seed)) 123L else NULL, tryCatch(
    clusterProfiler::GSEA(
      geneList = ranked_features,
      TERM2GENE = terms$term2gene,
      TERM2NAME = terms$term2name,
      pvalueCutoff = p_cutoff,
      pAdjustMethod = p_adjust_method,
      minGSSize = min_size,
      maxGSSize = max_size,
      seed = isTRUE(seed)
    ),
    error = function(e) {
      warning("GSEA failed for database '", database, "': ", conditionMessage(e),
              call. = FALSE)
      NULL
    }
  ))
}

# Build a ranked vector from a diff bundle and run GSEA against a database.
# `direction` filters the standardized output but does not affect the rank
# vector; clusterProfiler always returns both signs.
run_gsea_from_bundle <- function(
  diff_bundle,
  database,
  organism = "Hs",
  direction = c("both", "up", "down"),
  p_cutoff = 0.05,
  # GSEA never selects features -- the whole ranked list goes in -- so
  # here p_cutoff only ever bounded the output, and this is its name.
  output_p_cutoff = NULL,
  p_adjust_method = "BH",
  min_size = 10L,
  max_size = 500L
) {
  direction <- match.arg(direction)

  result_df <- diff_result_from_bundle(diff_bundle)
  feature_col <- if ("feature_symbol" %in% colnames(result_df)) "feature_symbol" else "feature_id"
  ranked <- make_ranked_features(
    result_df = result_df,
    feature_col = feature_col,
    rank_col = "effect"
  )
  comparison <- diff_bundle$params$comparison %||% "comparison"

  obj <- run_gsea_database(
    ranked_features = ranked,
    database = database,
    organism = organism,
    p_cutoff = output_p_cutoff %||% p_cutoff,
    p_adjust_method = p_adjust_method,
    min_size = min_size,
    max_size = max_size
  )

  std <- standardize_enrich_result(
    enrich_obj = obj,
    database = database,
    result_type = "gsea",
    comparison = comparison,
    source_label = paste0("gsea_", database)
  )

  if (direction != "both" && nrow(std) > 0L && "direction" %in% colnames(std)) {
    std <- std[!is.na(std$direction) & std$direction == direction, , drop = FALSE]
    rownames(std) <- NULL
  }

  list(object = obj, std = std)
}
