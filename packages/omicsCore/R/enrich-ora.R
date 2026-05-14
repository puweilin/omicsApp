# ORA backend. Always uses `clusterProfiler::enricher()` against an MSigDB-
# derived TERM2GENE table in symbol space so we don't need org.Hs.eg.db or
# ReactomePA. The public dispatcher is `run_enrichment()` (run-enrichment.R).

run_ora_database <- function(
  features,
  universe,
  database,
  organism = "Hs",
  p_cutoff = 0.05,
  p_adjust_method = "BH"
) {
  ensure_enrichment_deps()
  database <- normalize_enrich_database(database)
  organism <- normalize_organism(organism)

  features <- unique(stats::na.omit(features))
  universe <- unique(stats::na.omit(universe))
  if (length(features) == 0L) return(NULL)

  terms <- build_term_tables(database = database, organism = organism)

  tryCatch(
    clusterProfiler::enricher(
      gene = features,
      universe = universe,
      TERM2GENE = terms$term2gene,
      TERM2NAME = terms$term2name,
      pvalueCutoff = p_cutoff,
      pAdjustMethod = p_adjust_method,
      qvalueCutoff = 1
    ),
    error = function(e) {
      warning("ORA failed for database '", database, "': ", conditionMessage(e),
              call. = FALSE)
      NULL
    }
  )
}

# Direction-aware ORA driver used by `run_enrichment()`. Splits the input
# diff bundle by sign and runs ORA per direction (or as a single set when
# `direction = "both"`).
run_ora_from_bundle <- function(
  diff_bundle,
  database,
  organism = "Hs",
  direction = c("both", "up", "down"),
  p_cutoff = 0.05,
  effect_cutoff = NULL,
  p_preference = c("adjusted", "raw"),
  p_adjust_method = "BH"
) {
  direction <- match.arg(direction)
  p_preference <- match.arg(p_preference)

  result_df <- diff_result_from_bundle(diff_bundle)
  universe <- unique(stats::na.omit(result_df$feature_symbol))
  if (length(universe) == 0L) {
    universe <- unique(stats::na.omit(result_df$feature_id))
  }

  sig_df <- filter_diff_results(
    result_df = result_df,
    p_cutoff = p_cutoff,
    p_preference = p_preference,
    effect_cutoff = effect_cutoff
  )
  feature_col <- if ("feature_symbol" %in% colnames(sig_df)) "feature_symbol" else "feature_id"
  comparison <- diff_bundle$params$comparison %||% "comparison"

  case_definitions <- switch(
    direction,
    both = list(both = sig_df),
    up   = list(up = sig_df[sig_df$effect > 0, , drop = FALSE]),
    down = list(down = sig_df[sig_df$effect < 0, , drop = FALSE])
  )

  per_direction <- lapply(names(case_definitions), function(dir_label) {
    sub <- case_definitions[[dir_label]]
    feats <- unique(stats::na.omit(sub[[feature_col]]))
    obj <- run_ora_database(
      features = feats,
      universe = universe,
      database = database,
      organism = organism,
      p_cutoff = p_cutoff,
      p_adjust_method = p_adjust_method
    )
    std <- standardize_enrich_result(
      enrich_obj = obj,
      database = database,
      result_type = "ora",
      comparison = comparison,
      source_label = paste0("ora_", dir_label, "_", database),
      default_direction = if (dir_label == "both") NA_character_ else dir_label
    )
    list(object = obj, std = std)
  })
  names(per_direction) <- names(case_definitions)

  list(
    objects = lapply(per_direction, `[[`, "object"),
    std = dplyr::bind_rows(lapply(per_direction, `[[`, "std"))
  )
}
