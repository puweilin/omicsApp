# ORA backend. Always uses `clusterProfiler::enricher()` against an MSigDB-
# derived TERM2GENE table in symbol space so we don't need org.Hs.eg.db or
# ReactomePA. The public dispatcher is `run_enrichment()` (run-enrichment.R).

run_ora_database <- function(
  features,
  universe,
  database,
  organism = "Hs",
  p_cutoff = 0.05,
  p_adjust_method = "BH",
  min_size = 10L,
  max_size = 500L
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
      qvalueCutoff = 1,
      minGSSize = min_size,
      maxGSSize = max_size
    ),
    error = function(e) {
      warning("ORA failed for database '", database, "': ", conditionMessage(e),
              call. = FALSE)
      NULL
    }
  )
}

# Direction-aware ORA driver used by `run_enrichment()`. Splits the input
# diff bundle by sign and runs ORA per direction: `"separate"` runs the
# up- and down-regulated lists one after the other, each with its own
# multiple-testing correction, exactly as two runs with `"up"` and
# `"down"` would; `"both"` pools the two into a single list.
run_ora_from_bundle <- function(
  diff_bundle,
  database,
  organism = "Hs",
  direction = c("separate", "up", "down", "both"),
  p_cutoff = 0.05,
  # Separate from p_cutoff, which for ORA also decides which features go
  # in. Bounding the stored result at the same number means a bundle can
  # only ever hold what one threshold admitted, so switching the display
  # between raw and adjusted p has nothing to switch to. NULL keeps the
  # old single-threshold behaviour.
  output_p_cutoff = NULL,
  effect_cutoff = NULL,
  p_preference = c("adjusted", "raw"),
  p_adjust_method = "BH",
  min_size = 10L,
  max_size = 500L
) {
  direction <- match.arg(direction)
  p_preference <- match.arg(p_preference)

  result_df <- diff_result_from_bundle(diff_bundle)
  # The universe is what was tested. A gene edgeR set aside as too low
  # to test, or one with no p-value, could never have been a hit, and
  # counting it in the background makes every set look enriched.
  tested <- result_df[!is.na(result_df$p_value), , drop = FALSE]
  if (nrow(tested) == 0L) tested <- result_df
  universe <- unique(stats::na.omit(tested$feature_symbol))
  if (length(universe) == 0L) {
    universe <- unique(stats::na.omit(tested$feature_id))
  }

  sig_df <- filter_diff_results(
    result_df = result_df,
    p_cutoff = p_cutoff,
    p_preference = p_preference,
    effect_cutoff = effect_cutoff
  )
  feature_col <- if ("feature_symbol" %in% colnames(sig_df)) "feature_symbol" else "feature_id"
  comparison <- diff_bundle$params$comparison %||% "comparison"

  up_df <- sig_df[!is.na(sig_df$effect) & sig_df$effect > 0, , drop = FALSE]
  down_df <- sig_df[!is.na(sig_df$effect) & sig_df$effect < 0, , drop = FALSE]
  case_definitions <- switch(
    direction,
    separate = list(up = up_df, down = down_df),
    both = list(both = sig_df),
    up   = list(up = up_df),
    down = list(down = down_df)
  )

  per_direction <- lapply(names(case_definitions), function(dir_label) {
    sub <- case_definitions[[dir_label]]
    feats <- unique(stats::na.omit(sub[[feature_col]]))
    obj <- run_ora_database(
      features = feats,
      universe = universe,
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
