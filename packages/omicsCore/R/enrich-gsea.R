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
  # One value per gene: the strongest, whichever its sign. Keeping the
  # first after a decreasing sort kept the most positive, so a gene
  # measured twice leaned "up".
  ranked_features <- ranked_features[order(abs(ranked_features), decreasing = TRUE)]
  ranked_features <- ranked_features[!duplicated(names(ranked_features))]
  ranked_features <- sort(ranked_features, decreasing = TRUE)

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
  ranked <- gsea_rank_vector(result_df, feature_col)
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

# The ranking GSEA walks: the test statistic with the sign of the
# effect, not the raw effect. A fold change of 3 from a gene seen in two
# samples out-ranked a fold change of 1.5 seen cleanly in every sample,
# and genes with no p-value (not tested, filtered) still took part.
# The metric used is kept as the vector's "metric" attribute.
gsea_rank_vector <- function(result_df, feature_col) {
  st <- unique(stats::na.omit(result_df$statistic_type))
  stat <- result_df$statistic
  eff <- result_df$effect
  if (length(st) == 1L && st %in% c("t", "wald") && any(is.finite(stat))) {
    metric <- sign(eff) * abs(stat)
    label <- "signed test statistic"
  } else if (length(st) == 1L && st == "F" && any(is.finite(stat))) {
    # A one-degree-of-freedom F (edgeR's QL test) is a squared t.
    metric <- sign(eff) * sqrt(abs(stat))
    label <- "signed sqrt(F)"
  } else {
    metric <- sign(eff) * -log10(pmax(result_df$p_value, .Machine$double.xmin))
    label <- "sign(effect) * -log10(p)"
  }
  keep <- !is.na(result_df$p_value) & is.finite(metric) &
    !is.na(result_df[[feature_col]]) & nzchar(result_df[[feature_col]])
  out <- metric[keep]
  names(out) <- result_df[[feature_col]][keep]
  out <- out[order(abs(out), decreasing = TRUE)]
  out <- out[!duplicated(names(out))]
  out <- sort(out, decreasing = TRUE)
  attr(out, "metric") <- label
  out
}

gsea_rank_metric <- function(diff_bundle) {
  df <- diff_result_from_bundle(diff_bundle)
  col <- if ("feature_symbol" %in% colnames(df)) "feature_symbol" else "feature_id"
  attr(gsea_rank_vector(df, col), "metric")
}
