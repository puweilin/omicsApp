# ActivePathways combined-p pathway enrichment across two omics layers.
# Heavy Suggests-gated; requires both `ActivePathways` and our enrichment
# stack (msigdbr → MSigDB gene sets). Builds a per-feature scores matrix
# (rows = gene symbols, columns = experiment tags) from two diff_bundles,
# then runs ActivePathways against an in-memory GMT-like list.

ensure_active_pathways <- function() {
  if (!is_installed("ActivePathways")) {
    stop(
      "Package 'ActivePathways' is required for method = 'active_pathways'. ",
      "Install with: install.packages('ActivePathways').",
      call. = FALSE
    )
  }
}

# Convert a named list of character vectors (the format returned by
# `get_gene_set_list()`) into the `GMT` structure that
# `ActivePathways::ActivePathways()` expects:
#   list( id = ..., name = ..., genes = c(...) ), with class "GMT".
to_active_pathways_gmt <- function(gene_sets) {
  ids <- names(gene_sets)
  out <- lapply(seq_along(gene_sets), function(i) {
    list(id = ids[[i]], name = ids[[i]], genes = as.character(gene_sets[[i]]))
  })
  names(out) <- ids
  class(out) <- "GMT"
  out
}

run_integration_active_pathways <- function(
  project,
  experiments,
  diff_bundles,
  database = "hallmark",
  organism = "Hs",
  by = "feature_symbol",
  p_preference = c("raw", "adjusted"),
  significant = 0.05,
  geneset_filter = c(5L, 1000L),
  merge_method = "Brown"
) {
  experiments <- resolve_experiment_pair(project, experiments)
  p_preference <- match.arg(p_preference)
  validate_diff_bundles(diff_bundles, experiments)
  ensure_active_pathways()
  ensure_enrichment_deps()

  tag_a <- experiments[[1L]]
  tag_b <- experiments[[2L]]
  res_a <- diff_bundles[[tag_a]]$results$diff_result_df
  res_b <- diff_bundles[[tag_b]]$results$diff_result_df
  check_diff_result_schema(res_a)
  check_diff_result_schema(res_b)
  if (!by %in% colnames(res_a) || !by %in% colnames(res_b)) {
    stop("`", by, "` must be a column in both diff_result_df's.")
  }

  for (side in list(list(res_a, tag_a), list(res_b, tag_b))) {
    n_cmp <- length(unique(stats::na.omit(side[[1L]]$comparison)))
    if (n_cmp > 1L) {
      stop("The diff result for '", side[[2L]], "' holds ", n_cmp,
           " comparisons; pick one with `select_comparison()` before ",
           "integrating.", call. = FALSE)
    }
  }

  p_col <- if (p_preference == "adjusted") "adj_p_value" else "p_value"

  # Same key rules as concordance: case-insensitive symbols, one row per
  # gene (the most abundant), so the two methods integrate the same genes.
  build_score <- function(df) {
    key <- integration_join_key(df[[by]])
    keep <- !is.na(key)
    df <- df[keep, , drop = FALSE]
    key <- key[keep]
    idx <- dedupe_by_key(key, df$base_mean)
    data.frame(key = key[idx], p = as.numeric(df[[p_col]][idx]),
               stringsAsFactors = FALSE)
  }
  s_a <- build_score(res_a)
  s_b <- build_score(res_b)
  all_keys <- union(s_a$key, s_b$key)
  if (length(all_keys) == 0L) {
    stop("No features found in either diff bundle.")
  }

  scores <- matrix(
    1, nrow = length(all_keys), ncol = 2L,
    dimnames = list(all_keys, c(tag_a, tag_b))
  )
  # A gene that was not tested in a layer (no p-value) carries no
  # evidence there, the same as one that was not measured.
  na_to_1 <- function(p) ifelse(is.na(p), 1, pmax(p, .Machine$double.xmin))
  scores[s_a$key, tag_a] <- na_to_1(s_a$p)
  scores[s_b$key, tag_b] <- na_to_1(s_b$p)

  gene_sets <- get_gene_set_list(
    database = database,
    organism = organism,
    min_size = geneset_filter[[1L]],
    max_size = geneset_filter[[2L]]
  )
  if (length(gene_sets) == 0L) {
    stop("No gene sets available for ActivePathways after size filter.")
  }
  # The gene sets in the same key space as the scores. The scores are
  # upper-cased (integration_join_key()), the GMT was not, so a mouse
  # gene set ("Trp53") matched nothing.
  gene_sets <- lapply(gene_sets, function(g) unique(stats::na.omit(integration_join_key(g))))
  gmt <- to_active_pathways_gmt(gene_sets)
  # The background is what was measured, not every gene in the GMT:
  # ActivePathways' default background is the GMT's genes, which counts
  # unmeasured genes as tested-and-null and inflates every pathway.
  background <- intersect(rownames(scores), unique(unlist(gene_sets, use.names = FALSE)))
  if (length(background) < 2L) {
    stop("Too few measured genes are in the '", database, "' gene sets for ",
         "ActivePathways. Check the organism and that features carry gene symbols.",
         call. = FALSE)
  }
  scores <- scores[background, , drop = FALSE]

  ap_raw <- ActivePathways::ActivePathways(
    scores = scores,
    gmt = gmt,
    background = background,
    geneset_filter = geneset_filter,
    significant = significant,
    merge_method = merge_method,
    cytoscape_file_tag = NA
  )

  if (is.null(ap_raw) || nrow(ap_raw) == 0L) {
    out <- new_integration_result_template()
    return(list(
      std = out,
      raw = NULL,
      info = list(
        experiments = experiments,
        database = normalize_enrich_database(database),
        organism = normalize_organism(organism),
        n_features = nrow(scores),
        n_pathways_significant = 0L
      )
    ))
  }

  ap_df <- as.data.frame(ap_raw, stringsAsFactors = FALSE)
  # ActivePathways returns columns: term_id, term_name, adjusted_p_val,
  # term_size, overlap, evidence (list), Genes_<col> (list-cols).
  n <- nrow(ap_df)
  evidence_str <- vapply(ap_df$evidence, function(x) {
    if (is.null(x) || length(x) == 0L) return(NA_character_)
    paste(as.character(unlist(x)), collapse = ",")
  }, character(1))
  # `evidence` names the layers whose own p-values found the pathway, or
  # "combined" when only the merged p-value did. Two layers is "shared";
  # one layer is "unique"; "combined" alone is its own class -- the
  # pathway exists only because the layers were merged, which is the
  # finding ActivePathways is for, and it used to be filed as "unique".
  ev_sets <- strsplit(ifelse(is.na(evidence_str), "", evidence_str), ",")
  evidence_class <- vapply(ev_sets, function(x) {
    x <- trimws(x)
    x <- x[nzchar(x)]
    if (!length(x)) return(NA_character_)
    layers <- intersect(x, c(tag_a, tag_b))
    if (length(layers) >= 2L) "shared"
    else if (length(layers) == 1L) "unique"
    else if ("combined" %in% x) "combined"
    else NA_character_
  }, character(1))

  out <- data.frame(
    feature_id = as.character(ap_df$term_id),
    feature_symbol = as.character(ap_df$term_name),
    result_type = "active_pathways",
    experiments = paste(tag_a, "vs", tag_b),
    comparison = paste(
      diff_bundles[[tag_a]]$params$comparison %||% "comparison",
      diff_bundles[[tag_b]]$params$comparison %||% "comparison",
      sep = " | "
    ),
    effect = -log10(pmax(as.numeric(ap_df$adjusted_p_val), .Machine$double.xmin)),
    effect_type = "neg_log10_padj",
    statistic = as.numeric(ap_df$adjusted_p_val),
    statistic_type = "adjusted_p_val",
    p_value = as.numeric(ap_df$adjusted_p_val),
    adj_p_value = as.numeric(ap_df$adjusted_p_val),
    direction = evidence_class,
    quadrant = evidence_str,
    is_significant = !is.na(ap_df$adjusted_p_val) & ap_df$adjusted_p_val < significant,
    source_label = paste0("integration_active_pathways_", tag_a, "_", tag_b),
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL

  list(
    std = out,
    raw = ap_df,
    info = list(
      experiments = experiments,
      database = normalize_enrich_database(database),
      organism = normalize_organism(organism),
      n_features = nrow(scores),
      n_pathways_significant = sum(out$is_significant, na.rm = TRUE),
      merge_method = merge_method
    )
  )
}
