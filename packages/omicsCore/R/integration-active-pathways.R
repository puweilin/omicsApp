# ActivePathways combined-p pathway enrichment across two omics layers.
# Heavy Suggests-gated; requires both `ActivePathways` and our enrichment
# stack (msigdbr → MSigDB gene sets). Builds a per-feature scores matrix
# (rows = gene symbols, columns = experiment tags) from two diff_bundles,
# then runs ActivePathways against an in-memory GMT-like list.

ensure_active_pathways <- function() {
  if (!requireNamespace("ActivePathways", quietly = TRUE)) {
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

  p_col <- if (p_preference == "adjusted") "adj_p_value" else "p_value"

  build_score <- function(df) {
    sub <- df[, c(by, p_col), drop = FALSE]
    names(sub) <- c("key", "p")
    sub <- sub[!is.na(sub$key) & nzchar(sub$key), , drop = FALSE]
    sub <- sub[!duplicated(sub$key), , drop = FALSE]
    sub
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
  scores[s_a$key, tag_a] <- pmax(s_a$p, .Machine$double.xmin)
  scores[s_b$key, tag_b] <- pmax(s_b$p, .Machine$double.xmin)

  gene_sets <- get_gene_set_list(
    database = database,
    organism = organism,
    min_size = geneset_filter[[1L]],
    max_size = geneset_filter[[2L]]
  )
  if (length(gene_sets) == 0L) {
    stop("No gene sets available for ActivePathways after size filter.")
  }
  gmt <- to_active_pathways_gmt(gene_sets)

  ap_raw <- ActivePathways::ActivePathways(
    scores = scores,
    gmt = gmt,
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
  shared <- !is.na(evidence_str) & vapply(strsplit(evidence_str, ","),
                                          function(x) length(x) >= 2L, logical(1))

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
    direction = ifelse(shared, "shared", "unique"),
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
