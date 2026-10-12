# ActivePathways combined-p pathway enrichment across two omics layers.
# Heavy Suggests-gated; requires both `ActivePathways` and our enrichment
# stack (msigdbr → MSigDB gene sets). Builds a per-gene scores matrix
# (rows = gene symbols, columns = experiment tags) from two diff_bundles,
# then runs ActivePathways against an in-memory GMT-like list.
#
# Direction. Brown's method merges two p-values with no regard for sign,
# so a gene up in the RNA and down in the protein counted as much as one
# up in both, and a pathway came back "significant" with no word on which
# way it went. ActivePathways 2.0 added the directional merge (DPM,
# Slobodyanyuk et al. 2024): each gene's merged score is penalised when
# its layers' fold-change signs break a stated expectation. The
# expectation used here is that the two layers agree -- RNA and protein
# of a gene move together -- so the constraints vector is c(1, 1).
# Genes measured in one layer only carry no sign in the other (effect 0,
# p = 1), and neither help nor hurt.
#
# DPM penalises at the gene level; it says nothing of a pathway's
# direction. That is read off afterwards from the genes that drove the
# pathway (ActivePathways' `overlap`): each layer's direction is the sign
# most of those genes share in it, and the pathway's direction is "up"
# or "down" when the layers agree and "mixed" when they do not.

# The merge methods that take fold-change signs.
DIRECTIONAL_MERGE_METHODS <- c("DPM", "Fisher_directional",
                               "Stouffer_directional", "Strube_directional")

# Whether the installed ActivePathways can merge directionally (>= 2.0).
active_pathways_directional <- function() {
  "scores_direction" %in% names(formals(ActivePathways::ActivePathways))
}

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
  merge_method = "DPM",
  constraints_vector = c(1, 1),
  link = NULL
) {
  experiments <- resolve_experiment_pair(project, experiments)
  p_preference <- match.arg(p_preference)
  validate_diff_bundles(diff_bundles, experiments)
  ensure_active_pathways()
  ensure_enrichment_deps()
  assert_string(merge_method, "merge_method")
  if (!is.numeric(constraints_vector) || length(constraints_vector) != 2L ||
      !all(constraints_vector %in% c(-1, 0, 1))) {
    stop("`constraints_vector` must be two values from -1, 0, 1: the ",
         "expected sign relation of the two layers (c(1, 1): they agree).",
         call. = FALSE)
  }
  notes <- character(0)
  directional <- merge_method %in% DIRECTIONAL_MERGE_METHODS
  if (directional && !active_pathways_directional()) {
    notes <- c(notes, sprintf(paste(
      "The installed ActivePathways (%s) cannot weigh the direction of",
      "change; pathways were found with Brown's method, which ignores it.",
      "ActivePathways 2.0 or later adds the directional method."),
      as.character(utils::packageVersion("ActivePathways"))))
    merge_method <- "Brown"
    directional <- FALSE
  }

  tag_a <- experiments[[1L]]
  tag_b <- experiments[[2L]]
  res_a <- diff_bundles[[tag_a]]$results$diff_result_df
  res_b <- diff_bundles[[tag_b]]$results$diff_result_df
  check_diff_result_schema(res_a)
  check_diff_result_schema(res_b)
  if (is.null(link) && (!by %in% colnames(res_a) || !by %in% colnames(res_b))) {
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

  # Genes, not features, are what ActivePathways scores, so each layer is
  # brought to one row per gene. A feature paired across the layers (by
  # symbol or by the feature link, as concordance pairs them) stands for
  # its pair's gene; one measured in a single layer for its own symbol.
  pairs <- link_features(res_a$feature_id, res_a[[by]],
                         res_b$feature_id, res_b[[by]], link = link)
  gene_keys <- function(df, paired_rows) {
    own <- if (by %in% names(df)) integration_join_key(df[[by]])
           else rep(NA_character_, nrow(df))
    alone <- setdiff(seq_len(nrow(df)), paired_rows)
    keys <- rbind(
      data.frame(i = paired_rows, key = integration_join_key(pairs$label),
                 stringsAsFactors = FALSE),
      data.frame(i = alone, key = own[alone], stringsAsFactors = FALSE))
    unique(keys[!is.na(keys$key), , drop = FALSE])
  }
  # Where several features of a layer stand for one gene (isoforms,
  # several protein groups), the gene takes its most significant
  # feature's p-value and fold change. The minimum of k p-values is not a
  # p-value -- it is small by chance far more often than one feature's --
  # so it is corrected for the number of features on offer (Sidak:
  # 1 - (1 - p_min)^k). Before isoforms were kept, the most abundant
  # feature stood for the gene; that wasted the evidence of an isoform
  # that moved while the abundant one did not. A gene with one feature is
  # scored exactly as before.
  build_score <- function(df, paired_rows) {
    keys <- gene_keys(df, paired_rows)
    p <- as.numeric(df[[p_col]][keys$i])
    eff <- as.numeric(df$effect[keys$i])
    k <- stats::ave(as.numeric(!is.na(p)), keys$key, FUN = sum)
    ord <- order(keys$key, is.na(p), p, keys$i)
    best <- ord[!duplicated(keys$key[ord])]
    # In the order the genes first appear in the table, as before.
    best <- best[order(stats::ave(keys$i, keys$key, FUN = min)[best])]
    kk <- pmax(k[best], 1)
    pb <- p[best]
    data.frame(key = keys$key[best],
               p = ifelse(kk > 1, -expm1(kk * log1p(-pb)), pb),
               effect = eff[best], n_features = kk,
               stringsAsFactors = FALSE)
  }
  s_a <- build_score(res_a, pairs$i_a)
  s_b <- build_score(res_b, pairs$i_b)
  for (sc in list(list(s_a, tag_a), list(s_b, tag_b))) {
    n_multi <- sum(sc[[1L]]$n_features > 1)
    if (n_multi > 0L) {
      notes <- c(notes, sprintf(paste(
        "%d gene(s) are measured by more than one feature in '%s'. Each",
        "enters ActivePathways with its most significant feature's p-value,",
        "corrected for the number of its features, and that feature's fold change."),
        n_multi, sc[[2L]]))
    }
  }
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
  # Each gene's fold change per layer; 0 where a layer did not measure
  # or test it, which DPM reads as "no direction here".
  effects <- matrix(0, nrow = length(all_keys), ncol = 2L,
                    dimnames = list(all_keys, c(tag_a, tag_b)))
  na_to_0 <- function(x) ifelse(is.finite(x), x, 0)
  effects[s_a$key, tag_a] <- na_to_0(s_a$effect)
  effects[s_b$key, tag_b] <- na_to_0(s_b$effect)
  effects[scores == 1] <- 0

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
  effects <- effects[background, , drop = FALSE]

  ap_args <- list(
    scores = scores,
    gmt = gmt,
    background = background,
    geneset_filter = geneset_filter,
    significant = significant,
    merge_method = merge_method,
    cytoscape_file_tag = NA
  )
  if (directional) {
    # Only the sign enters DPM; a column ActivePathways is told carries
    # no direction (constraint 0) must hold zeros.
    dir_mat <- sign(effects)
    dir_mat[, constraints_vector == 0] <- 0
    ap_args$scores_direction <- dir_mat
    ap_args$constraints_vector <- as.numeric(constraints_vector)
  }
  ap_raw <- do.call(ActivePathways::ActivePathways, ap_args)

  info_base <- list(
    experiments = experiments,
    database = normalize_enrich_database(database),
    organism = normalize_organism(organism),
    n_features = nrow(scores),
    merge_method = merge_method,
    directional = directional,
    constraints_vector = if (directional) as.numeric(constraints_vector),
    notes = notes,
    feature_pairing = attr(pairs, "info"),
    n_genes_multi = c(sum(s_a$n_features > 1), sum(s_b$n_features > 1))
  )

  if (is.null(ap_raw) || nrow(ap_raw) == 0L) {
    out <- new_integration_result_template()
    return(list(
      std = out,
      raw = NULL,
      info = c(info_base, list(n_pathways_significant = 0L))
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
    direction = NA_character_,
    quadrant = evidence_str,
    is_significant = !is.na(ap_df$adjusted_p_val) & ap_df$adjusted_p_val < significant,
    source_label = paste0("integration_active_pathways_", tag_a, "_", tag_b),
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  # Which way each pathway went, from the genes that drove it.
  pdir <- pathway_directions(ap_df$overlap, effects, tag_a, tag_b)
  out$direction <- pdir$direction
  out$direction_a <- pdir$direction_a
  out$direction_b <- pdir$direction_b
  out$layers_agree <- pdir$layers_agree
  out$n_genes_agree <- pdir$n_agree
  out$n_genes_disagree <- pdir$n_disagree
  # Which layers found the pathway on their own: "shared" (both),
  # "unique" (one), "combined" (only the merged evidence). Before the
  # directional rework this was the `direction` column.
  out$evidence <- evidence_class

  list(
    std = out,
    raw = ap_df,
    info = c(info_base, list(
      n_pathways_significant = sum(out$is_significant, na.rm = TRUE)
    ))
  )
}

# Direction of each pathway from its driving genes (`overlap`, a list of
# gene vectors) and the per-layer fold changes. A layer's direction is
# "up" or "down" when more than half of the driving genes it measured
# moved that way, "mixed" otherwise, NA when it measured none of them.
# The pathway is "up"/"down" when the layers that say something agree,
# and "mixed" when they point opposite ways or either is itself mixed.
pathway_directions <- function(overlap, effects, tag_a, tag_b) {
  layer_dir <- function(e) {
    e <- e[e != 0]
    if (!length(e)) return(NA_character_)
    up <- mean(e > 0)
    if (up > 0.5) "up" else if (up < 0.5) "down" else "mixed"
  }
  rows <- lapply(overlap, function(genes) {
    genes <- intersect(as.character(unlist(genes)), rownames(effects))
    ea <- effects[genes, tag_a]
    eb <- effects[genes, tag_b]
    da <- layer_dir(ea)
    db <- layer_dir(eb)
    both <- ea != 0 & eb != 0
    known <- stats::na.omit(c(da, db))
    dir <- if (!length(known)) NA_character_
           else if (length(unique(known)) == 1L && known[[1L]] != "mixed") known[[1L]]
           else "mixed"
    list(direction = dir, direction_a = da, direction_b = db,
         layers_agree = if (is.na(da) || is.na(db)) NA
                        else da == db && da != "mixed",
         n_agree = sum(both & sign(ea) == sign(eb)),
         n_disagree = sum(both & sign(ea) != sign(eb)))
  })
  data.frame(
    direction = vapply(rows, `[[`, "", "direction"),
    direction_a = vapply(rows, `[[`, "", "direction_a"),
    direction_b = vapply(rows, `[[`, "", "direction_b"),
    layers_agree = vapply(rows, `[[`, NA, "layers_agree"),
    n_agree = vapply(rows, `[[`, 0L, "n_agree"),
    n_disagree = vapply(rows, `[[`, 0L, "n_disagree"),
    stringsAsFactors = FALSE
  )
}
