# Per-feature correlation across paired samples between two omics layers.
# For each gene/protein pair that appears in both experiments, computes a
# Pearson or Spearman correlation across the donors that have data in
# both layers (after `build_sample_pairs()` alignment). Returns the full
# integration schema with `effect = r`, `statistic_type = "spearman"` /
# `"pearson"`, and p-values from `stats::cor.test()`.

run_integration_correlation <- function(
  project,
  experiments,
  method = c("spearman", "pearson"),
  by = "feature_symbol",
  p_adjust_method = "BH",
  min_samples = 4L,
  p_cutoff = 0.05
) {
  experiments <- resolve_experiment_pair(project, experiments)
  method <- match.arg(method)
  tag_a <- experiments[[1L]]
  tag_b <- experiments[[2L]]

  sample_pairs <- build_sample_pairs(project, tag_a, tag_b)
  if (nrow(sample_pairs) < min_samples) {
    stop("Need at least ", min_samples, " paired samples for correlation; got ",
         nrow(sample_pairs), ".")
  }

  feature_pairs <- build_feature_pairs(project, tag_a, tag_b, by = by)

  mat_a <- project$experiments[[tag_a]]$expr_mat[feature_pairs$feature_a,
                                                  sample_pairs[[tag_a]], drop = FALSE]
  mat_b <- project$experiments[[tag_b]]$expr_mat[feature_pairs$feature_b,
                                                  sample_pairs[[tag_b]], drop = FALSE]

  # Pre-log RNA-seq raw counts to put them on a comparable scale to the
  # other layer. We respect each input's assay_type independently.
  mat_a <- coerce_to_continuous(mat_a, project$experiments[[tag_a]]$assay_type)
  mat_b <- coerce_to_continuous(mat_b, project$experiments[[tag_b]]$assay_type)

  n_feat <- nrow(feature_pairs)
  effect <- statistic <- p_value <- rep(NA_real_, n_feat)
  for (i in seq_len(n_feat)) {
    x <- as.numeric(mat_a[i, ])
    y <- as.numeric(mat_b[i, ])
    keep <- is.finite(x) & is.finite(y)
    if (sum(keep) < min_samples) next
    res <- tryCatch(
      suppressWarnings(stats::cor.test(x[keep], y[keep], method = method)),
      error = function(e) NULL
    )
    if (is.null(res)) next
    effect[i] <- unname(res$estimate)
    statistic[i] <- unname(res$statistic)
    p_value[i] <- res$p.value
  }
  adj <- stats::p.adjust(p_value, method = p_adjust_method)

  direction <- rep(NA_character_, n_feat)
  direction[!is.na(effect) & effect >= 0] <- "positive"
  direction[!is.na(effect) & effect <  0] <- "negative"

  out <- new_integration_result_template()
  out <- out[rep(1L, n_feat), , drop = FALSE]  # keep zero-row structure when n_feat == 0
  if (n_feat == 0L) return(out[FALSE, , drop = FALSE])

  out <- data.frame(
    feature_id = feature_pairs$feature_id,
    feature_symbol = feature_pairs$feature_symbol,
    result_type = "correlation",
    experiments = paste(tag_a, "vs", tag_b),
    comparison = NA_character_,
    effect = effect,
    effect_type = paste0(method, "_r"),
    statistic = statistic,
    statistic_type = method,
    p_value = p_value,
    adj_p_value = adj,
    direction = direction,
    quadrant = NA_character_,
    is_significant = !is.na(adj) & adj < p_cutoff,
    source_label = paste0("integration_correlation_", tag_a, "_", tag_b),
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL

  list(
    std = out,
    info = list(
      method = method,
      experiments = experiments,
      n_features = n_feat,
      n_samples = nrow(sample_pairs)
    )
  )
}

coerce_to_continuous <- function(mat, assay_type) {
  if (identical(assay_type, "raw_count")) {
    log2(mat + 1)
  } else {
    mat
  }
}
