# Per-feature correlation across paired samples between two omics layers.
# For each protein-gene pair the two experiments share (by symbol or by a
# feature link; a gene measured by two proteins gives two pairs, each
# correlated on its own and each a row of the result), computes a
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
  p_cutoff = 0.05,
  link = NULL
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

  feature_pairs <- build_feature_pairs(project, tag_a, tag_b, by = by, link = link)

  mat_a <- project$experiments[[tag_a]]$expr_mat[feature_pairs$feature_a,
                                                  sample_pairs[[tag_a]], drop = FALSE]
  mat_b <- project$experiments[[tag_b]]$expr_mat[feature_pairs$feature_b,
                                                  sample_pairs[[tag_b]], drop = FALSE]

  # Pre-log RNA-seq raw counts to put them on a comparable scale to the
  # other layer. We respect each input's assay_type independently.
  # Library sizes are taken over the whole layer, not the matched subset
  # of features, or a count layer's CPM would depend on which genes the
  # other layer happened to measure.
  mat_a <- coerce_to_continuous(mat_a, project$experiments[[tag_a]]$assay_type,
                                lib_size = lib_sizes(project$experiments[[tag_a]],
                                                     sample_pairs[[tag_a]]))
  mat_b <- coerce_to_continuous(mat_b, project$experiments[[tag_b]]$assay_type,
                                lib_size = lib_sizes(project$experiments[[tag_b]],
                                                     sample_pairs[[tag_b]]))

  n_feat <- nrow(feature_pairs)
  effect <- statistic <- p_value <- rep(NA_real_, n_feat)
  for (i in seq_len(n_feat)) {
    x <- as.numeric(mat_a[i, ])
    y <- as.numeric(mat_b[i, ])
    keep <- is.finite(x) & is.finite(y)
    if (sum(keep) < min_samples) next
    # exact = FALSE for Spearman: above nine pairs cor.test()'s "exact"
    # p-value is an Edgeworth-series approximation that truncates to 0
    # for strong correlations -- twelve samples at rho = 0.96 came back
    # as p = 0, i.e. -log10(p) = 308 on the plot. The t approximation
    # is what cor.test() already uses whenever there are ties.
    res <- tryCatch(
      suppressWarnings(
        if (method == "spearman") {
          stats::cor.test(x[keep], y[keep], method = method, exact = FALSE)
        } else {
          stats::cor.test(x[keep], y[keep], method = method)
        }),
      error = function(e) NULL
    )
    if (is.null(res)) next
    effect[i] <- unname(res$estimate)
    statistic[i] <- unname(res$statistic)
    p_value[i] <- res$p.value
    # Identical rankings (rho = 1) give t = Inf and p = 0. Over n pairs
    # the smallest p a rank correlation can reach is 2 / n!, so that is
    # the floor.
    if (method == "spearman") {
      p_value[i] <- max(p_value[i], 2 / factorial(sum(keep)))
    }
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
    # Beyond the schema: the two features each row correlates.
    feature_id_a = feature_pairs$feature_a,
    feature_id_b = feature_pairs$feature_b,
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL

  list(
    std = out,
    info = list(
      method = method,
      experiments = experiments,
      n_features = n_feat,
      n_samples = nrow(sample_pairs),
      pairing_source = attr(sample_pairs, "source") %||% NA_character_,
      n_ambiguous_samples = attr(sample_pairs, "n_ambiguous") %||% 0L,
      feature_pairing = attr(feature_pairs, "info")
    )
  )
}

# Put a layer on a log scale before correlating it. Counts become
# log2(CPM + 1) -- without the library-size step a deeply sequenced sample
# is "high" in every gene, and that shared offset shows up as a positive
# correlation for everything. Other linear assays (raw intensities, TPM,
# FPKM) are log2(x + 1). Log-scale assays, and labels outside the
# vocabulary, are left as they are.
coerce_to_continuous <- function(mat, assay_type, lib_size = NULL) {
  if (identical(assay_type, "raw_count")) {
    if (is.null(lib_size)) lib_size <- colSums(mat, na.rm = TRUE)
    lib_size[!is.finite(lib_size) | lib_size <= 0] <- NA_real_
    log2(sweep(mat, 2L, lib_size / 1e6, "/") + 1)
  } else if (isTRUE(assay_type %in% c("raw_intensity", "tpm", "fpkm"))) {
    log2(pmax(mat, 0) + 1)
  } else {
    mat
  }
}

lib_sizes <- function(input, samples) {
  if (!identical(input$assay_type, "raw_count")) return(NULL)
  colSums(input$expr_mat[, samples, drop = FALSE], na.rm = TRUE)
}
