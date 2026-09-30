# Per-feature concordance between two diff_bundles. For each gene/protein
# that appears in both layers, joins the standardized effect, direction,
# and (adjusted) p-value, classifies the (dir_a, dir_b) pair into a
# quadrant, and computes a combined-p (Fisher's method).

run_integration_concordance <- function(
  project,
  experiments,
  diff_bundles,
  by = "feature_symbol",
  p_preference = c("adjusted", "raw"),
  p_cutoff = 0.05,
  p_adjust_method = "BH",
  effect_cutoff = 0
) {
  experiments <- resolve_experiment_pair(project, experiments)
  p_preference <- match.arg(p_preference)
  validate_diff_bundles(diff_bundles, experiments)
  if (is.null(effect_cutoff)) effect_cutoff <- 0

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

  build_side <- function(df, side) {
    key <- integration_join_key(df[[by]])
    keep <- !is.na(key)
    df <- df[keep, , drop = FALSE]
    key <- key[keep]
    idx <- dedupe_by_key(key, df$base_mean)
    df <- df[idx, , drop = FALSE]
    out <- data.frame(
      key = key[idx],
      symbol = as.character(df[[by]]),
      feature_id = as.character(df$feature_id),
      effect = as.numeric(df$effect),
      direction = direction_sign(df$direction, df$effect),
      p = as.numeric(df$p_value),
      padj = as.numeric(df$adj_p_value),
      stringsAsFactors = FALSE
    )
    names(out)[-1L] <- paste0(names(out)[-1L], "_", side)
    out
  }

  joined <- merge(build_side(res_a, "a"), build_side(res_b, "b"), by = "key")
  if (nrow(joined) == 0L) {
    return(list(
      std = new_integration_result_template(),
      info = list(experiments = experiments, n_features = 0L,
                  p_preference = p_preference)
    ))
  }
  joined <- joined[order(joined$key), , drop = FALSE]

  quadrant <- classify_concordance_quadrant(joined$direction_a, joined$direction_b)
  direction <- ifelse(
    is.na(quadrant), NA_character_,
    ifelse(quadrant %in% c("up_up", "down_down"), "concordant", "discordant")
  )

  # Effect = difference of (signed) effects. The two effects themselves
  # are kept as `effect_a` / `effect_b`: the effect-pair plot needs them,
  # and the difference alone cannot give them back.
  effect <- joined$effect_a - joined$effect_b

  # Combined p via Fisher's method on the *raw* p-values. Adjusted
  # p-values are not uniform under the null, so a Fisher statistic built
  # from them is not chi-squared and its "p" means nothing; the
  # combination is corrected once, below, across features. A p-value that
  # underflowed to 0 (DESeq2 does this for strong hits) is clamped rather
  # than dropped -- dropping it lost exactly the strongest features.
  pa <- joined$p_a
  pb <- joined$p_b
  ok <- !is.na(pa) & !is.na(pb)
  tiny <- .Machine$double.xmin
  chisq <- rep(NA_real_, nrow(joined))
  chisq[ok] <- -2 * (log(pmax(pa[ok], tiny)) + log(pmax(pb[ok], tiny)))
  combined_p <- stats::pchisq(chisq, df = 4L, lower.tail = FALSE)
  combined_adj <- stats::p.adjust(combined_p, method = p_adjust_method)

  # Per-layer significance, at the same thresholds the Differential view
  # reads its hit table at.
  sel_a <- if (p_preference == "adjusted") joined$padj_a else pa
  sel_b <- if (p_preference == "adjusted") joined$padj_b else pb
  sig_a <- !is.na(sel_a) & sel_a < p_cutoff &
    !is.na(joined$effect_a) & abs(joined$effect_a) >= effect_cutoff
  sig_b <- !is.na(sel_b) & sel_b < p_cutoff &
    !is.na(joined$effect_b) & abs(joined$effect_b) >= effect_cutoff

  is_sig <- sig_a & sig_b & !is.na(direction) & direction == "concordant"

  out <- data.frame(
    feature_id = joined$symbol_a,
    feature_symbol = joined$symbol_a,
    result_type = "concordance",
    experiments = paste(tag_a, "vs", tag_b),
    comparison = paste(
      diff_bundles[[tag_a]]$params$comparison %||% "comparison",
      diff_bundles[[tag_b]]$params$comparison %||% "comparison",
      sep = " | "
    ),
    effect = effect,
    effect_type = "effect_diff",
    statistic = chisq,
    statistic_type = "fisher_chisq",
    p_value = combined_p,
    adj_p_value = combined_adj,
    direction = direction,
    quadrant = quadrant,
    is_significant = is_sig,
    source_label = paste0("integration_concordance_", tag_a, "_", tag_b),
    # Beyond the schema: what each layer said on its own.
    feature_id_a = joined$feature_id_a,
    feature_id_b = joined$feature_id_b,
    effect_a = joined$effect_a,
    effect_b = joined$effect_b,
    p_value_a = pa,
    p_value_b = pb,
    adj_p_value_a = joined$padj_a,
    adj_p_value_b = joined$padj_b,
    significant_a = sig_a,
    significant_b = sig_b,
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL

  both <- sig_a & sig_b
  list(
    std = out,
    info = list(
      experiments = experiments,
      n_features = nrow(out),
      p_preference = p_preference,
      effect_cutoff = effect_cutoff,
      quadrant_counts = as.list(table(quadrant, useNA = "no")),
      n_significant_a = sum(sig_a),
      n_significant_b = sum(sig_b),
      n_significant_both = sum(both),
      n_concordant = sum(is_sig),
      n_discordant = sum(both & direction %in% "discordant")
    )
  )
}
