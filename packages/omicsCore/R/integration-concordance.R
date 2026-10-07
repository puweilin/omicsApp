# Per-feature concordance between two diff_bundles. For each protein-gene
# pair the two layers share, joins the standardized effect, direction,
# and (adjusted) p-value, classifies the (dir_a, dir_b) pair into a
# quadrant, and computes a combined-p (Fisher's method).
#
# Pairs come from link_features(): by symbol, or by a feature link. Where
# a gene is measured by several proteins each protein-gene pair is a row
# of its own, with its own quadrant and combined p, and the combined p is
# corrected across pairs. Keeping only one protein per gene (as before)
# threw away isoforms that can move differently; reporting each lets the
# reader see that they do.

run_integration_concordance <- function(
  project,
  experiments,
  diff_bundles,
  by = "feature_symbol",
  p_preference = c("adjusted", "raw"),
  p_cutoff = 0.05,
  p_adjust_method = "BH",
  effect_cutoff = 0,
  link = NULL
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

  pairs <- link_features(res_a$feature_id, res_a[[by]],
                         res_b$feature_id, res_b[[by]], link = link)
  side_cols <- function(df, idx) {
    df <- df[idx, , drop = FALSE]
    data.frame(
      feature_id = as.character(df$feature_id),
      effect = as.numeric(df$effect),
      direction = direction_sign(df$direction, df$effect),
      p = as.numeric(df$p_value),
      padj = as.numeric(df$adj_p_value),
      stringsAsFactors = FALSE
    )
  }
  sa <- side_cols(res_a, pairs$i_a)
  sb <- side_cols(res_b, pairs$i_b)
  names(sa) <- paste0(names(sa), "_a")
  names(sb) <- paste0(names(sb), "_b")
  joined <- cbind(data.frame(id = pairs$feature_id, symbol = pairs$label,
                             stringsAsFactors = FALSE), sa, sb)
  if (nrow(joined) == 0L) {
    return(list(
      std = new_integration_result_template(),
      info = list(experiments = experiments, n_features = 0L,
                  p_preference = p_preference,
                  feature_pairing = attr(pairs, "info"))
    ))
  }

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
    feature_id = joined$id,
    feature_symbol = joined$symbol,
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
      # Among the features significant in both layers -- the ones whose
      # agreement means anything. `quadrant_counts_all` keeps the count
      # over every feature with a sign in both.
      quadrant_counts = as.list(table(factor(quadrant[both],
        levels = c("up_up", "down_down", "up_down", "down_up")))),
      quadrant_counts_all = as.list(table(quadrant, useNA = "no")),
      n_significant_a = sum(sig_a),
      n_significant_b = sum(sig_b),
      n_significant_both = sum(both),
      n_concordant = sum(is_sig),
      n_discordant = sum(both & direction %in% "discordant"),
      feature_pairing = attr(pairs, "info")
    )
  )
}
