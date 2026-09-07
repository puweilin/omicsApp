# Per-feature concordance between two diff_bundles. For each gene/protein
# that appears in both layers, joins the standardized effect, direction,
# and (adjusted) p-value, classifies the (dir_a, dir_b) pair into a
# quadrant, and computes a combined-p (Fisher's method).

# Where each side came from. Not in `joined`, which holds per-feature
# columns; these are per-analysis, and they are what lets a reader tell
# whether the two effects are even the same kind of quantity -- a log2 fold
# change against a correlation is not a comparison.
side_effect_type <- function(bundle) {
  et <- bundle$results$diff_result_df$effect_type
  if (length(et) == 0L) return(NA_character_)
  as.character(et[[1L]])
}

side_analysis_id <- function(bundle) {
  bundle$provenance_binding$analysis_id %||% NA_character_
}

run_integration_concordance <- function(
  project,
  experiments,
  diff_bundles,
  by = "feature_symbol",
  p_preference = c("adjusted", "raw"),
  p_cutoff = 0.05,
  p_adjust_method = "BH",
  combine_p = FALSE
) {
  experiments <- resolve_experiment_pair(project, experiments)
  p_preference <- match.arg(p_preference)
  validate_diff_bundles(diff_bundles, experiments)

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

  build_side <- function(df, side) {
    sub <- df[, c(by, "feature_id", "effect", "direction", p_col),
              drop = FALSE]
    sub <- sub[!is.na(sub[[by]]) & nzchar(sub[[by]]), , drop = FALSE]
    sub <- sub[!duplicated(sub[[by]]), , drop = FALSE]
    names(sub) <- c("key", paste0("feature_id_", side),
                    paste0("effect_", side),
                    paste0("direction_", side),
                    paste0("p_", side))
    sub
  }

  joined <- merge(build_side(res_a, "a"), build_side(res_b, "b"), by = "key")
  if (nrow(joined) == 0L) {
    return(list(
      std = new_integration_result_template(),
      info = list(experiments = experiments, n_features = 0L)
    ))
  }

  quadrant <- classify_concordance_quadrant(joined$direction_a, joined$direction_b)
  direction <- ifelse(
    is.na(quadrant), NA_character_,
    ifelse(quadrant %in% c("up_up", "down_down"), "concordant", "discordant")
  )

  # Effect = difference of (signed) effects. Useful for plotting and a
  # rough magnitude signal.
  effect <- joined$effect_a - joined$effect_b

  pa <- joined$p_a
  pb <- joined$p_b

  # Fisher's method combines p-values that are independent and uniform
  # under the null. Neither holds here as it is being used: `p_preference`
  # defaults to the adjusted p, which is not a p-value in that sense, and
  # two omics layers measured on the same subjects are not independent.
  # Adjusting the combination afterwards does not repair either.
  #
  # So the combination is computed but not reported as significance, and
  # nothing downstream reads it as one, until the statistical basis is
  # settled. `combine_p` is the switch for that decision, not a preference:
  # turning it on is asserting the assumptions hold for the data at hand.
  #
  # What is reported instead is the thing that needs no such assumption:
  # each layer's own p, and whether the feature is significant in both and
  # moving the same way. That is a description of agreement rather than a
  # test of a joint hypothesis, and it is what the quadrant plot already
  # showed.
  combined_p <- rep(NA_real_, nrow(joined))
  combined_adj <- rep(NA_real_, nrow(joined))
  if (isTRUE(combine_p)) {
    ok <- !is.na(pa) & !is.na(pb) & pa > 0 & pb > 0
    if (any(ok)) {
      chisq <- -2 * (log(pa[ok]) + log(pb[ok]))
      combined_p[ok] <- stats::pchisq(chisq, df = 4L, lower.tail = FALSE)
    }
    combined_adj <- stats::p.adjust(combined_p, method = p_adjust_method)
  }

  # Significant in both layers, in the same direction. Named for what it
  # is: `is_significant` alone read as a verdict on the pair, and it was a
  # different judgement from the combined p sitting beside it -- two
  # columns, two criteria, one name to confuse them under.
  concordant_sig <- !is.na(pa) & !is.na(pb) &
                    pa < p_cutoff & pb < p_cutoff &
                    !is.na(direction) & direction == "concordant"
  is_sig <- concordant_sig

  out <- data.frame(
    feature_id = joined$key,
    feature_symbol = joined$key,
    result_type = "concordance",
    experiments = paste(tag_a, "vs", tag_b),
    comparison = paste(
      diff_bundles[[tag_a]]$params$comparison %||% "comparison",
      diff_bundles[[tag_b]]$params$comparison %||% "comparison",
      sep = " | "
    ),
    effect = effect,
    effect_type = "effect_diff",
    # Each layer's own effect and p, kept rather than left to be worked
    # back out of the difference. The plotting helper used to do exactly
    # that -- with one number and two unknowns, so it set the first axis to
    # the difference and the second to zero, and labelled them as the two
    # layers' effects. Every point sat on a line the data never described.
    effect_a = joined$effect_a,
    effect_b = joined$effect_b,
    p_a = pa,
    p_b = pb,
    # Which p these are. `p_preference` decides whether the columns above
    # hold raw or adjusted values, and a reader cannot tell from the name.
    p_type = p_preference,
    effect_type_a = side_effect_type(diff_bundles[[tag_a]]),
    effect_type_b = side_effect_type(diff_bundles[[tag_b]]),
    analysis_id_a = side_analysis_id(diff_bundles[[tag_a]]),
    analysis_id_b = side_analysis_id(diff_bundles[[tag_b]]),
    # The Fisher statistic is part of the combination, so it goes with it.
    # Reporting the chi-square while withholding the p it produces would be
    # withholding the label rather than the claim.
    statistic = if (isTRUE(combine_p)) {
      -2 * (log(pmax(pa, .Machine$double.xmin)) +
            log(pmax(pb, .Machine$double.xmin)))
    } else {
      rep(NA_real_, nrow(joined))
    },
    statistic_type = if (isTRUE(combine_p)) "fisher_chisq" else NA_character_,
    p_value = combined_p,
    adj_p_value = combined_adj,
    direction = direction,
    quadrant = quadrant,
    is_significant = is_sig,
    # The same judgement under a name that says which one it is. Kept
    # alongside rather than instead of, so a reader of an older result and
    # a reader of a new one are looking at the same column.
    is_concordant_significant = concordant_sig,
    is_combined_significant = if (isTRUE(combine_p)) {
      !is.na(combined_adj) & combined_adj < p_cutoff
    } else {
      rep(NA, nrow(joined))
    },
    source_label = paste0("integration_concordance_", tag_a, "_", tag_b),
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL

  list(
    std = out,
    info = list(
      experiments = experiments,
      n_features = nrow(out),
      p_preference = p_preference,
      quadrant_counts = as.list(table(quadrant, useNA = "no"))
    )
  )
}
