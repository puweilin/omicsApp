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
  p_adjust_method = "BH"
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

  # Combined p via Fisher's method, only when both inputs are present.
  pa <- joined$p_a
  pb <- joined$p_b
  combined_p <- rep(NA_real_, nrow(joined))
  ok <- !is.na(pa) & !is.na(pb) & pa > 0 & pb > 0
  if (any(ok)) {
    chisq <- -2 * (log(pa[ok]) + log(pb[ok]))
    combined_p[ok] <- stats::pchisq(chisq, df = 4L, lower.tail = FALSE)
  }
  combined_adj <- stats::p.adjust(combined_p, method = p_adjust_method)

  is_sig <- !is.na(pa) & !is.na(pb) &
            pa < p_cutoff & pb < p_cutoff &
            !is.na(direction) & direction == "concordant"

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
    statistic = -2 * (log(pmax(pa, .Machine$double.xmin)) +
                      log(pmax(pb, .Machine$double.xmin))),
    statistic_type = "fisher_chisq",
    p_value = combined_p,
    adj_p_value = combined_adj,
    direction = direction,
    quadrant = quadrant,
    is_significant = is_sig,
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
