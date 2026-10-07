#' Filter standardized differential results
#'
#' Subsets a standardized diff result table to features that pass a
#' significance cutoff (either adjusted or raw p-value), plus optional
#' absolute-effect and model-fit thresholds. Returns the filtered table with
#' `is_significant = TRUE` set on every retained row.
#'
#' @param result_df A standardized diff result `data.frame`, as produced by
#'   [run_diff()] or any of the backend functions.
#' @param p_cutoff P-value threshold; defaults to 0.05.
#' @param p_preference Whether to threshold on the adjusted or raw p-value.
#' @param effect_cutoff Optional absolute effect cutoff (e.g., 1 for
#'   |log2FC| >= 1).
#' @param model_fit_cutoff Optional model-fit cutoff (e.g., adjusted R^2
#'   threshold for continuous lm/limma fits).
#'
#' @return Filtered standardized diff result `data.frame`.
#' @export
#' @family diff
filter_diff_results <- function(
  result_df,
  p_cutoff = 0.05,
  p_preference = c("adjusted", "raw"),
  effect_cutoff = NULL,
  model_fit_cutoff = NULL
) {
  check_diff_result_schema(result_df)
  p_preference <- match.arg(p_preference)
  assert_number(p_cutoff, "p_cutoff", lower = 0, upper = 1)
  assert_number(effect_cutoff, "effect_cutoff", lower = 0, allow_null = TRUE)
  assert_number(model_fit_cutoff, "model_fit_cutoff", allow_null = TRUE)

  p_col <- resolve_p_col(result_df, p_preference = p_preference)
  out <- result_df[!is.na(result_df[[p_col]]) & result_df[[p_col]] < p_cutoff, , drop = FALSE]

  if (!is.null(effect_cutoff)) {
    out <- out[!is.na(out$effect) & abs(out$effect) >= effect_cutoff, , drop = FALSE]
  }
  if (!is.null(model_fit_cutoff)) {
    out <- out[!is.na(out$model_fit) & out$model_fit >= model_fit_cutoff, , drop = FALSE]
  }

  out$is_significant <- rep(TRUE, nrow(out))
  rownames(out) <- NULL
  out
}

#' Build a named ranked feature vector for preranked GSEA
#'
#' Produces a named numeric vector sorted in decreasing order, suitable as
#' input to preranked GSEA (e.g., `fgsea::fgsea(stats = ...)`). Duplicate
#' feature labels are dropped, keeping the first occurrence after the sort.
#'
#' @param result_df A standardized diff result `data.frame`.
#' @param feature_col Column used to name the vector (defaults to
#'   `"feature_symbol"`).
#' @param rank_col Numeric column used to rank features (defaults to
#'   `"effect"`).
#'
#' @return Named numeric vector sorted decreasing.
#' @export
#' @family diff
make_ranked_features <- function(
  result_df,
  feature_col = "feature_symbol",
  rank_col = "effect"
) {
  check_diff_result_schema(result_df)
  assert_string(feature_col, "feature_col")
  assert_string(rank_col, "rank_col")

  if (!feature_col %in% colnames(result_df)) {
    stop("Feature column not found: ", feature_col)
  }
  if (!rank_col %in% colnames(result_df)) {
    stop("Rank column not found: ", rank_col)
  }

  ranked_df <- result_df[, c(feature_col, rank_col), drop = FALSE]
  ranked_df <- ranked_df[stats::complete.cases(ranked_df), , drop = FALSE]
  ranked_df <- ranked_df[order(ranked_df[[rank_col]], decreasing = TRUE), , drop = FALSE]
  ranked_df <- ranked_df[!duplicated(ranked_df[[feature_col]]), , drop = FALSE]

  out <- ranked_df[[rank_col]]
  names(out) <- ranked_df[[feature_col]]
  out
}

# ---- internal helpers --------------------------------------------------

resolve_p_col <- function(result_df, p_preference = c("adjusted", "raw")) {
  p_preference <- match.arg(p_preference)
  if (p_preference == "adjusted") "adj_p_value" else "p_value"
}

# Coerce a metadata column to numeric. Used by continuous-mode backends to
# accept character-formatted age / dose / etc. columns without a hard
# dependency on readr::parse_number.
coerce_continuous_col <- function(x, col_name) {
  if (is.numeric(x)) return(x)
  parsed <- suppressWarnings(as.numeric(as.character(x)))
  if (all(is.na(parsed))) {
    stop("`", col_name, "` must be numeric or coercible to numeric.")
  }
  parsed
}

# ---- shared backend set-up ----------------------------------------------
#
# limma, DESeq2 and edgeR choose a group comparison's samples, assemble
# its model terms and stack its contrasts in the same way. Each used to
# carry its own copy, so a fix to one (the pairing check, a covariate
# that is missing) had to be made three times.

# The samples a group comparison fits: every group any contrast names,
# in `control_group`-first order, with the pairing checked. Returns the
# contrast specs, the group levels (the first is the reference) and the
# metadata of the samples kept, its group column a factor.
contrast_target_meta <- function(meta_df, group_col, control_group, case_group,
                                 contrasts = NULL, paired_col = NULL) {
  if (!group_col %in% colnames(meta_df)) {
    stop("`group_col` not found in `meta_df`: ", group_col)
  }
  check_paired_col(meta_df, paired_col, object_name = "meta_df")

  specs <- contrasts %||% case_control_contrasts(control_group, case_group)
  group_levels <- contrast_levels(specs, order = control_group)
  meta_df[[group_col]] <- as.character(meta_df[[group_col]])
  target_meta <- meta_df[!is.na(meta_df[[group_col]]) &
                           meta_df[[group_col]] %in% group_levels, , drop = FALSE]
  target_meta[[group_col]] <- factor(target_meta[[group_col]], levels = group_levels)
  validate_two_group_pairing(
    target_meta,
    group_col = group_col,
    paired_col = paired_col,
    control_group = if (is.null(contrasts)) control_group,
    case_group = if (is.null(contrasts)) case_group else group_levels,
    object_name = "target_meta"
  )
  list(specs = specs, group_levels = group_levels, target_meta = target_meta)
}

# The terms of a count model, in the order DESeq2 and edgeR fit them: the
# pairing block first (as a factor), then what is tested, then the
# covariates. `paired_na` is the start of the message for a block with a
# missing value. Returns the metadata (block factored) and the terms.
count_model_terms <- function(meta, primary, paired_col = NULL, covariates = NULL,
                              paired_na = "`paired_col` contains missing values after group filtering: ") {
  terms <- primary
  if (!is.null(paired_col)) {
    if (anyNA(meta[[paired_col]])) stop(paired_na, paired_col)
    meta[[paired_col]] <- factor(meta[[paired_col]])
    terms <- c(paired_col, terms)
  }
  missing_cov <- setdiff(covariates, colnames(meta))
  if (length(missing_cov) > 0L) {
    stop("Missing covariates: ", paste(missing_cov, collapse = ", "))
  }
  list(meta = meta, terms = c(terms, covariates))
}

# `~ term + term`, each backticked: a column called "Treatment Group" was
# an "unexpected symbol" in the formula. No terms is the intercept alone.
backtick_formula <- function(terms) {
  stats::as.formula(paste("~", if (length(terms))
    paste0("`", terms, "`", collapse = " + ") else "1"))
}

# Every contrast read off one fit and stacked into one table under its
# `comparison` label. `raw_for(i)` returns the engine's table for the
# i-th contrast and `standardize(raw_df, comparison)` its standardized
# form. A single contrast keeps the engine's table without a
# `comparison` column, as it always had.
stack_contrast_results <- function(specs, raw_for, standardize) {
  comparisons <- vapply(specs, `[[`, character(1), "label")
  per <- lapply(seq_along(specs), function(i) {
    raw_df <- raw_for(i)
    std <- standardize(raw_df, comparisons[[i]])
    raw_df$comparison <- comparisons[[i]]
    list(raw = raw_df, std = std)
  })
  raw_df <- do.call(rbind, lapply(per, `[[`, "raw"))
  if (length(specs) == 1L) raw_df$comparison <- NULL
  results_std <- do.call(rbind, lapply(per, `[[`, "std"))
  rownames(raw_df) <- NULL
  rownames(results_std) <- NULL
  list(raw = raw_df, std = results_std, comparisons = comparisons)
}
