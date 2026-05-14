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
