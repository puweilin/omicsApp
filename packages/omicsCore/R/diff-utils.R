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

#' Refuse limma input that is not on the scale limma reports
#'
#' limma fits a normal linear model to the values as given and labels the
#' group effect `log2FC`. Hand it linear intensities and it does not complain:
#' the difference of two means comes back as a fold change, six orders of
#' magnitude wide, and the volcano plot, the effect threshold and everything
#' downstream reads it as one.
#'
#' Scoped to limma on purpose. DESeq2 and edgeR model counts on the count
#' scale and their `log2FC` is earned; `lm` transforms `raw_count` itself and
#' reports `beta`; `ttest` reports `mean_diff`, which is what it computes
#' whatever the scale. None of those label the effect as something the input
#' cannot support, so none of them are narrowed here.
#'
#' Where magnitude can tell the two scales apart -- the intensity family, per
#' [SCALE_CHECKED_ASSAY_TYPES] -- the values get the last word over the label.
#' A project saved before the vocabulary existed carries `"intensity"`, which
#' the import view stamped on every proteomics upload whether or not it had
#' been transformed, so that label asserts nothing and refusing on it alone
#' would make old log2 projects unanalysable. TPM and FPKM are excluded from
#' that rescue for the reason counts are: a gene at 5 TPM is ordinary, so a
#' small maximum is not evidence of a log scale.
#'
#' @param input A validated `omics_input`.
#' @param method The backend actually chosen, i.e. after `"auto"` resolves.
#'
#' @return Invisibly `TRUE`; stops otherwise.
#' @keywords internal
#' @noRd
assert_limma_scale <- function(input, method) {
  if (!identical(method, "limma")) return(invisible(TRUE))
  at <- input$assay_type %||% NA_character_
  if (is.na(at) || at %in% LOG_SCALE_ASSAY_TYPES) return(invisible(TRUE))

  # Not canonical_assay_type(): that warns, and validate_omics_input() has
  # already said this once by the time a backend is chosen.
  canonical <- if (at %in% names(DEPRECATED_ASSAY_TYPE_ALIASES)) {
    unname(DEPRECATED_ASSAY_TYPE_ALIASES[[at]])
  } else {
    at
  }
  if (canonical %in% LOG_SCALE_ASSAY_TYPES) return(invisible(TRUE))

  if (canonical %in% SCALE_CHECKED_ASSAY_TYPES) {
    max_value <- suppressWarnings(max(input$expr_mat, na.rm = TRUE))
    if (!is.finite(max_value) || max_value <= MAX_PLAUSIBLE_LOG_SCALE_VALUE) {
      return(invisible(TRUE))
    }
  }

  fix <- if (identical(input$omics_type, "proteomics")) {
    "Normalize the layer first with normalize_omics(), which puts proteomics intensities on a log2 scale."
  } else {
    # normalize_omics() is proteomics-only, so pointing RNA users at it would
    # be advice they cannot follow.
    "Use a log-scale layer (vst or logcpm) for limma, or a backend that models this scale: deseq2 / edger for raw counts, lm for a linear predictor."
  }
  stop(sprintf(
    "limma reports a log2 fold change, but `assay_type = '%s'` is a linear scale, so the values it would report are differences of untransformed measurements. %s",
    at, fix), call. = FALSE)
}
