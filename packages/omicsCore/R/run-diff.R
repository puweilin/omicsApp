# Public entry-point for the differential-expression layer. Dispatches to the
# backend-specific functions in diff-ttest.R / diff-lm.R / diff-limma.R /
# diff-deseq2.R / diff-edger.R and wraps the result in an analysis_bundle.

#' Differential backends this package implements
#'
#' Every value [run_diff()] accepts for `method`. Which of them are
#' *appropriate* for a given dataset is a narrower question — see
#' [applicable_diff_methods()].
#'
#' @format Character vector.
#' @export
#' @family diff
SUPPORTED_DIFF_METHODS <- c("auto", "deseq2", "edger", "limma", "ttest", "lm")

# `SUPPORTED_DIFF_ANALYSIS_TYPES` lives in constants.R. It used to be
# defined here as well, with the same value; whichever file collates
# last won, which is not a way to hold a constant.

#' Differential methods that are valid for an input
#'
#' `SUPPORTED_DIFF_METHODS` lists every backend that exists;
#' this lists the ones whose assumptions the data actually meets.
#'
#' DESeq2 and edgeR model raw counts as negative binomial. Given
#' continuous intensities they do not refuse — DESeq2 rounds to integers
#' ("converting counts to integer mode") and reports p-values for a
#' model the data never fitted. Conversely, this package's limma backend
#' does not apply voom, so it has no business being handed raw counts.
#' Both mistakes produce a full, plausible result table and no error,
#' which is the failure mode worth engineering against.
#'
#' Written next to `auto_select_diff_method()` on purpose: one decides
#' what to run by default and the other what a caller may choose, and
#' the two disagreeing would be its own bug.
#'
#' @param input An [`omics_input`][omics_input()].
#' @param analysis_type One of [SUPPORTED_DIFF_ANALYSIS_TYPES].
#'
#' @return Character vector of method names, always including `"auto"`.
#' @export
#' @family diff
#' @examples
#' \dontrun{
#'   applicable_diff_methods(rnaseq_counts_input)  # deseq2, edger, ...
#'   applicable_diff_methods(proteomics_input)     # limma, ttest, lm
#' }
applicable_diff_methods <- function(input, analysis_type = "group") {
  # Validated rather than duck-typed: this answer gates which engines a
  # user may run, and quietly returning the continuous set for a NULL or
  # a mistyped variable would hand back a plausible answer to a question
  # that was never asked.
  if (!is_omics_input(input)) {
    stop("`input` must be an `omics_input`.")
  }
  is_counts <- identical(input$omics_type, "rnaseq") &&
    identical(input$assay_type, "raw_count")
  methods <- if (is_counts) c("deseq2", "edger") else c("limma", "ttest", "lm")
  # Only the regression backends carry a continuous predictor.
  if (identical(analysis_type, "continuous")) {
    methods <- intersect(methods, c("limma", "lm"))
  }
  c("auto", methods)
}

# Auto-pick a method given omics_type + analysis_type. If the preferred
# Bioconductor backend is not installed, fall back to base-R (ttest/lm) and
# emit a message so the caller knows.
auto_select_diff_method <- function(input, analysis_type) {
  prefer <- if (input$omics_type == "rnaseq" &&
                identical(input$assay_type, "raw_count")) {
    "deseq2"
  } else if (input$omics_type == "proteomics") {
    "limma"
  } else {
    NA_character_
  }

  ok <- switch(
    as.character(prefer),
    "deseq2" = is_installed("DESeq2"),
    "limma"  = is_installed("limma"),
    FALSE
  )

  if (isTRUE(ok)) return(prefer)

  fallback <- if (analysis_type == "continuous") "lm" else "ttest"
  if (!is.na(prefer)) {
    message(
      "`method = 'auto'`: preferred backend '", prefer,
      "' is not installed; falling back to '", fallback,
      "'. Install with: omicsCore::install_optional(",
      if (prefer == "deseq2") "'rnaseq'" else "'proteomics'", ")."
    )
  }
  fallback
}

dispatch_diff_backend <- function(input, method, analysis_type, args) {
  call_backend <- function(fun, args) do.call(fun, c(list(input = input), args))

  if (method == "limma" && analysis_type == "group") {
    return(call_backend(run_limma_group, args))
  }
  if (method == "limma" && analysis_type == "continuous") {
    return(call_backend(run_limma_continuous, args))
  }
  if (method == "limma" && analysis_type == "anova") {
    return(call_backend(run_limma_anova, args))
  }
  if (method == "deseq2" && analysis_type == "group") {
    return(call_backend(run_deseq2_group, args))
  }
  if (method == "deseq2" && analysis_type == "continuous") {
    return(call_backend(run_deseq2_continuous, args))
  }
  if (method == "edger" && analysis_type == "group") {
    return(call_backend(run_edger_group, args))
  }
  if (method == "ttest" && analysis_type == "group") {
    return(call_backend(run_ttest_group, args))
  }
  if (method == "lm" && analysis_type == "group") {
    return(call_backend(run_lm_group, args))
  }
  if (method == "lm" && analysis_type == "continuous") {
    return(call_backend(run_lm_continuous, args))
  }

  stop(
    "Unsupported combination: omics_type = ", input$omics_type,
    ", method = ", method, ", analysis_type = ", analysis_type
  )
}

#' Run a differential-expression analysis
#'
#' Single entry point for proteomics and RNA-seq differential analysis.
#' Dispatches to the appropriate backend (DESeq2, edgeR, limma, t-test, or
#' linear model) and returns an [`analysis_bundle`][is_analysis_bundle()]
#' wrapping the standardized result, the backend's raw table, and the fitted
#' model object.
#'
#' `method = "auto"` picks limma for proteomics and DESeq2 for raw-count
#' RNA-seq. When the preferred Bioconductor backend is not installed, it
#' silently falls back to t-test / lm (with a message) so analyses still run
#' in restricted environments without `omicsCore::install_optional()` being
#' invoked first.
#'
#' For paired designs supply `paired_col`; for ANOVA-style multi-group tests
#' set `analysis_type = "anova"` (currently limma-backed only). The
#' `continuous` analysis type requires `continuous_col` instead of
#' `group_col` + `control_group` + `case_group`.
#'
#' @param input A validated `omics_input`.
#' @param method Backend name. `"auto"` (default) lets `omicsCore` pick one
#'   based on `omics_type` and installed Suggests.
#' @param analysis_type One of `"group"`, `"continuous"`, or `"anova"`.
#' @param group_col Group column in sample metadata (group/anova).
#' @param control_group Control-group label (group only).
#' @param case_group Case-group label (group only).
#' @param continuous_col Continuous metadata column (continuous only).
#' @param covariates Optional character vector of covariate column names.
#' @param paired_col Optional pairing/block column.
#' @param selected_groups Optional subset of groups to retain (anova only).
#' @param ... Extra arguments forwarded to the backend, e.g. `var_equal` for
#'   t-test or `df = 3` for limma spline.
#'
#' @return An [`analysis_bundle`][is_analysis_bundle()] with
#'   `results$diff_result_df` (standardized schema), `results$diff_raw_df`
#'   (backend-native), and `results$diff_object` (fitted model, may be
#'   `NULL`).
#' @export
#' @family diff
#' @examples
#' \dontrun{
#'   res <- run_diff(input,
#'                   analysis_type = "group",
#'                   group_col = "treatment",
#'                   control_group = "DMSO",
#'                   case_group = "Drug")
#'   head(res$results$diff_result_df)
#' }
run_diff <- function(
  input,
  method = "auto",
  analysis_type = c("group", "continuous", "anova"),
  group_col = NULL,
  control_group = NULL,
  case_group = NULL,
  continuous_col = NULL,
  covariates = NULL,
  paired_col = NULL,
  selected_groups = NULL,
  ...
) {
  validate_omics_input(input)
  analysis_type <- match.arg(analysis_type)
  method <- match.arg(method, choices = SUPPORTED_DIFF_METHODS)

  if (method == "auto") {
    method <- auto_select_diff_method(input, analysis_type)
  }

  backend_args <- switch(
    analysis_type,
    group = list(
      group_col = group_col,
      control_group = control_group,
      case_group = case_group,
      covariates = covariates,
      paired_col = paired_col
    ),
    continuous = list(
      continuous_col = continuous_col,
      covariates = covariates,
      paired_col = paired_col
    ),
    anova = list(
      group_col = group_col,
      covariates = covariates,
      selected_groups = selected_groups,
      paired_col = paired_col
    )
  )

  validate_diff_args(analysis_type, method, backend_args)

  # Drop arguments the chosen backend doesn't accept (e.g. ttest has no
  # `covariates`, lm/ttest have no `paired_col`-via-limma corfit, ...).
  backend_args <- prune_backend_args(method, analysis_type, backend_args)

  extra_args <- list(...)
  backend_result <- dispatch_diff_backend(
    input = input,
    method = method,
    analysis_type = analysis_type,
    args = c(backend_args, extra_args)
  )

  new_analysis_bundle(
    analysis_name = "run_diff",
    input_info = list(
      omics_type = input$omics_type,
      assay_type = input$assay_type,
      n_samples = ncol(input$expr_mat),
      n_features = nrow(input$expr_mat)
    ),
    params = c(
      list(
        method = method,
        analysis_type = analysis_type,
        comparison = backend_result$analysis_info$comparison
      ),
      backend_args,
      extra_args
    ),
    results = list(
      diff_result_df = backend_result$results_std,
      diff_raw_df = backend_result$results_raw,
      diff_object = backend_result$model_object
    )
  )
}

#' Continuous-variable differential analysis
#'
#' Convenience wrapper around [run_diff()] that pins
#' `analysis_type = "continuous"`.
#'
#' @inheritParams run_diff
#'
#' @return An [`analysis_bundle`][is_analysis_bundle()].
#' @export
#' @family diff
run_diff_continuous <- function(
  input,
  method = "auto",
  continuous_col,
  covariates = NULL,
  paired_col = NULL,
  ...
) {
  run_diff(
    input = input,
    method = method,
    analysis_type = "continuous",
    continuous_col = continuous_col,
    covariates = covariates,
    paired_col = paired_col,
    ...
  )
}

# ---- internal helpers --------------------------------------------------

validate_diff_args <- function(analysis_type, method, args) {
  if (analysis_type %in% c("group", "anova")) {
    if (is.null(args$group_col)) {
      stop("`group_col` is required for analysis_type = '", analysis_type, "'.")
    }
  }
  if (analysis_type == "group") {
    if (is.null(args$control_group) || is.null(args$case_group)) {
      stop("`control_group` and `case_group` are required for analysis_type = 'group'.")
    }
  }
  if (analysis_type == "continuous") {
    if (is.null(args$continuous_col)) {
      stop("`continuous_col` is required for analysis_type = 'continuous'.")
    }
  }
  if (analysis_type == "anova" && method != "limma") {
    stop("ANOVA analysis_type currently only supports method = 'limma'.")
  }
  if (method == "edger" && analysis_type != "group") {
    stop("edgeR backend currently only supports analysis_type = 'group'.")
  }
  invisible(TRUE)
}

# Strip args the backend's signature doesn't accept. This keeps the public
# entry point uniform (callers can always pass `covariates`/`paired_col`)
# while keeping each backend's `...` honest.
prune_backend_args <- function(method, analysis_type, args) {
  drop <- character(0)
  if (method == "ttest") {
    drop <- c(drop, "covariates")
  }
  if (method == "lm") {
    drop <- c(drop, "paired_col")
  }
  args[setdiff(names(args), drop)]
}
