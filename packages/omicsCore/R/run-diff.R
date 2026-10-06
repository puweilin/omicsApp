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
  } else if (input$omics_type == "proteomics" ||
             isTRUE(input$assay_type %in% LOG_SCALE_ASSAY_TYPES)) {
    # Log-scale RNA-seq (vst, logCPM) is what limma is for. It used to
    # fall through to the t-test, which has no covariates: a user who
    # asked for an age-adjusted comparison got an unadjusted one.
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
  if (method == "edger" && analysis_type == "anova") {
    return(call_backend(run_edger_anova, args))
  }
  if (method == "deseq2" && analysis_type == "anova") {
    return(call_backend(run_deseq2_anova, args))
  }
  if (method == "deseq2" && analysis_type == "continuous") {
    return(call_backend(run_deseq2_continuous, args))
  }
  if (method == "edger" && analysis_type == "group") {
    return(call_backend(run_edger_group, args))
  }
  if (method == "ttest" && analysis_type == "group") {
    return(loop_case_groups(run_ttest_group, input, args))
  }
  if (method == "lm" && analysis_type == "group") {
    return(loop_case_groups(run_lm_group, input, args))
  }
  if (method == "lm" && analysis_type == "continuous") {
    return(call_backend(run_lm_continuous, args))
  }

  stop(
    "Unsupported combination: omics_type = ", input$omics_type,
    ", method = ", method, ", analysis_type = ", analysis_type
  )
}

# The t-test and per-feature lm have no shared model to fit several
# groups into, so several case groups are several two-group runs, stacked
# into one table under their `comparison` labels.
loop_case_groups <- function(fun, input, args) {
  specs <- args$contrasts
  args$contrasts <- NULL
  paired_t <- !is.null(args$paired_col) && identical(fun, run_ttest_group)
  if (is.null(specs)) {
    if (length(args$case_group) <= 1L) {
      if (!paired_t) return(do.call(fun, c(list(input = input), args)))
      # One comparison is paired the same way as several: the subjects
      # sampled in both groups, with a note naming any left out. It used
      # to refuse here, so the same Control vs TreatB comparison ran
      # inside a multi-group analysis and failed on its own.
      args$incomplete_pairs <- "drop"
      res <- do.call(fun, c(list(input = input), args))
      note <- paired_comparisons_note(
        list(list(case = args$case_group, control = args$control_group)),
        res$analysis_info$pairs_used, list(res$analysis_info$pairs_left_out),
        args$paired_col)
      res$analysis_info$warnings <- c(res$analysis_info$warnings, note)
      return(res)
    }
    specs <- case_control_contrasts(args$control_group, args$case_group)
  }
  not_pair <- vapply(specs, function(s) is.null(s$case), logical(1))
  if (any(not_pair)) {
    stop("The t-test and lm backends compare two groups at a time; ",
         "the contrast \"", specs[not_pair][[1L]]$spec, "\" needs limma, ",
         "DESeq2 or edgeR.", call. = FALSE)
  }
  # A paired t-test is paired per comparison. With several treatments
  # against one control, a subject who missed one treatment is the
  # ordinary unbalanced design, not an error: each comparison uses the
  # subjects sampled in both of its groups. Demanding every subject in
  # every comparison refused the whole run over one missing sample.
  per_pair <- paired_t
  if (per_pair) args$incomplete_pairs <- "drop"
  runs <- lapply(specs, function(s) {
    a <- args
    a$control_group <- s$control
    a$case_group <- s$case
    do.call(fun, c(list(input = input), a))
  })
  raw <- lapply(seq_along(runs), function(i) {
    r <- runs[[i]]$results_raw
    if (is.data.frame(r)) r$comparison <- runs[[i]]$analysis_info$comparison
    r
  })
  info <- runs[[1L]]$analysis_info
  info$comparison <- vapply(runs, function(r) r$analysis_info$comparison,
                            character(1))
  if (per_pair) {
    info$pairs_used <- stats::setNames(
      vapply(runs, function(r) as.integer(r$analysis_info$pairs_used), integer(1)),
      info$comparison)
    info$pairs_left_out <- stats::setNames(
      lapply(runs, function(r) r$analysis_info$pairs_left_out), info$comparison)
    note <- paired_comparisons_note(specs, info$pairs_used, info$pairs_left_out,
                                    args$paired_col)
    info$warnings <- c(info$warnings, note)
  }
  std <- do.call(rbind, lapply(runs, `[[`, "results_std"))
  rownames(std) <- NULL
  raw <- if (all(vapply(raw, is.data.frame, logical(1)))) do.call(rbind, raw) else raw
  list(
    results_raw = raw,
    results_std = std,
    model_object = lapply(runs, `[[`, "model_object"),
    analysis_info = info
  )
}

# How many pairs each comparison of a paired t-test used, in words. Said
# only when some comparison left pairs out: a balanced design used every
# pair everywhere, which needs no note.
paired_comparisons_note <- function(specs, used, left_out, paired_col) {
  if (!any(lengths(left_out) > 0L)) return(NULL)
  parts <- vapply(seq_along(specs), function(i) {
    out <- sprintf("%s vs %s used %d pairs", specs[[i]]$case, specs[[i]]$control,
                   used[[i]])
    lo <- left_out[[i]]
    if (length(lo)) {
      shown <- paste(utils::head(lo, 5L), collapse = ", ")
      if (length(lo) > 5L) shown <- paste0(shown, ", ...")
      out <- sprintf("%s (%d left out, without one sample in each group: %s)",
                     out, length(lo), shown)
    }
    out
  }, character(1))
  sprintf("Paired t-test by `%s`, pairing each comparison on its own: %s.",
          paired_col, paste(parts, collapse = "; "))
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
#' set `analysis_type = "anova"` (limma for continuous data, edgeR's QL F-test or DESeq2's LRT for raw counts). The
#' `continuous` analysis type requires `continuous_col` instead of
#' `group_col` + `control_group` + `case_group`.
#'
#' @param input A validated `omics_input`.
#' @param method Backend name. `"auto"` (default) lets `omicsCore` pick one
#'   based on `omics_type` and installed Suggests.
#' @param analysis_type One of `"group"`, `"continuous"`, or `"anova"`.
#' @param group_col Group column in sample metadata (group/anova).
#' @param control_group Control-group label (group only).
#' @param case_group Case-group label (group only). Several labels run
#'   every one of them against `control_group` in one call: limma, DESeq2
#'   and edgeR fit all the groups in a single model (so the contrasts share
#'   one variance / dispersion estimate), the t-test and lm backends run
#'   each pair in turn. The contrasts are stacked in `diff_result_df` under
#'   their `comparison` labels; use [select_comparison()] to take one out.
#' @param contrasts Optional comparisons between groups, instead of
#'   `control_group` / `case_group`: a character vector of expressions over
#'   the levels of `group_col` (`"B - A"`, `"(B + C)/2 - A"`; a level that
#'   is not a plain word goes in backticks, `` "`Drug A` - Control" ``), or
#'   `"pairwise"` for every pair of groups (see [pairwise_contrasts()]).
#'   All groups named are fitted in one model. Weighted contrasts need
#'   limma, DESeq2 or edgeR.
#' @param continuous_col Continuous metadata column (continuous only).
#' @param covariates Optional character vector of covariate column names.
#' @param paired_col Optional pairing/block column.
#' @param selected_groups Optional subset of groups to retain (anova only).
#' @param ... Extra arguments forwarded to the backend, e.g. `var_equal` for
#'   t-test, or `model = "spline", df = 3` for a limma spline fit.
#'
#' @return An [`analysis_bundle`][is_analysis_bundle()] with
#'   `results$diff_result_df` (standardized schema), `results$diff_raw_df`
#'   (backend-native), and `results$diff_object` (the fitted model for
#'   limma, t-test and lm; `NULL` for DESeq2 and edgeR unless
#'   `options(omicsCore.keep_count_models = TRUE)` -- see
#'   `keep_model_object()`).
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
  contrasts = NULL,
  ...
) {
  validate_omics_input(input)
  analysis_type <- match.arg(analysis_type)
  method <- match.arg(method, choices = SUPPORTED_DIFF_METHODS)

  # An empty selection means no covariate, which is what NULL says.
  if (length(covariates) == 0L) covariates <- NULL
  if (length(selected_groups) == 0L) selected_groups <- NULL
  assert_string(group_col, "group_col", allow_null = TRUE)
  assert_label(control_group, "control_group", allow_null = TRUE)
  assert_labels(case_group, "case_group", allow_null = TRUE)
  assert_string(continuous_col, "continuous_col", allow_null = TRUE)
  assert_names(covariates, "covariates", allow_null = TRUE)
  assert_string(paired_col, "paired_col", allow_null = TRUE)
  assert_names(selected_groups, "selected_groups", allow_null = TRUE)
  assert_character(contrasts, "contrasts", allow_null = TRUE)
  if (length(contrasts) == 0L) contrasts <- NULL
  if (!is.null(contrasts) && analysis_type != "group") {
    stop("`contrasts` applies to analysis_type = 'group' only.", call. = FALSE)
  }

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
      paired_col = paired_col,
      contrasts = contrasts
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
  # Contrasts are parsed against the levels actually in the column, and
  # recorded in params as written (so a script can repeat the call); the
  # parsed weights travel to the backend alongside.
  specs <- NULL
  if (!is.null(contrasts)) {
    if (!group_col %in% colnames(input$meta_df)) {
      stop("`group_col` not found in `meta_df`: ", group_col, call. = FALSE)
    }
    present <- sort(unique(as.character(stats::na.omit(input$meta_df[[group_col]]))))
    if (identical(contrasts, "pairwise")) {
      # The control first, so every comparison with it reads "X vs control".
      ctrl <- intersect(as.character(control_group), present)
      contrasts <- pairwise_contrasts(c(ctrl, setdiff(present, ctrl)))
    }
    specs <- parse_diff_contrasts(contrasts, present)
    backend_args$contrasts <- vapply(specs, `[[`, character(1), "spec")
  }
  validate_diff_design(input, analysis_type, backend_args, specs = specs,
                       method = method)
  pre <- preflight_diff_matrix(input, method, analysis_type, backend_args,
                               specs = specs)
  input <- pre$input
  original_assay <- input$assay_type
  scale <- prepare_diff_scale(input, method)
  input <- scale$input
  if (!is.null(scale$note)) pre$warnings <- c(pre$warnings, scale$note)

  # Drop arguments the chosen backend doesn't accept (e.g. ttest has no
  # `covariates`, lm/ttest have no `paired_col`-via-limma corfit, ...).
  # Said out loud: an adjustment that was asked for and not made is a
  # different analysis from the one the label on the result describes.
  pruned <- setdiff(names(Filter(Negate(is.null), backend_args)),
                    names(prune_backend_args(method, analysis_type, backend_args)))
  backend_args <- prune_backend_args(method, analysis_type, backend_args)
  if (length(pruned)) {
    note <- sprintf("method = '%s' does not support %s; it was ignored.",
                    method, paste(sprintf("`%s`", pruned), collapse = ", "))
    warning(note, call. = FALSE)
    pre$warnings <- c(pre$warnings, note)
  }

  extra_args <- list(...)
  dispatch_args <- backend_args
  dispatch_args$contrasts <- specs
  backend_result <- dispatch_diff_backend(
    input = input,
    method = method,
    analysis_type = analysis_type,
    args = c(dispatch_args, extra_args)
  )

  new_analysis_bundle(
    analysis_name = "run_diff",
    input_info = list(
      omics_type = input$omics_type,
      assay_type = original_assay,
      analysed_scale = input$assay_type,
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
      extra_args,
      # Which groups each comparison sets against which, so one can be
      # taken out and repeated elsewhere (select_comparison(), and the
      # Integration view repeating it on another layer).
      if (!is.null(specs)) list(contrast_table = data.frame(
        comparison = vapply(specs, `[[`, character(1), "label"),
        spec = vapply(specs, `[[`, character(1), "spec"),
        case = vapply(specs, function(s) s$case %||% NA_character_, character(1)),
        control = vapply(specs, function(s) s$control %||% NA_character_, character(1)),
        stringsAsFactors = FALSE))
    ),
    results = list(
      diff_result_df = backend_result$results_std,
      diff_raw_df = backend_result$results_raw,
      diff_object = keep_model_object(backend_result$model_object)
    ),
    warnings = c(pre$warnings, backend_result$analysis_info$warnings)
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

# The design, checked before any engine sees it.
#
# Every engine answers a degenerate design in its own way, and most of
# the answers are tables. limma and lm drop a covariate that is
# confounded with the group and return the unadjusted result labelled
# as adjusted; lm returns a table of NA for a constant covariate; both
# return all-NA for a continuous variable that does not vary; a case
# level that is not in the column reaches limma as "trying to take
# contrast of non-estimable coefficient". DESeq2 and edgeR refuse
# outright, which is the right answer, so it is now the answer
# everywhere, in words that name the column.
validate_diff_design <- function(input, analysis_type, args, specs = NULL,
                                 method = NULL) {
  meta <- input$meta_df

  if (analysis_type == "anova") {
    # The global test used to skip every check below, so a covariate with
    # a missing value reached the engines as "row dimension of design
    # doesn't match" (limma) or "nrow(design) disagrees with ncol(y)"
    # (edgeR), and a misspelt selected group was dropped without a word.
    group_col <- args$group_col
    if (!group_col %in% colnames(meta)) {
      stop("`group_col` not found in `meta_df`: ", group_col, call. = FALSE)
    }
    present <- unique(as.character(stats::na.omit(meta[[group_col]])))
    sel <- args$selected_groups
    unknown <- setdiff(as.character(sel), present)
    if (length(unknown)) {
      stop(sprintf("`selected_groups` names groups that are not in `%s`: %s. Levels present: %s.",
                   group_col, paste(sprintf("'%s'", unknown), collapse = ", "),
                   paste(sprintf("'%s'", present), collapse = ", ")), call. = FALSE)
    }
    keep <- !is.na(meta[[group_col]]) &
      (is.null(sel) | meta[[group_col]] %in% sel)
    sub <- droplevels(meta[keep, , drop = FALSE])
    primary <- factor(as.character(sub[[group_col]]))
    if (nlevels(primary) < 2L) {
      stop("ANOVA needs at least two groups in `", group_col, "`.", call. = FALSE)
    }
  } else if (analysis_type == "group" && !is.null(specs)) {
    used <- contrast_levels(specs)
    keep <- !is.na(meta[[args$group_col]]) & meta[[args$group_col]] %in% used
    sub <- droplevels(meta[keep, , drop = FALSE])
    primary <- factor(sub[[args$group_col]], levels = used)
  } else if (analysis_type == "group") {
    group_col <- args$group_col
    if (!group_col %in% colnames(meta)) {
      stop("`group_col` not found in `meta_df`: ", group_col, call. = FALSE)
    }
    levels_present <- unique(as.character(meta[[group_col]]))
    levels_present <- levels_present[!is.na(levels_present)]
    for (nm in c("control_group", "case_group")) {
      for (lv in args[[nm]]) {
        if (!lv %in% levels_present) {
          stop(sprintf(
            "`%s` '%s' is not a level of `%s`. Levels present: %s.",
            nm, lv, group_col,
            paste(sprintf("'%s'", levels_present), collapse = ", ")
          ), call. = FALSE)
        }
      }
    }
    if (as.character(args$control_group) %in% as.character(args$case_group)) {
      stop("`control_group` and `case_group` must be distinct.", call. = FALSE)
    }
    keep <- !is.na(meta[[group_col]]) &
      meta[[group_col]] %in% c(args$control_group, args$case_group)
    # droplevels(): a factor covariate level that only occurs in a group
    # outside the comparison is an all-zero column, which the rank check
    # below read as "confounded" on a perfectly balanced design.
    sub <- droplevels(meta[keep, , drop = FALSE])
    primary <- factor(sub[[group_col]],
                      levels = as.character(c(args$control_group, args$case_group)))
  } else {
    cont <- args$continuous_col
    if (!cont %in% colnames(meta)) {
      stop("`continuous_col` not found in `meta_df`: ", cont, call. = FALSE)
    }
    sub <- meta
    values <- coerce_continuous_col(sub[[cont]], cont)
    if (anyNA(values)) {
      stop(sprintf(
        "`%s` has missing values in %d sample(s): %s.",
        cont, sum(is.na(values)),
        paste(utils::head(rownames(sub)[is.na(values)], 5L), collapse = ", ")
      ), call. = FALSE)
    }
    if (length(unique(values)) < 2L) {
      stop(sprintf("`%s` has no variation across samples; there is nothing to regress on.",
                   cont), call. = FALSE)
    }
    primary <- values
  }

  covariates <- args$covariates
  if (is.null(covariates) || length(covariates) == 0L) return(invisible(TRUE))
  # A pairing block absorbs anything constant within a block (a subject's
  # sex, say). limma fits the block as a correlation and is unaffected;
  # DESeq2 and edgeR fit it as fixed effects, and a covariate it absorbs
  # stopped them with "the model matrix is not full rank" /
  # "coefficients not estimable".
  block <- args$paired_col
  block_fixed <- !is.null(block) && isTRUE(method %in% c("deseq2", "edger", "lm")) &&
    block %in% colnames(sub)
  missing_cov <- setdiff(covariates, colnames(sub))
  if (length(missing_cov) > 0L) {
    stop("Missing covariates: ", paste(missing_cov, collapse = ", "), call. = FALSE)
  }
  design_df <- data.frame(.primary = primary, sub[, covariates, drop = FALSE],
                          check.names = FALSE)
  base <- stats::model.matrix(~ .primary, data = design_df)
  for (cov in covariates) {
    x <- sub[[cov]]
    if (anyNA(x)) {
      stop(sprintf(
        "Covariate `%s` has missing values in %d of the samples being compared: %s.",
        cov, sum(is.na(x)),
        paste(utils::head(rownames(sub)[is.na(x)], 5L), collapse = ", ")
      ), call. = FALSE)
    }
    if (length(unique(x)) < 2L) {
      stop(sprintf(
        "Covariate `%s` is constant across the samples being compared, so there is nothing to adjust for.",
        cov), call. = FALSE)
    }
    # Confounded with what is being tested: adding the covariate to the
    # design adds no rank, so the engine would silently drop one of them.
    with_cov <- stats::model.matrix(
      stats::as.formula(paste0("~ .primary + `", cov, "`")), data = design_df)
    if (qr(with_cov)$rank < ncol(with_cov)) {
      what <- if (analysis_type %in% c("group", "anova")) args$group_col else args$continuous_col
      stop(sprintf(
        "Covariate `%s` is confounded with `%s`: the two cannot be separated, so the effect of `%s` is not estimable with it in the model.",
        cov, what, what), call. = FALSE)
    }
    if (block_fixed) {
      bdf <- data.frame(.primary = primary, .block = factor(sub[[block]]),
                        sub[, cov, drop = FALSE], check.names = FALSE)
      with_block <- stats::model.matrix(
        stats::as.formula(paste0("~ .block + .primary + `", cov, "`")), data = bdf)
      if (qr(with_block)$rank < ncol(with_block)) {
        stop(sprintf(
          "Covariate `%s` does not vary within the blocks of `%s` (e.g. a subject's sex in a paired design), so the pairing already accounts for it. Remove the covariate, or use method = 'limma', which models the pairing as a correlation.",
          cov, block), call. = FALSE)
      }
    }
  }
  invisible(TRUE)
}

# The scale the continuous engines need, made so rather than assumed.
#
# limma, the t-test and lm model values as they come. Handed linear
# intensities, TPM or raw counts they used to report the raw difference
# of means as "log2FC" -- a median effect of 1,102,463 for a true 2-fold
# change -- and on counts the t-test and lm compared log2(count + 1)
# without any library-size step, so a deeper-sequenced group came out
# "up" in 347 of 500 null genes. Linear assays are now log2-transformed
# and counts become log2-CPM (TMM-scaled when edgeR is available) before
# such an engine sees them, and the bundle says so.
prepare_diff_scale <- function(input, method) {
  if (!method %in% c("limma", "ttest", "lm")) return(list(input = input, note = NULL))
  at <- input$assay_type
  if (identical(at, "raw_count")) {
    m <- input$expr_mat
    lib <- colSums(m, na.rm = TRUE)
    nf <- rep(1, ncol(m))
    if (is_installed("edgeR")) {
      nf <- tryCatch(cached_tmm(m),
                     error = function(e) rep(1, ncol(m)))
    }
    eff <- lib * nf
    eff[!is.finite(eff) | eff <= 0] <- NA_real_
    # voom's log-CPM: half a read added to each count and one read to
    # each library. log2(CPM + 0.5) added half a CPM -- ten reads at 20M
    # depth -- and shrank a true 2-fold change on a 10-read gene to 1.5.
    input$expr_mat <- log2(sweep(m + 0.5, 2L, (eff + 1) / 1e6, "/"))
    input$assay_type <- "logcpm"
    return(list(input = input, note = sprintf(
      "Raw counts were converted to log2-CPM%s for method = '%s'; DESeq2 or edgeR model counts directly.",
      if (any(nf != 1)) " (TMM-scaled)" else "", method)))
  }
  if (isTRUE(at %in% c("raw_intensity", "tpm", "fpkm"))) {
    input$expr_mat <- log2(pmax(input$expr_mat, 0) + 1)
    input$assay_type <- if (identical(input$omics_type, "rnaseq")) "logcpm"
                        else "normalized_intensity"
    return(list(input = input, note = sprintf(
      "`%s` values were log2-transformed for method = '%s', so effects are log2 fold changes.",
      at, method)))
  }
  list(input = input, note = NULL)
}

# The values, checked before any engine sees them.
#
# A count engine given a sample with no counts answered "missing value
# where TRUE/FALSE needed" (edgeR) or "every gene contains at least one
# zero" (DESeq2), and one given a count past the integer range answered
# with an invalid-object error after a coercion warning. Neither names
# the sample or the value. The intensity engines take an infinite value
# as it comes and hand back an infinite effect with no p-value for the
# feature; here it becomes a missing value, and the bundle says so.
preflight_diff_matrix <- function(input, method, analysis_type, args,
                                  specs = NULL) {
  mat <- input$expr_mat
  samples <- colnames(mat)
  if (analysis_type == "group") {
    g <- input$meta_df[[args$group_col]]
    groups <- if (!is.null(specs)) contrast_levels(specs)
              else c(args$control_group, args$case_group)
    in_contrast <- rownames(input$meta_df)[!is.na(g) & g %in% groups]
    samples <- intersect(samples, in_contrast)
  }
  sub <- mat[, samples, drop = FALSE]

  if (method %in% c("deseq2", "edger")) {
    n_na <- sum(is.na(sub))
    n_inf <- sum(is.infinite(sub))
    n_neg <- sum(sub < 0, na.rm = TRUE)
    if (n_na + n_inf + n_neg > 0L) {
      stop(sprintf(
        "Counts must be finite and non-negative for method = '%s'; the matrix has %d missing, %d infinite and %d negative value(s).",
        method, n_na, n_inf, n_neg), call. = FALSE)
    }
    lib <- colSums(sub)
    empty <- names(lib)[lib == 0]
    if (length(empty) > 0L) {
      stop(sprintf(
        "%d sample(s) have no counts at all: %s. Remove them before running %s.",
        length(empty), paste(utils::head(empty, 5L), collapse = ", "), method),
        call. = FALSE)
    }
    if (method == "deseq2") {
      n_big <- sum(sub > .Machine$integer.max)
      if (n_big > 0L) {
        stop(sprintf(
          "%d count(s) exceed %d, the largest value DESeq2 can store; check the units of the matrix.",
          n_big, .Machine$integer.max), call. = FALSE)
      }
    }
    return(list(input = input, warnings = character(0)))
  }

  nonfinite <- is.infinite(mat)
  if (!any(nonfinite)) return(list(input = input, warnings = character(0)))
  note <- sprintf(
    "%d infinite value(s) in %d feature(s) were treated as missing.",
    sum(nonfinite), sum(rowSums(nonfinite) > 0L))
  input$expr_mat[nonfinite] <- NA
  warning(note, call. = FALSE)
  list(input = input, warnings = note)
}

validate_diff_args <- function(analysis_type, method, args) {
  if (analysis_type %in% c("group", "anova")) {
    if (is.null(args$group_col)) {
      stop("`group_col` is required for analysis_type = '", analysis_type, "'.")
    }
  }
  if (analysis_type == "group" && is.null(args$contrasts)) {
    if (is.null(args$control_group) || is.null(args$case_group)) {
      stop("`control_group` and `case_group` (or `contrasts`) are required for analysis_type = 'group'.")
    }
  }
  if (analysis_type == "continuous") {
    if (is.null(args$continuous_col)) {
      stop("`continuous_col` is required for analysis_type = 'continuous'.")
    }
  }
  if (analysis_type == "anova" && !method %in% c("limma", "edger", "deseq2")) {
    stop("ANOVA analysis_type needs method = 'limma' (continuous data) or ",
         "'edger' / 'deseq2' (raw counts).")
  }
  if (method == "edger" && !analysis_type %in% c("group", "anova")) {
    stop("edgeR backend currently only supports analysis_type = 'group' or 'anova'.")
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

# The fitted model kept in the bundle. DESeq2's DESeqDataSet and edgeR's
# fit are dropped: they were 56-314 MB per bundle at 30k-60k genes,
# nothing reads them, and every bundle is sent back from the worker and
# written by every autosave (a project of two layers went from 69 to
# 18 MB without them). limma's fit is a few MB and stays.
keep_model_object <- function(obj) {
  heavy <- inherits(obj, c("DESeqDataSet", "DGEGLM", "DGELRT", "DGEList")) ||
    (is.list(obj) && !is.object(obj) && length(obj) &&
       all(vapply(obj, inherits, logical(1), c("DESeqDataSet", "DGEGLM", "DGELRT", "DGEList"))))
  if (heavy && !isTRUE(getOption("omicsCore.keep_count_models", FALSE))) NULL else obj
}
