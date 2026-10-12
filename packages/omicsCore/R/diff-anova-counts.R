# "Does this feature differ between any of the groups?" for raw counts.
#
# The limma ANOVA (`run_limma_anova()`) only takes continuous data, so a
# counts layer with several groups had no global test at all. edgeR asks
# it as a quasi-likelihood F-test on every group coefficient at once;
# DESeq2 as a likelihood-ratio test of the model with the groups against
# the model without them. Both return one row per feature with no
# direction -- which groups differ is what the pairwise contrasts are for.
# Both normalise as the pairwise fits do, with tximport's gene-length
# offsets when the input carries them.

# Shared set-up: the samples of the groups tested, a design with the
# groups after any block and before any covariates, all checked.
anova_counts_design <- function(input, group_col, covariates, selected_groups,
                                paired_col, fn) {
  validate_omics_input(input)
  if (input$omics_type != "rnaseq" || !identical(input$assay_type, "raw_count")) {
    stop("`", fn, "()` requires RNA-seq raw counts.", call. = FALSE)
  }
  meta_df <- input$meta_df
  if (!group_col %in% colnames(meta_df)) {
    stop("`group_col` not found in `meta_df`: ", group_col)
  }
  check_paired_col(meta_df, paired_col, object_name = "meta_df")
  target <- meta_df[!is.na(meta_df[[group_col]]), , drop = FALSE]
  if (!is.null(selected_groups)) {
    target <- target[target[[group_col]] %in% selected_groups, , drop = FALSE]
  }
  target[[group_col]] <- factor(as.character(target[[group_col]]))
  if (nlevels(target[[group_col]]) < 2L) {
    stop("ANOVA needs at least two groups in `", group_col, "`.", call. = FALSE)
  }
  mt <- count_model_terms(target, group_col, paired_col, covariates)
  target <- mt$meta
  red_terms <- setdiff(mt$terms, group_col)
  list(target = target, counts = as.matrix(input$expr_mat)[, rownames(target), drop = FALSE],
       full = backtick_formula(mt$terms), reduced = backtick_formula(red_terms),
       terms = mt$terms)
}

#' edgeR global test across groups
#'
#' Quasi-likelihood F-test on all group coefficients together. Internal --
#' call via `run_diff(method = "edger", analysis_type = "anova")`.
#'
#' @inheritParams run_limma_anova
#' @param prefilter Whether to set aside genes with too few counts before
#'   fitting, with `edgeR::filterByExpr()`. On by default; see [run_diff()].
#' @return List with `results_raw`, `results_std`, `model_object`, and
#'   `analysis_info`.
#' @keywords internal
run_edger_anova <- function(input, group_col, covariates = NULL,
                            selected_groups = NULL, paired_col = NULL,
                            prefilter = TRUE) {
  ensure_edger()
  d <- anova_counts_design(input, group_col, covariates, selected_groups,
                           paired_col, "run_edger_anova")
  design <- stats::model.matrix(d$full, data = d$target)
  # The group's columns by term, not by a name pattern: a covariate
  # called "group_batch" matched "^group" and was tested along with the
  # groups.
  term <- match(group_col, d$terms)
  grp_cols <- which(attr(design, "assign") == term)
  # With tximport's gene lengths when the input has them, as the
  # pairwise fit uses (edger_ql_fit(), diff-edger.R).
  ef <- edger_ql_fit(d$counts, design, input, prefilter = prefilter)
  fit <- ef$fit
  filter_note <- ef$filter_note
  report_progress("Testing (4 of 4)")
  qlf <- edgeR::glmQLFTest(fit, coef = grp_cols)
  raw_df <- as.data.frame(edgeR::topTags(qlf, n = Inf, sort.by = "none")$table) |>
    tibble::rownames_to_column("feature_id")
  raw_df <- pad_untested(raw_df, rownames(d$counts))
  std <- standardize_global_test_results(
    raw_df, input$feature_df, group_col, input$omics_type, method = "edger",
    stat = "F", stat_type = "F", p = "PValue", padj = "FDR", base_mean = "logCPM")
  list(results_raw = raw_df, results_std = std, model_object = fit,
       analysis_info = list(omics_type = input$omics_type, method = "edger",
                            analysis_type = "anova", comparison = group_col,
                            covariates = covariates,
                            selected_groups = selected_groups,
                            paired_col = paired_col,
                            warnings = filter_note))
}

#' DESeq2 global test across groups
#'
#' Likelihood-ratio test of the model with the groups against the model
#' without them. Internal -- call via
#' `run_diff(method = "deseq2", analysis_type = "anova")`.
#'
#' @inheritParams run_limma_anova
#' @param prefilter Whether to set aside genes with too few reads before
#'   fitting (DESeq2's recommended pre-filter: at least 10 reads in at least
#'   as many samples as the smallest group). Off by default; see [run_diff()].
#' @return List with `results_raw`, `results_std`, `model_object`, and
#'   `analysis_info`.
#' @keywords internal
run_deseq2_anova <- function(input, group_col, covariates = NULL,
                             selected_groups = NULL, paired_col = NULL,
                             prefilter = FALSE) {
  ensure_deseq2()
  d <- anova_counts_design(input, group_col, covariates, selected_groups,
                           paired_col, "run_deseq2_anova")
  safe <- deseq2_safe_coldata(d$target, d$terms)
  red <- setdiff(unlist(safe$map), safe$map[[group_col]])
  reduced <- stats::as.formula(paste("~", if (length(red)) paste(red, collapse = " + ") else "1"))
  # The pairwise fit's dataset builder, so tximport's gene lengths become
  # normalisation factors here too (DESeqDataSetFromTximport).
  pf <- deseq2_prefilter(d$counts, safe, prefilter)
  dds <- build_deseq_dataset(input, pf$counts, safe$col_data, safe$formula)
  bp <- deseq2_bpparam()
  dds <- with_fixed_seed(1L, with_deseq2_progress(
    if (is.null(bp)) DESeq2::DESeq(dds, test = "LRT", reduced = reduced, quiet = FALSE)
    else DESeq2::DESeq(dds, test = "LRT", reduced = reduced, quiet = FALSE,
                       parallel = TRUE, BPPARAM = bp)))
  raw_df <- as.data.frame(DESeq2::results(dds)) |>
    tibble::rownames_to_column("feature_id")
  raw_df <- pad_untested(raw_df, rownames(d$counts))
  std <- standardize_global_test_results(
    raw_df, input$feature_df, group_col, input$omics_type, method = "deseq2",
    stat = "stat", stat_type = "LRT", p = "pvalue", padj = "padj", base_mean = "baseMean")
  list(results_raw = raw_df, results_std = std, model_object = dds,
       analysis_info = list(omics_type = input$omics_type, method = "deseq2",
                            analysis_type = "anova", comparison = group_col,
                            covariates = covariates,
                            selected_groups = selected_groups,
                            paired_col = paired_col,
                            warnings = pf$note))
}
