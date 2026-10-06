# "Does this feature differ between any of the groups?" for raw counts.
#
# The limma ANOVA (`run_limma_anova()`) only takes continuous data, so a
# counts layer with several groups had no global test at all. edgeR asks
# it as a quasi-likelihood F-test on every group coefficient at once;
# DESeq2 as a likelihood-ratio test of the model with the groups against
# the model without them. Both return one row per feature with no
# direction -- which groups differ is what the pairwise contrasts are for.

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
  block <- character(0)
  if (!is.null(paired_col)) {
    if (anyNA(target[[paired_col]])) {
      stop("`paired_col` contains missing values after group filtering: ", paired_col)
    }
    target[[paired_col]] <- factor(target[[paired_col]])
    block <- paired_col
  }
  missing_cov <- setdiff(covariates, colnames(target))
  if (length(missing_cov)) {
    stop("Missing covariates: ", paste(missing_cov, collapse = ", "))
  }
  bt <- function(x) if (length(x)) paste0("`", x, "`") else character(0)
  full <- stats::as.formula(paste("~", paste(c(bt(block), bt(group_col), bt(covariates)),
                                            collapse = " + ")))
  red_terms <- c(bt(block), bt(covariates))
  reduced <- stats::as.formula(paste("~", if (length(red_terms))
    paste(red_terms, collapse = " + ") else "1"))
  list(target = target, counts = as.matrix(input$expr_mat)[, rownames(target), drop = FALSE],
       full = full, reduced = reduced, terms = c(block, group_col, covariates))
}

standardize_anova_counts <- function(raw_df, feature_df, method, stat, stat_type,
                                     p, padj, base_mean, comparison, omics_type) {
  feature_df <- prep_feature_df_for_standardize(feature_df)
  out <- dplyr::left_join(raw_df, feature_df, by = "feature_id")
  out <- data.frame(
    feature_id = out$feature_id,
    feature_symbol = out$feature_symbol,
    feature_type = out$feature_type,
    omics_type = omics_type,
    method = method,
    analysis_type = "anova",
    comparison = comparison,
    effect = out[[stat]],
    effect_type = paste0(stat_type, "_statistic"),
    statistic = out[[stat]],
    statistic_type = stat_type,
    p_value = out[[p]],
    adj_p_value = out[[padj]],
    direction = "ns",
    base_mean = out[[base_mean]],
    model_fit = NA_real_,
    is_significant = NA,
    stringsAsFactors = FALSE
  )
  check_diff_result_schema(out)
  out
}

#' edgeR global test across groups
#'
#' Quasi-likelihood F-test on all group coefficients together. Internal --
#' call via `run_diff(method = "edger", analysis_type = "anova")`.
#'
#' @inheritParams run_limma_anova
#' @return List with `results_raw`, `results_std`, `model_object`, and
#'   `analysis_info`.
#' @keywords internal
run_edger_anova <- function(input, group_col, covariates = NULL,
                            selected_groups = NULL, paired_col = NULL) {
  ensure_edger()
  d <- anova_counts_design(input, group_col, covariates, selected_groups,
                           paired_col, "run_edger_anova")
  design <- stats::model.matrix(d$full, data = d$target)
  # The group's columns by term, not by a name pattern: a covariate
  # called "group_batch" matched "^group" and was tested along with the
  # groups.
  term <- match(group_col, d$terms)
  grp_cols <- which(attr(design, "assign") == term)
  y <- edger_filter(edgeR::DGEList(counts = d$counts), design)
  filter_note <- attr(y, "filter_note")
  report_progress("Normalising library sizes (1 of 4)")
  y <- edgeR::calcNormFactors(y)
  report_progress("Estimating dispersions (2 of 4)")
  y <- edgeR::estimateDisp(y, design = design)
  report_progress("Fitting the model (3 of 4)")
  fit <- edgeR::glmQLFit(y, design = design)
  report_progress("Testing (4 of 4)")
  qlf <- edgeR::glmQLFTest(fit, coef = grp_cols)
  raw_df <- as.data.frame(edgeR::topTags(qlf, n = Inf, sort.by = "none")$table) |>
    tibble::rownames_to_column("feature_id")
  raw_df <- pad_untested(raw_df, rownames(d$counts))
  std <- standardize_anova_counts(raw_df, input$feature_df, "edger", "F", "F",
                                  "PValue", "FDR", "logCPM", group_col,
                                  input$omics_type)
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
#' @return List with `results_raw`, `results_std`, `model_object`, and
#'   `analysis_info`.
#' @keywords internal
run_deseq2_anova <- function(input, group_col, covariates = NULL,
                             selected_groups = NULL, paired_col = NULL) {
  ensure_deseq2()
  d <- anova_counts_design(input, group_col, covariates, selected_groups,
                           paired_col, "run_deseq2_anova")
  safe <- deseq2_safe_coldata(d$target, d$terms)
  red <- setdiff(unlist(safe$map), safe$map[[group_col]])
  reduced <- stats::as.formula(paste("~", if (length(red)) paste(red, collapse = " + ") else "1"))
  dds <- build_deseq_dataset(input, d$counts, safe$col_data, safe$formula)
  bp <- deseq2_bpparam()
  dds <- with_fixed_seed(1L, with_deseq2_progress(
    if (is.null(bp)) DESeq2::DESeq(dds, test = "LRT", reduced = reduced, quiet = FALSE)
    else DESeq2::DESeq(dds, test = "LRT", reduced = reduced, quiet = FALSE,
                       parallel = TRUE, BPPARAM = bp)))
  raw_df <- as.data.frame(DESeq2::results(dds)) |>
    tibble::rownames_to_column("feature_id")
  std <- standardize_anova_counts(raw_df, input$feature_df, "deseq2", "stat",
                                  "LRT", "pvalue", "padj", "baseMean", group_col,
                                  input$omics_type)
  list(results_raw = raw_df, results_std = std, model_object = dds,
       analysis_info = list(omics_type = input$omics_type, method = "deseq2",
                            analysis_type = "anova", comparison = group_col,
                            covariates = covariates,
                            selected_groups = selected_groups,
                            paired_col = paired_col))
}
