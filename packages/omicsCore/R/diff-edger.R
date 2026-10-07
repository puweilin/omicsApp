# edgeR QL F-test backend, Suggests-gated. Uses an optional tximport length
# matrix (from input$misc$tximport$length) to scale offsets when available;
# otherwise falls back to standard library-size normalization.

ensure_edger <- function() {
  if (!is_installed("edgeR")) {
    stop(
      "Package 'edgeR' is required for the edgeR differential backend. ",
      "Install with: omicsCore::install_optional('rnaseq').",
      call. = FALSE
    )
  }
}

#' edgeR QL F-test for a two-group comparison
#'
#' @param input A validated RNA-seq `omics_input` built from raw counts.
#' @param group_col Group column in sample metadata.
#' @param control_group Control-group label.
#' @param case_group Case-group label.
#' @param covariates Optional covariate column names.
#' @param paired_col Optional pairing column.
#' @param contrasts Parsed contrast specs (from `run_diff(contrasts = )`); when
#'   given, every group they name is fitted and each contrast read off the fit.
#'
#' @return List with `results_raw`, `results_std`, `model_object`
#'   (`DGEGLM`), and `analysis_info`.
#' @keywords internal
run_edger_group <- function(
  input,
  group_col,
  control_group,
  case_group,
  covariates = NULL,
  paired_col = NULL,
  contrasts = NULL
) {
  validate_omics_input(input)
  if (input$omics_type != "rnaseq") {
    stop("`run_edger_group()` requires `omics_type = 'rnaseq'`.")
  }
  if (!identical(input$assay_type, "raw_count")) {
    stop("`run_edger_group()` requires `assay_type = 'raw_count'`.")
  }
  ensure_edger()

  count_mat <- as.matrix(input$expr_mat)
  meta_df <- input$meta_df
  feature_df <- input$feature_df

  tm <- contrast_target_meta(meta_df, group_col, control_group, case_group,
                             contrasts = contrasts, paired_col = paired_col)
  specs <- tm$specs
  group_levels <- tm$group_levels
  count_sub <- count_mat[, rownames(tm$target_meta), drop = FALSE]
  mt <- count_model_terms(tm$target_meta, group_col, paired_col, covariates)
  target_meta <- mt$meta
  design_terms <- mt$terms

  design_mat <- stats::model.matrix(backtick_formula(design_terms), data = target_meta)
  # The group's columns, found by which term they belong to rather than
  # by name: a non-syntactic column name puts backticks into the column
  # names, and a covariate whose name starts with the group column's
  # matched the old pattern.
  grp_term <- match(group_col, design_terms)
  grp_cols <- which(attr(design_mat, "assign") == grp_term)
  names(grp_cols) <- levels(target_meta[[group_col]])[-1L]

  all_features <- rownames(count_sub)
  ef <- edger_ql_fit(count_sub, design_mat, input)
  fit <- ef$fit
  filter_note <- ef$filter_note

  # model.matrix() names a factor's columns `<column><level>` verbatim,
  # so the coefficients are found by exact name. They used to be found by
  # substring, which picked "B" out of "groupAB" when both were present.
  # The design is treatment coding against the first level, so each
  # other level's column is its difference from it and a contrast is its
  # weights placed on those columns.
  ref <- group_levels[[1L]]
  st <- stack_contrast_results(
    specs,
    raw_for = function(i) {
      w <- specs[[i]]$weights
      cvec <- stats::setNames(rep(0, ncol(design_mat)), colnames(design_mat))
      for (lv in setdiff(names(w), ref)) {
        if (abs(w[[lv]]) < 1e-12) next
        if (!lv %in% names(grp_cols)) {
          stop("Could not locate group coefficient in design matrix for: ", lv)
        }
        cvec[[grp_cols[[lv]]]] <- w[[lv]]
      }
      report_progress(if (length(specs) > 1L) {
        sprintf("Testing comparison %d of %d (4 of 4)", i, length(specs))
      } else "Testing (4 of 4)")
      qlf <- edgeR::glmQLFTest(fit, contrast = unname(cvec))
      tt <- edgeR::topTags(qlf, n = Inf, sort.by = "none")
      raw_df <- as.data.frame(tt$table) |>
        tibble::rownames_to_column("feature_id")
      pad_untested(raw_df, all_features)
    },
    standardize = function(raw_df, comparison) {
      standardize_edger_group_results(
        raw_df = raw_df,
        feature_df = feature_df,
        comparison = comparison,
        omics_type = input$omics_type
      )
    }
  )
  raw_df <- st$raw
  results_std <- st$std
  comparison <- st$comparisons

  list(
    results_raw = raw_df,
    results_std = results_std,
    model_object = fit,
    analysis_info = list(
      omics_type = input$omics_type,
      method = "edger",
      analysis_type = "group",
      comparison = comparison,
      covariates = covariates,
      paired_col = paired_col,
      warnings = filter_note
    )
  )
}


# The quasi-likelihood fit both edgeR tests share (the pairwise
# contrasts and the global test): low-count genes set aside, library
# sizes (with tximport's gene lengths when present), dispersions, fit.
# Returns the fit and the note naming the genes set aside, if any.
edger_ql_fit <- function(counts, design, input) {
  y <- edgeR::DGEList(counts = counts)
  y <- edger_filter(y, design)
  filter_note <- attr(y, "filter_note")
  report_progress("Normalising library sizes (1 of 4)")
  y <- edger_normalise(y, input)
  report_progress("Estimating dispersions (2 of 4)")
  y <- edgeR::estimateDisp(y, design = design)
  report_progress("Fitting the model (3 of 4)")
  list(fit = edgeR::glmQLFit(y, design = design), filter_note = filter_note)
}

# Library sizes for a (filtered) DGEList, with tximport's gene lengths
# when the input carries them. Shared by the pairwise fit and the global
# test, which used to normalise by library size alone: a gene whose
# dominant isoform is four times longer in one group, at the same
# expression, has four times the reads there, and the global test called
# it changed where the pairwise test (with the offsets) did not.
#
# Length offsets only for counts that still carry the length bias.
# Counts from abundance ("scaledTPM", "lengthScaledTPM") have it removed
# already, and offsetting them again corrected twice.
edger_normalise <- function(y, input) {
  txi_info <- get_tximport_info(input)
  if (is.null(txi_info$length) ||
      !identical(txi_info$counts_from_abundance %||% "no", "no")) {
    return(edgeR::calcNormFactors(y))
  }
  counts <- y$counts
  length_sub <- txi_info$length[rownames(counts), colnames(counts), drop = FALSE]
  # The tximport vignette's recipe: length factors centred per gene,
  # then TMM on the length-corrected counts. The offsets used to be
  # log(length) + log(library size) with no TMM at all, so a group in
  # which a tenth of the genes went up 8-fold had nearly every other
  # gene called "down" (1449 of 1774 unchanged genes).
  norm_mat <- length_sub / exp(rowMeans(log(length_sub)))
  norm_cts <- counts / norm_mat
  eff_lib <- edgeR::calcNormFactors(norm_cts) * colSums(norm_cts)
  norm_mat <- log(sweep(norm_mat, 2L, eff_lib, "*"))
  edgeR::scaleOffset(y, offset = norm_mat)
}

# Low-count genes out before the model sees them (edgeR's own
# filterByExpr()). Without it thousands of near-zero genes flattened the
# dispersion trend and widened the multiple-testing burden: on a
# simulated 3-vs-3 set with 100 true changes, edgeR found 1 of them
# unfiltered and 24 filtered. The number removed is recorded; a filter
# that would leave fewer than two genes is not applied.
edger_filter <- function(y, design) {
  keep <- tryCatch(edgeR::filterByExpr(y, design = design),
                   error = function(e) rep(TRUE, nrow(y)))
  if (sum(keep) < 2L || all(keep)) return(y)
  out <- y[keep, , keep.lib.sizes = FALSE]
  attr(out, "filter_note") <- sprintf(
    "%d of %d genes with too few counts to test were set aside (edgeR::filterByExpr).",
    sum(!keep), length(keep))
  out
}

# The genes set aside by the filter, back in the table as untested rows
# (no effect, no p-value) in the matrix's order -- every feature keeps a
# row, as it does from every other engine.
pad_untested <- function(raw_df, all_features) {
  missing <- setdiff(all_features, raw_df$feature_id)
  if (!length(missing)) return(raw_df)
  pad <- raw_df[rep(NA_integer_, length(missing)), , drop = FALSE]
  pad$feature_id <- missing
  out <- rbind(raw_df, pad)
  out <- out[match(all_features, out$feature_id), , drop = FALSE]
  rownames(out) <- NULL
  out
}
