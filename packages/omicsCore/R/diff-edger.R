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

  keep_samples <- rownames(target_meta)
  count_sub <- count_mat[, keep_samples, drop = FALSE]

  design_terms <- c(group_col)
  if (!is.null(paired_col)) {
    if (any(is.na(target_meta[[paired_col]]))) {
      stop("`paired_col` contains missing values after group filtering: ", paired_col)
    }
    target_meta[[paired_col]] <- factor(target_meta[[paired_col]])
    design_terms <- c(paired_col, design_terms)
  }
  if (!is.null(covariates)) {
    missing_cov <- setdiff(covariates, colnames(target_meta))
    if (length(missing_cov) > 0L) {
      stop("Missing covariates: ", paste(missing_cov, collapse = ", "))
    }
    design_terms <- c(design_terms, covariates)
  }

  # Backticked: a column called "Treatment Group" was an "unexpected
  # symbol" in the formula.
  design_formula <- stats::as.formula(paste("~", paste0("`", design_terms, "`", collapse = " + ")))
  design_mat <- stats::model.matrix(design_formula, data = target_meta)
  # The group's columns, found by which term they belong to rather than
  # by name: a non-syntactic column name puts backticks into the column
  # names, and a covariate whose name starts with the group column's
  # matched the old pattern.
  grp_term <- match(group_col, design_terms)
  grp_cols <- which(attr(design_mat, "assign") == grp_term)
  names(grp_cols) <- levels(target_meta[[group_col]])[-1L]

  all_features <- rownames(count_sub)
  y <- edgeR::DGEList(counts = count_sub)
  y <- edger_filter(y, design_mat)
  filter_note <- attr(y, "filter_note")
  report_progress("Normalising library sizes (1 of 4)")
  y <- edger_normalise(y, input)

  report_progress("Estimating dispersions (2 of 4)")
  y <- edgeR::estimateDisp(y, design = design_mat)
  report_progress("Fitting the model (3 of 4)")
  fit <- edgeR::glmQLFit(y, design = design_mat)

  # model.matrix() names a factor's columns `<column><level>` verbatim,
  # so the coefficients are found by exact name. They used to be found by
  # substring, which picked "B" out of "groupAB" when both were present.
  # The design is treatment coding against the first level, so each
  # other level's column is its difference from it and a contrast is its
  # weights placed on those columns.
  ref <- group_levels[[1L]]
  comparisons <- vapply(specs, `[[`, character(1), "label")
  per <- lapply(seq_along(specs), function(i) {
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
    raw_df <- pad_untested(raw_df, all_features)
    std <- standardize_edger_group_results(
      raw_df = raw_df,
      feature_df = feature_df,
      comparison = comparisons[[i]],
      omics_type = input$omics_type
    )
    raw_df$comparison <- comparisons[[i]]
    list(raw = raw_df, std = std)
  })
  raw_df <- do.call(rbind, lapply(per, `[[`, "raw"))
  if (length(specs) == 1L) raw_df$comparison <- NULL
  results_std <- do.call(rbind, lapply(per, `[[`, "std"))
  rownames(raw_df) <- NULL
  rownames(results_std) <- NULL
  comparison <- comparisons

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
