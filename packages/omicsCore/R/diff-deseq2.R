# DESeq2 backend, Suggests-gated. Hooks into an optional tximport length
# matrix at input$misc$tximport$length (so downstream Salmon-derived counts
# can carry gene-length offsets) without making tximport a hard dep.

ensure_deseq2 <- function() {
  if (!is_installed("DESeq2")) {
    stop(
      "Package 'DESeq2' is required for the DESeq2 differential backend. ",
      "Install with: omicsCore::install_optional('rnaseq').",
      call. = FALSE
    )
  }
}

# Extract tximport metadata stored on an omics_input. Returns a list with
# `length` (matrix or NULL), `counts_from_abundance` (one of "no",
# "scaledTPM", "lengthScaledTPM", or NULL), and `has_tximport` (logical).
get_tximport_info <- function(input) {
  misc <- input$misc
  if (!is.list(misc) || !is.list(misc$tximport)) {
    return(list(length = NULL, counts_from_abundance = NULL, has_tximport = FALSE))
  }
  txi <- misc$tximport
  length_mat <- if (!is.null(txi$length)) as.matrix(txi$length) else NULL
  cfa <- txi$counts_from_abundance
  if (!is.null(cfa)) cfa <- as.character(cfa)[[1]]
  if (is.null(cfa) && !is.null(length_mat)) cfa <- "no"
  list(
    length = length_mat,
    counts_from_abundance = cfa,
    has_tximport = !is.null(length_mat) || !is.null(cfa)
  )
}

build_deseq_dataset <- function(input, count_mat, col_data, design) {
  count_mat <- as.matrix(count_mat)
  storage.mode(count_mat) <- "numeric"
  if (any(is.na(count_mat))) {
    stop("DESeq2 count input contains missing values.")
  }
  if (any(!is.finite(count_mat))) {
    stop("DESeq2 count input must be finite.")
  }
  if (any(count_mat < 0)) {
    stop("DESeq2 count input must be non-negative.")
  }

  txi_info <- get_tximport_info(input)
  allowed_cfa <- c("no", "scaledTPM", "lengthScaledTPM")
  if (!is.null(txi_info$counts_from_abundance) &&
      !txi_info$counts_from_abundance %in% allowed_cfa) {
    stop(
      "Unsupported `counts_from_abundance`: ", txi_info$counts_from_abundance,
      ". Expected one of: ", paste(allowed_cfa, collapse = ", ")
    )
  }

  if (isTRUE(txi_info$has_tximport) && !is.null(txi_info$length)) {
    length_mat <- txi_info$length
    if (is.null(rownames(length_mat)) || is.null(colnames(length_mat))) {
      stop("Stored tximport length matrix must have row and column names.")
    }
    missing_feat <- setdiff(rownames(count_mat), rownames(length_mat))
    missing_samp <- setdiff(colnames(count_mat), colnames(length_mat))
    if (length(missing_feat) > 0L || length(missing_samp) > 0L) {
      stop("Stored tximport length matrix is not aligned to the requested subset.")
    }
    txi <- list(
      counts = count_mat,
      length = length_mat[rownames(count_mat), colnames(count_mat), drop = FALSE],
      countsFromAbundance = if (is.null(txi_info$counts_from_abundance)) "no" else txi_info$counts_from_abundance
    )
    return(DESeq2::DESeqDataSetFromTximport(txi = txi, colData = col_data, design = design))
  }

  if (!all(abs(count_mat - round(count_mat)) < .Machine$double.eps^0.5, na.rm = TRUE)) {
    stop("DESeq2 requires integer-like raw count values when no tximport metadata is available.")
  }
  DESeq2::DESeqDataSetFromMatrix(
    countData = round(count_mat),
    colData = col_data,
    design = design
  )
}

# DESeq2 draws from the random stream while it fits, so a differential
# run used to leave the caller's stream somewhere else -- and in a fresh
# session, to create a `.Random.seed` that was not there before. The
# stream is pinned for the fit and put back.
deseq_with_dispersion_fallback <- function(dds) {
  with_fixed_seed(1L, tryCatch(
    DESeq2::DESeq(dds, quiet = TRUE),
    error = function(e) {
      msg <- conditionMessage(e)
      if (!grepl("dispersion", msg, ignore.case = TRUE)) stop(e)
      # Toy / tiny datasets occasionally fail the default dispersion fit.
      # Gene-wise dispersions without shrinkage are anti-conservative, so
      # the fallback is said out loud rather than taken quietly.
      warning("DESeq2's dispersion trend could not be fitted (", msg,
              "); gene-wise dispersion estimates were used, which makes ",
              "p-values optimistic.", call. = FALSE)
      dds_fb <- DESeq2::estimateSizeFactors(dds)
      dds_fb <- DESeq2::estimateDispersionsGeneEst(dds_fb, quiet = TRUE)
      DESeq2::dispersions(dds_fb) <- S4Vectors::mcols(dds_fb)$dispGeneEst
      DESeq2::nbinomWaldTest(dds_fb, quiet = TRUE)
    }
  ))
}

#' DESeq2 two-group differential test
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
#' @return List with `results_raw`, `results_std`, `model_object` (`DESeqDataSet`),
#'   and `analysis_info`.
#' @keywords internal
run_deseq2_group <- function(
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
    stop("`run_deseq2_group()` requires `omics_type = 'rnaseq'`.")
  }
  if (!identical(input$assay_type, "raw_count")) {
    stop("`run_deseq2_group()` requires `assay_type = 'raw_count'`.")
  }
  ensure_deseq2()

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

  safe <- deseq2_safe_coldata(target_meta, design_terms)
  dds <- build_deseq_dataset(input, count_sub, safe$col_data, safe$formula)
  ref <- group_levels[[1L]]
  gvar <- safe$map[[group_col]]
  dds[[gvar]] <- stats::relevel(dds[[gvar]], ref = ref)
  dds <- deseq_with_dispersion_fallback(dds)

  # One fit, one dispersion estimate over every group in the model; each
  # contrast is then read off it. A pair goes in by name, which DESeq2
  # resolves for any two levels; a weighted contrast goes in as a numeric
  # vector over the model's coefficients.
  comparisons <- vapply(specs, `[[`, character(1), "label")
  per <- lapply(seq_along(specs), function(i) {
    s <- specs[[i]]
    res <- if (!is.null(s$case)) {
      DESeq2::results(dds, contrast = c(gvar, s$case, s$control))
    } else {
      DESeq2::results(dds, contrast = deseq2_contrast_vector(dds, gvar, s$weights, ref))
    }
    raw_df <- as.data.frame(res) |>
      tibble::rownames_to_column("feature_id")
    std <- standardize_deseq2_group_results(
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
    model_object = dds,
    analysis_info = list(
      omics_type = input$omics_type,
      method = "deseq2",
      analysis_type = "group",
      comparison = comparison,
      covariates = covariates,
      paired_col = paired_col
    )
  )
}

# The design columns under placeholder names. DESeq2 rewrites colData
# names that are not R names ("Treatment Group" -> "Treatment.Group"), so
# a design formula naming the original column could not find it.
deseq2_safe_coldata <- function(meta, terms) {
  map <- stats::setNames(paste0("v", seq_along(terms)), terms)
  col_data <- meta[, terms, drop = FALSE]
  names(col_data) <- unname(map)
  rownames(col_data) <- rownames(meta)
  list(col_data = col_data, map = as.list(map),
       formula = stats::as.formula(paste("~", paste(map, collapse = " + "))))
}

# A contrast's weights as DESeq2's numeric contrast. The model is in
# treatment coding against `ref`, so each other level's coefficient is its
# difference from the reference and the reference's own weight drops out
# (the weights sum to zero). DESeq2 renames levels that are not plain
# words ("Drug A" -> "Drug.A"), so the coefficient is found by matching
# either spelling.
deseq2_contrast_vector <- function(dds, group_col, weights, ref) {
  rn <- DESeq2::resultsNames(dds)
  clean <- function(x) gsub("[^[:alnum:]_.]", ".", x)
  v <- stats::setNames(rep(0, length(rn)), rn)
  for (lv in setdiff(names(weights), ref)) {
    if (abs(weights[[lv]]) < 1e-12) next
    cand <- unique(c(
      paste0(group_col, "_", lv, "_vs_", ref),
      paste0(clean(group_col), "_", clean(lv), "_vs_", clean(ref)),
      paste0(make.names(group_col), "_", make.names(lv), "_vs_", make.names(ref))
    ))
    hit <- intersect(cand, rn)
    if (!length(hit)) {
      stop("Could not find the DESeq2 coefficient for group '", lv, "'.",
           call. = FALSE)
    }
    v[[hit[[1L]]]] <- weights[[lv]]
  }
  unname(v)
}

#' DESeq2 continuous-variable differential test
#'
#' @param input A validated RNA-seq `omics_input` built from raw counts.
#' @param continuous_col Continuous metadata column.
#' @param covariates Optional covariate column names.
#' @param paired_col Optional pairing column.
#'
#' @return List with `results_raw`, `results_std`, `model_object` (`DESeqDataSet`),
#'   and `analysis_info`.
#' @keywords internal
run_deseq2_continuous <- function(
  input,
  continuous_col,
  covariates = NULL,
  paired_col = NULL
) {
  validate_omics_input(input)
  if (input$omics_type != "rnaseq") {
    stop("`run_deseq2_continuous()` requires `omics_type = 'rnaseq'`.")
  }
  if (!identical(input$assay_type, "raw_count")) {
    stop("`run_deseq2_continuous()` requires `assay_type = 'raw_count'`.")
  }
  ensure_deseq2()

  count_mat <- as.matrix(input$expr_mat)
  meta_df <- input$meta_df
  feature_df <- input$feature_df

  if (!continuous_col %in% colnames(meta_df)) {
    stop("`continuous_col` not found in `meta_df`: ", continuous_col)
  }
  check_paired_col(meta_df, paired_col, object_name = "meta_df")
  validate_continuous_pairing(meta_df, paired_col, object_name = "meta_df")

  cont_vals <- coerce_continuous_col(meta_df[[continuous_col]], continuous_col)
  meta_df[[continuous_col]] <- cont_vals

  design_terms <- c(continuous_col)
  if (!is.null(paired_col)) {
    if (any(is.na(meta_df[[paired_col]]))) {
      stop("`paired_col` contains missing values: ", paired_col)
    }
    meta_df[[paired_col]] <- factor(meta_df[[paired_col]])
    design_terms <- c(paired_col, design_terms)
  }
  if (!is.null(covariates)) {
    missing_cov <- setdiff(covariates, colnames(meta_df))
    if (length(missing_cov) > 0L) {
      stop("Missing covariates: ", paste(missing_cov, collapse = ", "))
    }
    design_terms <- c(design_terms, covariates)
  }

  safe <- deseq2_safe_coldata(meta_df, design_terms)
  dds <- build_deseq_dataset(input, count_mat, safe$col_data, safe$formula)
  dds <- deseq_with_dispersion_fallback(dds)

  coef_names <- DESeq2::resultsNames(dds)
  cvar <- safe$map[[continuous_col]]
  if (!cvar %in% coef_names) {
    stop(
      "Failed to resolve DESeq2 coefficient for `continuous_col = '",
      continuous_col, "'`. Available coefficients: ",
      paste(coef_names, collapse = ", ")
    )
  }

  res <- DESeq2::results(dds, name = cvar)
  raw_df <- as.data.frame(res) |>
    tibble::rownames_to_column("feature_id")

  results_std <- standardize_deseq2_continuous_results(
    raw_df = raw_df,
    feature_df = feature_df,
    comparison = continuous_col,
    omics_type = input$omics_type
  )

  list(
    results_raw = raw_df,
    results_std = results_std,
    model_object = dds,
    analysis_info = list(
      omics_type = input$omics_type,
      method = "deseq2",
      analysis_type = "continuous_linear",
      comparison = continuous_col,
      covariates = covariates,
      paired_col = paired_col
    )
  )
}
