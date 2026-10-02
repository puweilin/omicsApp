#' Run the standard quality-control pipeline
#'
#' End-to-end QC for an [omics_input()]. Computes missingness, detects
#' sample-level outliers, optionally imputes the expression matrix, and
#' returns an [analysis_bundle][is_analysis_bundle()] containing the cleaned
#' input plus the QC summary so the result can be plotted by [plot_qc()] or
#' consumed directly by downstream analysis functions.
#'
#' Sensible defaults:
#'
#' * `omics_type == "proteomics"` → `outlier_method = "pca"`,
#'   `impute_method = "MinProb"`. Missingness in DIA/DDA is mostly
#'   left-censored -- a protein is absent because it fell below the
#'   detection limit -- and leaving `NA` is not the neutral choice it
#'   looks like: limma drops what it cannot fit, so "none" is
#'   complete-case analysis taken silently.
#' * `omics_type == "rnaseq"` → `outlier_method = "connectivity"`,
#'   `impute_method = "none"`. A zero count is an observation, and
#'   imputing it feeds a negative-binomial model numbers it never saw.
#'
#' Pass explicit arguments to override the defaults.
#'
#' Samples flagged as outliers are reported, not removed, unless
#' `remove_outliers = TRUE`: dropping a sample changes the design -- a
#' group of three becomes a group of two -- and that is the analyst's
#' call, made after looking at the PCA. Samples above an explicit
#' `sample_missing_threshold` are removed, since that threshold was asked
#' for.
#'
#' Imputation on a linear scale (raw intensities) runs on log2 values and
#' is transformed back: the left-censored methods draw from a normal
#' distribution, which only fits intensities after logging, and on the
#' raw scale they returned negative intensities.
#'
#' @param input An `omics_input`.
#' @param missing_threshold Feature missing-rate cutoff in `[0, 1]`. Features
#'   above this are flagged and removed from `cleaned_input`. Default `0.5`.
#' @param sample_missing_threshold Optional sample missing-rate cutoff.
#'   Samples above this are flagged and removed.
#' @param impute_method One of [IMPUTE_METHODS] -- DEP's method set.
#'   Applied to the expression matrix of `cleaned_input` after filtering.
#'   `NULL` (the default) resolves per modality via
#'   [resolve_impute_method()]: `"MinProb"` for proteomics, `"none"` for
#'   counts.
#' @param outlier_method One of `"none"`, `"pca"`, `"connectivity"`, `"iqr"`,
#'   or a vector of those (other than `"none"`) to union their flags.
#' @param outlier_sd_threshold Z-score / IQR multiplier passed to
#'   [qc_outliers()]. Default `3`.
#' @param remove_outliers Remove the samples [qc_outliers()] flags from
#'   `cleaned_input`. Default `FALSE`: they are listed in
#'   `recommended_filters$remove_samples` and kept.
#' @param ... Forwarded to the imputation backend.
#'
#' @return An `analysis_bundle` with the following fields under `results`:
#'   \describe{
#'     \item{`qc_summary`}{List with `missingness`, `depth`, `outliers`,
#'       `recommended_filters` (sample/feature IDs to remove), and, when
#'       values were imputed, `imputation` (the method and the positions of
#'       the imputed cells in `cleaned_input$expr_mat`).}
#'     \item{`cleaned_input`}{`omics_input` with flagged features (and
#'       samples, see above) removed and an (optionally) imputed expression
#'       matrix. `raw_mat` carries the pre-imputation matrix when
#'       imputation occurred and the input had none.}
#'   }
#'   Notes about what was done to the data are in `warnings`.
#' @export
#' @family qc
#' @examples
#' set.seed(1)
#' expr <- matrix(rnorm(60), nrow = 6,
#'                dimnames = list(paste0("g", 1:6), paste0("s", 1:10)))
#' expr[, 1] <- NA  # entirely missing sample
#' meta <- data.frame(group = rep(c("A", "B"), each = 5),
#'                    row.names = colnames(expr))
#' feat <- data.frame(feature_id = rownames(expr),
#'                    row.names = rownames(expr))
#' input <- omics_input(expr, meta, feat, omics_type = "proteomics")
#' bundle <- run_qc(input, sample_missing_threshold = 0.9,
#'                  outlier_method = "iqr")
#' bundle
run_qc <- function(
  input,
  missing_threshold = 0.5,
  sample_missing_threshold = NULL,
  impute_method = NULL,
  outlier_method = NULL,
  outlier_sd_threshold = 3,
  remove_outliers = FALSE,
  ...
) {
  assert_number(missing_threshold, "missing_threshold", lower = 0, upper = 1)
  assert_number(sample_missing_threshold, "sample_missing_threshold",
                lower = 0, upper = 1, allow_null = TRUE)
  assert_subset(outlier_method, "outlier_method",
                c("none", "pca", "connectivity", "iqr"), allow_null = TRUE)
  assert_number(outlier_sd_threshold, "outlier_sd_threshold", lower = 0)
  assert_flag(remove_outliers, "remove_outliers")
  validate_omics_input(input)

  # NULL rather than a fixed default, resolved per modality like
  # outlier_method below. Proteomics gets MinProb because its
  # missingness is left-censored; counts get "none" because a zero is an
  # observation. One global default would be wrong for one of them
  # whichever way it went.
  if (is.null(impute_method)) {
    impute_method <- resolve_impute_method(input$omics_type)
  }
  impute_method <- match.arg(impute_method, IMPUTE_METHODS)

  # Resolve outlier defaults per omics_type.
  if (is.null(outlier_method)) {
    outlier_method <- switch(
      input$omics_type %||% "",
      proteomics = "pca",
      rnaseq     = "connectivity",
      "pca"
    )
  }

  # ---- missingness ----
  missingness <- qc_missingness(
    input,
    sample_missing_cutoff = sample_missing_threshold,
    feature_missing_cutoff = missing_threshold
  )

  # ---- depth ----
  # Always, for both modalities: it is one pass over the matrix, and
  # deciding in advance which panel someone will want is how the RNA-seq
  # view ended up with nothing to show.
  depth <- qc_depth(input)

  # ---- outliers ----
  run_outliers <- !identical(outlier_method, "none") &&
                  !(length(outlier_method) == 1L && is.na(outlier_method))
  outliers <- if (run_outliers) {
    qc_outliers(
      input,
      method = outlier_method,
      sd_threshold = outlier_sd_threshold
    )
  } else {
    list(method = "none", stats = data.frame(), flagged_samples = character(0))
  }

  remove_samples <- unique(c(
    missingness$flagged_samples,
    outliers$flagged_samples
  ))
  remove_features <- unique(missingness$flagged_features)
  notes <- outliers$note

  # ---- build cleaned input ----
  drop_samples <- unique(c(missingness$flagged_samples,
                           if (remove_outliers) outliers$flagged_samples))
  kept_outliers <- setdiff(outliers$flagged_samples, drop_samples)
  if (length(kept_outliers)) {
    notes <- c(notes, sprintf(
      "Flagged as outlier%s and kept: %s. Look at the PCA, then exclude with remove_outliers = TRUE or subset_omics() if it is not biology.",
      if (length(kept_outliers) > 1L) "s" else "",
      paste(kept_outliers, collapse = ", ")))
  }
  keep_samples <- setdiff(colnames(input$expr_mat), drop_samples)
  keep_features <- setdiff(rownames(input$expr_mat), remove_features)
  if (length(keep_samples) == 0L || length(keep_features) == 0L) {
    stop("QC would remove all samples or features; loosen the thresholds.")
  }

  cleaned <- subset_omics(input, samples = keep_samples, features = keep_features)

  # ---- imputation ----
  imputation <- NULL
  if (impute_method != "none" && anyNA(cleaned$expr_mat)) {
    na_cells <- which(is.na(cleaned$expr_mat))
    cleaned$raw_mat <- cleaned$raw_mat %||% cleaned$expr_mat
    linear <- !is.null(cleaned$assay_type) &&
      !cleaned$assay_type %in% LOG_SCALE_ASSAY_TYPES &&
      !impute_method %in% c("zero", "min")
    if (linear) {
      m <- cleaned$expr_mat
      m[m <= 0] <- NA_real_
      imputed <- 2^impute_matrix(log2(m), method = impute_method, ...)
      # Only the cells that were missing take the imputed value; a zero
      # that was observed stays zero.
      out <- cleaned$expr_mat
      out[na_cells] <- imputed[na_cells]
      cleaned$expr_mat <- out
      notes <- c(notes, sprintf(
        "`%s` values were imputed on a log2 scale and transformed back.",
        cleaned$assay_type))
    } else {
      cleaned$expr_mat <- impute_matrix(cleaned$expr_mat, method = impute_method, ...)
      if (identical(cleaned$omics_type, "proteomics") &&
          isTRUE(cleaned$assay_type %in% c("normalized_intensity", "filtered_intensity"))) {
        cleaned$assay_type <- "imputed_intensity"
      }
    }
    imputation <- list(method = impute_method, imputed_cells = na_cells,
                       n_imputed = length(na_cells))
  }

  bundle <- new_analysis_bundle(
    analysis_name = "run_qc",
    input_info = list(
      omics_type = input$omics_type,
      assay_type = input$assay_type,
      n_samples_in = ncol(input$expr_mat),
      n_features_in = nrow(input$expr_mat),
      n_samples_out = ncol(cleaned$expr_mat),
      n_features_out = nrow(cleaned$expr_mat)
    ),
    params = list(
      missing_threshold = missing_threshold,
      sample_missing_threshold = sample_missing_threshold,
      impute_method = impute_method,
      outlier_method = outlier_method,
      outlier_sd_threshold = outlier_sd_threshold,
      remove_outliers = remove_outliers
    ),
    results = list(
      qc_summary = list(
        missingness = missingness,
        depth = depth,
        outliers = outliers,
        recommended_filters = list(
          remove_samples = remove_samples,
          remove_features = remove_features
        ),
        imputation = imputation
      ),
      cleaned_input = cleaned
    ),
    warnings = as.character(notes)
  )
  bundle
}
