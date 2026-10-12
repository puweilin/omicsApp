#' Run the standard quality-control pipeline
#'
#' End-to-end QC for an [omics_input()]. Computes missingness, detects
#' sample-level outliers, optionally imputes the expression matrix, and
#' returns an [analysis_bundle][is_analysis_bundle()] holding the QC
#' summary and a record of how the input was cleaned, so the result can be
#' plotted by [plot_qc()] and the cleaned input rebuilt by
#' [qc_cleaned_input()] for downstream analysis functions.
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
#' raw scale they returned negative intensities. Only missing cells take
#' an imputed value. When MinProb cannot be estimated (fewer than two
#' features measured in more than half of the samples) MinDet is used
#' instead, and a note says so; see [impute_matrix()].
#'
#' The feature filter is global by default. `missing_filter = "any_group"`
#' applies `missing_threshold` within each group of `group_col` and keeps
#' a feature measured well enough in at least one of them -- the usual
#' proteomics rule, which keeps a protein present in one condition and
#' absent from the other; `"all_groups"` asks for every group. See
#' [qc_missingness()].
#'
#' @param input An `omics_input`.
#' @param missing_threshold Feature missing-rate cutoff in `[0, 1]`. Features
#'   above this are flagged and removed from the cleaned input. Default `0.5`.
#' @param missing_filter Where `missing_threshold` is applied: `"global"`
#'   (default; over all samples), `"any_group"` (within each group; kept
#'   when at least one group passes) or `"all_groups"` (kept when every
#'   group passes).
#' @param group_col Sample-metadata column holding the groups, for the
#'   group filters. `NULL` uses the layer's recorded study design
#'   ([study_design()]). Ignored by `"global"`.
#' @param sample_missing_threshold Optional sample missing-rate cutoff.
#'   Samples above this are flagged and removed.
#' @param impute_method One of [IMPUTE_METHODS] -- DEP's method set.
#'   Applied to the cleaned expression matrix after filtering.
#'   `NULL` (the default) resolves per modality via
#'   [resolve_impute_method()]: `"MinProb"` for proteomics, `"none"` for
#'   counts.
#' @param outlier_method One of `"none"`, `"pca"`, `"connectivity"`, `"iqr"`,
#'   `"loo"`, or a vector of those (other than `"none"`) to union their
#'   flags. `"loo"`, the leave-one-out test of [qc_outliers()], is the one
#'   that can flag a sample in a study of ten or fewer.
#' @param outlier_sd_threshold Z-score / IQR multiplier passed to
#'   [qc_outliers()]. Default `3`.
#' @param remove_outliers Remove the samples [qc_outliers()] flags from
#'   the cleaned input. Default `FALSE`: they are listed in
#'   `recommended_filters$remove_samples` and kept.
#' @param ... Forwarded to the imputation backend.
#'
#' @return An `analysis_bundle` with the following fields under `results`:
#'   \describe{
#'     \item{`qc_summary`}{List with `missingness`, `depth`, `outliers`,
#'       `recommended_filters` (sample/feature IDs to remove), and, when
#'       values were imputed, `imputation` (the method used, the positions
#'       of the imputed cells in the cleaned expression matrix, and
#'       `requested_method` when MinProb fell back to MinDet).}
#'     \item{`cleaning`}{What was done to the input: the samples kept, the
#'       features dropped, the values put into the imputed cells, the
#'       resulting `assay_type`, and a checksum of the input's expression
#'       matrix. [qc_cleaned_input()] rebuilds the cleaned input from it.}
#'     \item{`plot_data`}{What [plot_qc()] draws, computed from the
#'       cleaned input.}
#'   }
#'   The cleaned input is not stored: it is as large as the input, and
#'   [qc_cleaned_input()] rebuilds it exactly -- flagged features (and
#'   samples, see above) removed, the expression matrix (optionally)
#'   imputed, and `raw_mat` holding the pre-imputation matrix when
#'   imputation occurred and the input had none. Results from before this
#'   change carry it as `results$cleaned_input`, and
#'   [qc_cleaned_input()] returns that.
#'
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
  missing_filter = c("global", "any_group", "all_groups"),
  group_col = NULL,
  ...
) {
  assert_number(missing_threshold, "missing_threshold", lower = 0, upper = 1)
  assert_number(sample_missing_threshold, "sample_missing_threshold",
                lower = 0, upper = 1, allow_null = TRUE)
  missing_filter <- match.arg(missing_filter)
  assert_subset(outlier_method, "outlier_method",
                c("none", "pca", "connectivity", "iqr", "loo"), allow_null = TRUE)
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
  report_progress("Counting missing values")
  # The group column is recorded as resolved (from the study design when
  # not given), so the exported script names it rather than depending on
  # a design the reader's copy of the layer may not carry.
  if (!identical(missing_filter, "global")) {
    group_col <- resolve_missing_group_col(input, group_col, missing_filter)
  } else {
    group_col <- NULL
  }
  missingness <- qc_missingness(
    input,
    sample_missing_cutoff = sample_missing_threshold,
    feature_missing_cutoff = missing_threshold,
    missing_filter = missing_filter,
    group_col = group_col
  )

  # ---- depth ----
  # Always, for both modalities: it is one pass over the matrix, and
  # deciding in advance which panel someone will want is how the RNA-seq
  # view ended up with nothing to show.
  depth <- qc_depth(input)

  # ---- outliers ----
  report_progress("Checking for outlier samples")
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
  imputed_values <- NULL
  if (impute_method != "none" && anyNA(cleaned$expr_mat)) {
    report_progress("Imputing missing values")
    na_cells <- which(is.na(cleaned$expr_mat))
    used_method <- impute_method
    # MinProb cannot always be estimated (too few well-measured features);
    # impute_matrix() then uses MinDet and says why. The reason is kept as
    # a note on the result rather than left as a console warning.
    impute <- function(m) withCallingHandlers(
      impute_matrix(m, method = impute_method, ...),
      omics_impute_fallback = function(w) {
        used_method <<- w$used_method %||% "MinDet"
        notes <<- c(notes, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    linear <- !is.null(cleaned$assay_type) &&
      !cleaned$assay_type %in% LOG_SCALE_ASSAY_TYPES &&
      !impute_method %in% c("zero", "min")
    if (linear) {
      m <- cleaned$expr_mat
      m[m <= 0] <- NA_real_
      imputed_values <- (2^impute(log2(m)))[na_cells]
      notes <- c(notes, sprintf(
        "`%s` values were imputed on a log2 scale and transformed back.",
        cleaned$assay_type))
    } else {
      imputed_values <- impute(cleaned$expr_mat)[na_cells]
    }
    imputation <- list(method = used_method, imputed_cells = na_cells,
                       n_imputed = length(na_cells))
    if (!identical(used_method, impute_method)) {
      imputation$requested_method <- impute_method
    }
  }
  # Only the cells that were missing take the imputed value -- a zero
  # that was observed stays zero -- and the matrix is put together the way
  # qc_cleaned_input() puts it together again later, so the two cannot
  # drift apart.
  cleaning <- list(
    version = 1L,
    input_hash = rlang::hash(input$expr_mat),
    kept_samples = keep_samples,
    dropped_features = setdiff(rownames(input$expr_mat), keep_features),
    imputed_values = imputed_values,
    assay_type = qc_cleaned_assay_type(cleaned, imputation)
  )
  cleaned <- qc_apply_cleaning(cleaned, cleaning, imputation)

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
      remove_outliers = remove_outliers,
      missing_filter = missing_filter,
      group_col = group_col
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
      cleaning = cleaning,
      plot_data = qc_plot_data(cleaned, imputation, outliers,
                               outlier_sd_threshold)
    ),
    warnings = as.character(notes)
  )
  bundle
}

# ---- the cleaned input, rebuilt -----------------------------------------
#
# The bundle used to carry the cleaned input itself: a second copy of the
# layer's matrices (expression, raw and normalised, plus feature
# annotation and tximport lengths), 14 MB in memory for 8,000 x 60 and
# 25 MB on real data -- in every saved project and every autosave, for a
# result whose only job downstream is to be looked at. Everything in it
# follows from the layer and a short record of what QC did: which samples
# were kept, which features were dropped, and the values that went into
# the missing cells. So the record is what is kept (`results$cleaning`),
# what the plots draw is computed once (`results$plot_data`), and
# qc_cleaned_input() rebuilds the rest from the layer when it is asked
# for.

#' The cleaned input of a QC result
#'
#' Rebuilds the `omics_input` that [run_qc()] cleaned -- flagged features
#' (and samples) removed, missing values imputed -- from the input it was
#' run on and the record kept in the bundle. The bundle does not carry a
#' copy of the cleaned matrix, which on real data is as large as the
#' layer itself.
#'
#' The input must be the one QC ran on: the bundle remembers a checksum of
#' its expression matrix and refuses any other.
#'
#' Results saved before the record existed carry the cleaned input itself
#' (`results$cleaned_input`); it is returned as it is, and `input` is then
#' not needed.
#'
#' @param bundle An `analysis_bundle` from [run_qc()].
#' @param input The `omics_input` passed to [run_qc()].
#'
#' @return An `omics_input`: the cleaned input.
#' @export
#' @family qc
#' @examples
#' set.seed(1)
#' expr <- matrix(rnorm(60, 20), nrow = 10,
#'                dimnames = list(paste0("p", 1:10), paste0("s", 1:6)))
#' expr[1, 1:5] <- NA
#' meta <- data.frame(group = rep(c("A", "B"), each = 3),
#'                    row.names = colnames(expr))
#' feat <- data.frame(feature_id = rownames(expr), row.names = rownames(expr))
#' input <- omics_input(expr, meta, feat, omics_type = "proteomics",
#'                      assay_type = "normalized_intensity")
#' qc <- run_qc(input, impute_method = "min", outlier_method = "none")
#' cleaned <- qc_cleaned_input(qc, input)
#' dim(cleaned$expr_mat)
qc_cleaned_input <- function(bundle, input = NULL) {
  if (!is_analysis_bundle(bundle) || !identical(bundle$analysis_name, "run_qc")) {
    stop("`bundle` must be an analysis_bundle from run_qc().", call. = FALSE)
  }
  old <- bundle$results$cleaned_input
  if (!is.null(old)) return(old)
  rec <- bundle$results$cleaning
  if (is.null(rec)) {
    stop("This QC result carries neither its cleaned input nor the record ",
         "to rebuild it; run run_qc() again.", call. = FALSE)
  }
  if (is.null(input)) {
    stop("The cleaned input is rebuilt from the input QC ran on; pass it as ",
         "`input` (qc_cleaned_input(bundle, input)).", call. = FALSE)
  }
  validate_omics_input(input)
  if (!identical(rlang::hash(input$expr_mat), rec$input_hash)) {
    stop("`input` is not the data this QC result was computed on (its ",
         "values, features or samples differ). Pass the same input, or run ",
         "run_qc() again on this one.", call. = FALSE)
  }
  features <- setdiff(rownames(input$expr_mat), rec$dropped_features)
  cleaned <- subset_omics(input, samples = rec$kept_samples, features = features)
  qc_apply_cleaning(cleaned, rec, bundle$results$qc_summary$imputation)
}

# What an imputed layer is labelled. Proteomics log intensities that had
# their gaps filled are `imputed_intensity`; anything else keeps its
# label (a linear scale was imputed on log2 and transformed back).
qc_cleaned_assay_type <- function(cleaned, imputation) {
  at <- cleaned$assay_type
  if (!is.null(imputation) && identical(cleaned$omics_type, "proteomics") &&
      isTRUE(at %in% c("normalized_intensity", "filtered_intensity"))) {
    return("imputed_intensity")
  }
  at
}

# `cleaned` is the input with the kept samples and features. The record
# says what went into its missing cells.
qc_apply_cleaning <- function(cleaned, rec, imputation) {
  if (!is.null(imputation)) {
    # The matrix before imputation, unless the input already carried one
    # (after normalize_omics() it holds the linear values).
    cleaned$raw_mat <- cleaned$raw_mat %||% cleaned$expr_mat
    m <- cleaned$expr_mat
    m[imputation$imputed_cells] <- rec$imputed_values
    cleaned$expr_mat <- m
  }
  cleaned$assay_type <- rec$assay_type
  cleaned
}

# What plot_qc() draws, computed once from the cleaned input, so the plots
# need neither the cleaned input nor the layer. Each part is small: a few
# numbers per sample, and for the imputation view a summary of the
# observed values (the imputed ones are in the cleaning record).
qc_plot_data <- function(cleaned, imputation, outliers, sd_threshold) {
  log_mat <- qc_log_scale(cleaned)$mat
  out <- list(meta_df = cleaned$meta_df)

  out$pca <- tryCatch({
    mat <- mean_impute_rows(log_mat)
    if (ncol(mat) < 2L) stop("Need at least 2 samples to draw a PCA scatter.")
    pca <- pca_over_samples(mat)
    list(scores = pca$x[, 1:2, drop = FALSE],
         var_pct = (pca$sdev^2) / sum(pca$sdev^2) * 100,
         n_features = nrow(mat),
         n_dropped = attr(pca, "n_dropped") %||% 0L)
  }, error = function(e) list(error = conditionMessage(e)))

  # Connectivity is drawn from the outlier test when it ran; otherwise it
  # is worked out here, as the plot used to do on the cleaned input.
  # Quietly: a sample with no spread (every value imputed to one number)
  # makes cor() warn, and this is a panel nobody may open -- it shows the
  # missing correlation as a missing bar.
  if (is.null(qc_connectivity_from_outliers(outliers))) {
    out$connectivity <- tryCatch(
      suppressWarnings(
        qc_outliers_connectivity(log_mat, sd_threshold = sd_threshold)$stats),
      error = function(e) list(error = conditionMessage(e)))
  }

  if (!is.null(imputation)) {
    # Observed against imputed values, on the log scale. Every observed
    # value would be the matrix again; a thousand quantiles draw the same
    # curve, with the smoothing the full set would have had.
    mat <- cleaned$expr_mat
    if (!cleaned$assay_type %in% LOG_SCALE_ASSAY_TYPES) mat <- log2(mat)
    observed <- as.numeric(mat[-imputation$imputed_cells])
    observed <- observed[is.finite(observed)]
    out$imputation <- if (length(observed) >= 2L) {
      list(observed_quantiles = stats::quantile(
             observed, probs = seq(0, 1, length.out = 1001L), names = FALSE),
           observed_bw = stats::bw.nrd0(observed),
           n_observed = length(observed))
    } else {
      list(observed_quantiles = observed, observed_bw = NA_real_,
           n_observed = length(observed))
    }
  }
  out
}

qc_connectivity_from_outliers <- function(out) {
  if (!is.null(out$by_method) && "connectivity" %in% names(out$by_method)) {
    out$by_method$connectivity$stats
  } else if (identical(out$method, "connectivity")) {
    out$stats
  }
}
