#' Read the imported layer used for quality control
#'
#' After applying QC, analyses read the cleaned experiment while QC continues
#' to read the imported layer. This allows a filter to be relaxed later.
#' @param project An `omics_project`.
#' @param layer Experiment name.
#' @return An `omics_input`.
#' @export
qc_baseline_input <- function(project, layer) {
  assert_project(project)
  assert_string(layer, "layer")
  x <- project$experiments[[layer]]
  if (is.null(x)) stop("Unknown `layer`: ", layer, call. = FALSE)
  if (!is.null(x$applied_qc) && !is.null(project$qc_baselines[[layer]])) {
    base <- project$qc_baselines[[layer]]
    if (length(x$excluded_samples)) {
      base <- subset_omics_samples(base, setdiff(colnames(base$expr_mat), x$excluded_samples))
      base$excluded_samples <- x$excluded_samples
    }
    return(base)
  }
  x
}

#' Apply a QC result to downstream analyses
#'
#' Replaces the selected experiment with its cleaned input, retains its
#' imported baseline, and clears results computed before this change.
#' Merely running QC does not change downstream data.
#' @param project An `omics_project`.
#' @param layer Experiment name.
#' @param bundle A QC `analysis_bundle` computed from this layer's baseline.
#' @return The updated `omics_project`.
#' @export
apply_qc <- function(project, layer, bundle) {
  base <- qc_baseline_input(project, layer)
  if (!is_analysis_bundle(bundle) || !identical(bundle$analysis_name, "run_qc")) {
    stop("`bundle` must be a QC analysis_bundle.", call. = FALSE)
  }
  if (!is.null(bundle$input_info$layer) && !identical(bundle$input_info$layer, layer)) {
    stop("`bundle` belongs to a different layer.", call. = FALSE)
  }
  cleaned <- qc_cleaned_input(bundle, base)
  cleaned$applied_qc <- bundle$params
  cleaned$source_fingerprint <- paste0(base$source_fingerprint %||% "", ":qc=",
                                      rlang::hash(list(cleaned$expr_mat, cleaned$meta_df)))
  project$qc_baselines[[layer]] <- base
  project$experiments[[layer]] <- cleaned
  bundle$input_info$layer <- layer
  project$bundles <- list(qc = bundle)
  project
}

# The 0.1 desktop edition kept cleaned layers outside experiments. Preserve
# their meaning when opening these files with the shared 0.2 implementation.
# Stale preprocessing (a subsequently replaced import) must never be applied.
migrate_desktop_qc <- function(project) {
  for (tag in intersect(names(project$preprocessing), names(project$experiments))) {
    pp <- project$preprocessing[[tag]]
    base <- project$experiments[[tag]]
    if (is.null(base$applied_qc) && !is.null(base$version) &&
        identical(pp$base_input_version, base$version) && is_omics_input(pp$cleaned)) {
      project$qc_baselines[[tag]] <- base
      cleaned <- pp$cleaned
      cleaned$applied_qc <- pp$params
      cleaned$source_fingerprint <- paste0(base$source_fingerprint %||% "", ":qc=",
                                          rlang::hash(list(cleaned$expr_mat, cleaned$meta_df)))
      project$experiments[[tag]] <- cleaned
    }
  }
  project
}
