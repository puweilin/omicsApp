#' Subset an `omics_input`
#'
#' Subsets an `omics_input` by sample IDs and/or feature IDs. Either argument
#' may be `NULL` to retain all entries on that axis. The matching matrices
#' (`raw_mat`, `normalized_mat`) and metadata are subset accordingly.
#'
#' @param input An `omics_input`.
#' @param samples Character vector of sample IDs to retain, or `NULL`.
#' @param features Character vector of feature IDs to retain, or `NULL`.
#'
#' @return A new `omics_input`.
#' @export
#' @family omics_input
subset_omics <- function(input, samples = NULL, features = NULL) {
  assert_character(samples, "samples", allow_null = TRUE)
  assert_character(features, "features", allow_null = TRUE)
  validate_omics_input(input)
  out <- input
  if (!is.null(samples)) {
    out <- subset_omics_samples(out, samples)
  }
  if (!is.null(features)) {
    out <- subset_omics_features(out, features)
  }
  out
}

#' Subset an `omics_input` by sample IDs
#'
#' @param omics_input A validated `omics_input`.
#' @param sample_ids Sample identifiers to retain.
#'
#' @return A new `omics_input` containing the selected samples.
#' @export
#' @family omics_input
subset_omics_samples <- function(omics_input, sample_ids) {
  assert_character(sample_ids, "sample_ids")
  validate_omics_input(omics_input)

  sample_ids <- intersect(sample_ids, colnames(omics_input$expr_mat))
  if (length(sample_ids) == 0) {
    stop("No requested sample IDs were found in `omics_input`.")
  }

  expr_mat <- omics_input$expr_mat[, sample_ids, drop = FALSE]
  raw_mat <- if (!is.null(omics_input$raw_mat)) {
    omics_input$raw_mat[, sample_ids, drop = FALSE]
  } else NULL
  normalized_mat <- if (!is.null(omics_input$normalized_mat)) {
    omics_input$normalized_mat[, sample_ids, drop = FALSE]
  } else NULL
  meta_df <- omics_input$meta_df[sample_ids, , drop = FALSE]

  new_omics_input(
    omics_type = omics_input$omics_type,
    assay_type = omics_input$assay_type,
    expr_mat = expr_mat,
    raw_mat = raw_mat,
    normalized_mat = normalized_mat,
    meta_df = meta_df,
    feature_df = omics_input$feature_df,
    raw_object = omics_input$raw_object
  )
}

#' Subset an `omics_input` by feature IDs
#'
#' @param omics_input A validated `omics_input`.
#' @param feature_ids Feature identifiers to retain.
#'
#' @return A new `omics_input` containing the selected features.
#' @export
#' @family omics_input
subset_omics_features <- function(omics_input, feature_ids) {
  assert_character(feature_ids, "feature_ids")
  validate_omics_input(omics_input)

  feature_ids <- intersect(feature_ids, rownames(omics_input$expr_mat))
  if (length(feature_ids) == 0) {
    stop("No requested feature IDs were found in `omics_input`.")
  }

  expr_mat <- omics_input$expr_mat[feature_ids, , drop = FALSE]
  raw_mat <- if (!is.null(omics_input$raw_mat)) {
    omics_input$raw_mat[feature_ids, , drop = FALSE]
  } else NULL
  normalized_mat <- if (!is.null(omics_input$normalized_mat)) {
    omics_input$normalized_mat[feature_ids, , drop = FALSE]
  } else NULL
  feature_df <- omics_input$feature_df[
    match(feature_ids, omics_input$feature_df$feature_id), , drop = FALSE
  ]
  rownames(feature_df) <- feature_df$feature_id

  new_omics_input(
    omics_type = omics_input$omics_type,
    assay_type = omics_input$assay_type,
    expr_mat = expr_mat,
    raw_mat = raw_mat,
    normalized_mat = normalized_mat,
    meta_df = omics_input$meta_df,
    feature_df = feature_df,
    raw_object = omics_input$raw_object
  )
}

#' Drop samples with missing metadata
#'
#' Retains only samples whose `meta_df` rows are complete on the requested
#' columns.
#'
#' @param omics_input A validated `omics_input`.
#' @param cols Character vector of metadata columns that must be non-`NA`.
#'
#' @return A new `omics_input` with incomplete samples removed.
#' @export
#' @family omics_input
drop_meta_na <- function(omics_input, cols) {
  assert_names(cols, "cols")
  validate_omics_input(omics_input)
  check_required_cols(omics_input$meta_df, cols, object_name = "meta_df")

  keep_samples <- rownames(omics_input$meta_df)[
    stats::complete.cases(omics_input$meta_df[, cols, drop = FALSE])
  ]
  subset_omics_samples(omics_input, keep_samples)
}

#' Select features passing a missingness threshold
#'
#' @param omics_input A validated `omics_input`.
#' @param feature_missing_cutoff Maximum allowed per-feature missing rate
#'   (0–1).
#'
#' @return A new `omics_input` filtered by feature missingness.
#' @export
#' @family omics_input
select_complete_cases <- function(omics_input, feature_missing_cutoff = 1) {
  assert_number(feature_missing_cutoff, "feature_missing_cutoff", lower = 0, upper = 1)
  validate_omics_input(omics_input)

  feature_missing <- rowMeans(is.na(omics_input$expr_mat))
  keep_features <- names(feature_missing)[feature_missing <= feature_missing_cutoff]
  subset_omics_features(omics_input, keep_features)
}
