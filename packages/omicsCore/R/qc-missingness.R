#' Per-sample / per-feature missingness summary
#'
#' Computes per-sample and per-feature missing rates for an
#' [omics_input()] and flags samples / features whose missing rate exceeds
#' the supplied thresholds.
#'
#' @param input An `omics_input`.
#' @param sample_missing_cutoff Optional sample missing-rate threshold in
#'   `[0, 1]`. `NULL` (default) leaves all samples unflagged.
#' @param feature_missing_cutoff Feature missing-rate threshold in `[0, 1]`.
#'
#' @return A list with:
#'   \describe{
#'     \item{`sample_metrics`}{`data.frame` with `sample_id`, `missing_rate`.}
#'     \item{`feature_metrics`}{`data.frame` with `feature_id`, `missing_rate`.}
#'     \item{`flagged_samples`}{Character vector of sample IDs exceeding
#'       `sample_missing_cutoff`.}
#'     \item{`flagged_features`}{Character vector of feature IDs exceeding
#'       `feature_missing_cutoff`.}
#'     \item{`settings`}{Echo of the thresholds used.}
#'   }
#' @export
#' @family qc
qc_missingness <- function(
  input,
  sample_missing_cutoff = NULL,
  feature_missing_cutoff = 0.5
) {
  assert_number(sample_missing_cutoff, "sample_missing_cutoff",
                lower = 0, upper = 1, allow_null = TRUE)
  assert_number(feature_missing_cutoff, "feature_missing_cutoff", lower = 0, upper = 1)
  validate_omics_input(input)
  expr_mat <- input$expr_mat

  sample_metrics <- data.frame(
    sample_id = colnames(expr_mat),
    missing_rate = colMeans(is.na(expr_mat)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  feature_metrics <- data.frame(
    feature_id = rownames(expr_mat),
    missing_rate = rowMeans(is.na(expr_mat)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  flagged_samples <- if (is.null(sample_missing_cutoff)) {
    character(0)
  } else {
    sample_metrics$sample_id[sample_metrics$missing_rate > sample_missing_cutoff]
  }

  flagged_features <- feature_metrics$feature_id[
    feature_metrics$missing_rate > feature_missing_cutoff
  ]

  list(
    sample_metrics = sample_metrics,
    feature_metrics = feature_metrics,
    flagged_samples = flagged_samples,
    flagged_features = flagged_features,
    settings = list(
      sample_missing_cutoff = sample_missing_cutoff,
      feature_missing_cutoff = feature_missing_cutoff
    )
  )
}
