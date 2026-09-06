# Per-sample sequencing depth and detection.
#
# The missingness panels answer the question a proteomics run raises: a
# peptide that was not detected is a hole in the matrix, and how many
# holes each sample has is the first thing to look at.
#
# A counts matrix has no holes. Every gene has a number for every sample
# and most of those numbers are zero, so the missingness panel reports
# "63,241 features, all at 0%" -- a true statement with no information
# in it, and the panel it fills is the one that should have been showing
# whether a library was under-sequenced.
#
# What replaces it, per sample:
#
#   library_size    total counts. A library an order of magnitude below
#                   its neighbours is under-sequenced, and every
#                   downstream comparison involving it is weaker for it.
#   n_detected      genes with any signal. Falls when a library is
#                   shallow, and also when a sample is degraded -- two
#                   different problems, and the pair of numbers tells
#                   them apart: shallow drops both, degraded drops
#                   detection while the total holds.
#
# Computed for proteomics too, where they read as total intensity and
# features quantified. Both are one pass over the matrix, so there is
# nothing to gain by deciding in advance which one gets looked at.

#' Per-sample depth and detection
#'
#' @param input An [omics_input].
#' @return A data frame with one row per sample: `sample_id`,
#'   `library_size`, `n_detected`, `detection_rate`, and
#'   `library_size_ratio` (each library over the median, so "half the
#'   depth of a typical sample" is readable without arithmetic).
#' @export
#' @family qc
qc_depth <- function(input) {
  validate_omics_input(input)
  mat <- input$expr_mat
  n_feat <- nrow(mat)

  lib <- colSums(mat, na.rm = TRUE)
  # Zero and NA both mean "nothing seen here": a counts matrix says it
  # with 0 and an intensity matrix with NA, and the question -- how much
  # of the assay did this sample return -- is the same one.
  detected <- colSums(!is.na(mat) & mat > 0)

  med <- stats::median(lib[is.finite(lib)])
  data.frame(
    sample_id          = colnames(mat),
    library_size       = as.numeric(lib),
    n_detected         = as.integer(detected),
    detection_rate     = as.numeric(detected) / max(n_feat, 1L),
    library_size_ratio = if (is.finite(med) && med > 0) {
      as.numeric(lib) / med
    } else {
      rep(NA_real_, ncol(mat))
    },
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

#' Samples whose depth is far from the rest
#'
#' Flagged on the ratio to the median rather than an absolute count: what
#' counts as a shallow library depends entirely on the experiment, and a
#' fixed threshold would be wrong for every study but one.
#'
#' @param depth_df Output of [qc_depth()].
#' @param min_ratio Flag samples below this fraction of the median.
#' @return Character vector of sample ids.
#' @export
#' @family qc
qc_depth_outliers <- function(depth_df, min_ratio = 0.3) {
  assert_data_frame(depth_df, "depth_df")
  check_required_cols(depth_df, c("sample_id", "library_size_ratio"),
                      object_name = "depth_df")
  assert_number(min_ratio, "min_ratio", lower = 0)
  r <- depth_df$library_size_ratio
  depth_df$sample_id[!is.na(r) & r < min_ratio]
}
