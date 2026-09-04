# Shared preparation for every PCA in the package.
#
# `prcomp(scale. = TRUE)` divides each column by its standard deviation,
# so one feature that never varies stops the whole decomposition with
#
#   cannot rescale a constant/zero column to unit variance
#
# In RNA-seq that is the normal case rather than a corrupt file: a counts
# matrix carries every annotated gene, and thousands of them are zero in
# every sample of a given tissue. The error names neither the features
# nor the reason, so it reads as "the data is broken".
#
# Three call sites need this -- outlier detection, the QC scatter and the
# Differential scatter -- and they had drifted: one filtered, two did
# not. One implementation is what keeps that from happening again.

#' Row variances, without a per-row function call
#'
#' `apply(mat, 1, var)` calls `var()` once per feature, which on a 63k
#' gene matrix is tens of thousands of calls to answer a question that is
#' one pass over the data.
#'
#' The denominator is `ncol - 1` regardless of how many values were NA,
#' which is wrong as a variance and does not matter here: every caller
#' asks only whether the result is greater than zero, and no amount of
#' denominator changes that sign.
#'
#' @param mat Numeric matrix, features in rows.
#' @return Numeric vector, one per row. Rows that are entirely NA get 0.
#' @keywords internal
#' @noRd
row_variance <- function(mat) {
  n <- ncol(mat)
  if (n < 2L) return(rep(0, nrow(mat)))
  mu <- rowMeans(mat, na.rm = TRUE)
  mu[!is.finite(mu)] <- 0
  centred <- mat - mu
  rowSums(centred * centred, na.rm = TRUE) / (n - 1L)
}

#' Drop features a scaled PCA cannot use
#'
#' @param mat Numeric matrix, features in rows, samples in columns.
#' @param min_features Fewest features that still make a PCA meaningful.
#' @return The matrix without its constant features, carrying an
#'   `n_dropped` attribute so a caller can say what happened.
#' @keywords internal
#' @noRd
drop_constant_features <- function(mat, min_features = 2L) {
  v <- row_variance(mat)
  keep <- is.finite(v) & v > 0
  out <- mat[keep, , drop = FALSE]
  if (nrow(out) < min_features) {
    stop("Need at least ", min_features,
         " features that vary across samples for PCA; ",
         nrow(out), " of ", nrow(mat), " do.", call. = FALSE)
  }
  attr(out, "n_dropped") <- sum(!keep)
  out
}

#' PCA over samples, with the constant features removed first
#'
#' @param mat Numeric matrix, features in rows, samples in columns. NAs
#'   should already be imputed; `prcomp()` cannot take them.
#' @return The `prcomp` object, carrying `n_dropped`.
#' @keywords internal
#' @noRd
pca_over_samples <- function(mat) {
  if (ncol(mat) < 2L) {
    stop("Need at least 2 samples for PCA.", call. = FALSE)
  }
  kept <- drop_constant_features(mat)
  pca <- stats::prcomp(t(kept), scale. = TRUE)
  attr(pca, "n_dropped") <- attr(kept, "n_dropped")
  pca
}
