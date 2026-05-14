#' Impute missing values in an expression matrix
#'
#' Dispatcher over several imputation backends:
#'
#' * `"none"` — return the matrix unchanged.
#' * `"min"` — replace `NA` per feature with that feature's minimum observed
#'   value. Pure base R, always available.
#' * `"half_min"` — replace `NA` per feature with half its minimum observed
#'   value. Pure base R, common proteomics default.
#' * `"mean"` — replace `NA` per feature with its row mean. Pure base R.
#' * `"knn"` — k-nearest-neighbor imputation via the `impute` Bioconductor
#'   package (Suggests).
#' * `"missforest"` — Random-forest imputation via `missForest` (Suggests).
#' * `"bpca"` — Bayesian PCA via `pcaMethods` (Suggests).
#'
#' Backends listed in Suggests raise an instructive error pointing at
#' [install_optional()] when the package is missing.
#'
#' @param mat A numeric matrix.
#' @param method One of the methods listed above.
#' @param ... Additional arguments forwarded to the backend.
#'
#' @return A numeric matrix with the same dimensions and names as `mat`.
#' @export
#' @family qc
impute_matrix <- function(
  mat,
  method = c("none", "min", "half_min", "mean", "knn", "missforest", "bpca"),
  ...
) {
  method <- match.arg(method)
  mat <- as.matrix(mat)
  if (method == "none" || !anyNA(mat)) {
    return(mat)
  }

  switch(method,
    min       = impute_row_min(mat, scale = 1),
    half_min  = impute_row_min(mat, scale = 0.5),
    mean      = mean_impute_rows(mat),
    knn       = impute_knn(mat, ...),
    missforest= impute_missforest(mat, ...),
    bpca      = impute_bpca(mat, ...)
  )
}

# ---- backends ----------------------------------------------------------

impute_row_min <- function(mat, scale = 1) {
  row_min <- apply(mat, 1L, min, na.rm = TRUE)
  row_min[!is.finite(row_min)] <- 0
  idx <- which(is.na(mat), arr.ind = TRUE)
  mat[idx] <- scale * row_min[idx[, 1L]]
  mat
}

impute_knn <- function(mat, k = 10, ...) {
  if (!requireNamespace("impute", quietly = TRUE)) {
    stop(
      "Package 'impute' is required for method = 'knn'. ",
      "Install it with: omicsCore::install_optional('imputation') ",
      "or BiocManager::install('impute').",
      call. = FALSE
    )
  }
  res <- impute::impute.knn(mat, k = k, ...)
  out <- res$data
  dimnames(out) <- dimnames(mat)
  out
}

impute_missforest <- function(mat, maxiter = 10, ntree = 100, ...) {
  if (!requireNamespace("missForest", quietly = TRUE)) {
    stop(
      "Package 'missForest' is required for method = 'missforest'. ",
      "Install it with: omicsCore::install_optional('imputation').",
      call. = FALSE
    )
  }
  # missForest imputes columns; transpose so each row is one observation.
  df <- as.data.frame(t(mat))
  res <- missForest::missForest(df, maxiter = maxiter, ntree = ntree, ...)
  out <- t(as.matrix(res$ximp))
  dimnames(out) <- dimnames(mat)
  out
}

impute_bpca <- function(mat, n_pcs = 3, ...) {
  if (!requireNamespace("pcaMethods", quietly = TRUE)) {
    stop(
      "Package 'pcaMethods' is required for method = 'bpca'. ",
      "Install it with: omicsCore::install_optional('imputation').",
      call. = FALSE
    )
  }
  pc <- pcaMethods::pca(t(mat), method = "bpca", nPcs = n_pcs, ...)
  out <- t(pcaMethods::completeObs(pc))
  dimnames(out) <- dimnames(mat)
  out
}
