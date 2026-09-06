#' Detect sample-level outliers
#'
#' Outlier detection on an [omics_input()]. Three methods are supported:
#'
#' * `"pca"` — flags samples whose absolute z-score on PC1 or PC2 exceeds
#'   `sd_threshold`.
#' * `"connectivity"` — flags samples whose mean inter-sample correlation
#'   is more than `sd_threshold` standard deviations below the cohort mean
#'   (i.e. poorly-connected samples).
#' * `"iqr"` — flags samples whose mean log-intensity falls outside
#'   `[Q1 - k*IQR, Q3 + k*IQR]`; `sd_threshold` is reused as `k`.
#'
#' All three methods are pure base R; no Suggests packages required.
#'
#' @param input An `omics_input`.
#' @param method One of `"pca"`, `"connectivity"`, `"iqr"`. Multiple methods
#'   may be supplied to run them in parallel and take the union of the
#'   flagged sets.
#' @param sd_threshold Z-score or IQR multiplier used to flag outliers
#'   (default `3`).
#'
#' @return When `method` is length 1, a list with `method`, `stats`,
#'   `flagged_samples`. When length > 1, the same shape with `stats` row-bound
#'   across methods and a `by_method` field holding per-method results.
#' @export
#' @family qc
qc_outliers <- function(
  input,
  method = c("pca", "connectivity", "iqr"),
  sd_threshold = 3
) {
  assert_number(sd_threshold, "sd_threshold", lower = 0)
  validate_omics_input(input)
  method <- match.arg(method, several.ok = TRUE)
  expr_mat <- input$expr_mat

  per_method <- lapply(method, function(m) {
    switch(m,
      pca          = qc_outliers_pca(expr_mat, sd_threshold),
      connectivity = qc_outliers_connectivity(expr_mat, sd_threshold),
      iqr          = qc_outliers_iqr(expr_mat, sd_threshold)
    )
  })
  names(per_method) <- method

  if (length(per_method) == 1L) {
    return(per_method[[1L]])
  }

  all_stats <- do.call(rbind, lapply(per_method, function(r) {
    df <- r$stats[, c("sample_id", "is_outlier")]
    df$method <- r$method
    df
  }))
  rownames(all_stats) <- NULL

  list(
    method = method,
    stats = all_stats,
    flagged_samples = unique(unlist(lapply(per_method, `[[`, "flagged_samples"))),
    by_method = per_method
  )
}

# ---- internal per-method implementations -------------------------------

qc_outliers_pca <- function(expr_mat, sd_threshold) {
  # prcomp() cannot handle NA values; mean-impute per feature first.
  # pca_over_samples() then drops the features that never vary, which a
  # scaled PCA cannot use -- see pca-utils.R.
  mat_for_pca <- mean_impute_rows(expr_mat)
  pca_res <- pca_over_samples(mat_for_pca)
  coords <- as.data.frame(pca_res$x[, 1:2, drop = FALSE])
  coords$sample_id <- rownames(coords)
  rownames(coords) <- NULL

  coords$z_pc1 <- safe_z(coords$PC1)
  coords$z_pc2 <- safe_z(coords$PC2)
  coords$is_outlier <- abs(coords$z_pc1) > sd_threshold |
    abs(coords$z_pc2) > sd_threshold

  list(
    method = "pca",
    stats = coords[, c("sample_id", "PC1", "PC2", "z_pc1", "z_pc2", "is_outlier")],
    flagged_samples = coords$sample_id[coords$is_outlier]
  )
}

qc_outliers_connectivity <- function(expr_mat, sd_threshold) {
  cor_mat <- stats::cor(expr_mat, use = "pairwise.complete.obs", method = "pearson")
  diag(cor_mat) <- NA_real_
  mean_cor <- colMeans(cor_mat, na.rm = TRUE)
  z_score <- safe_z(mean_cor)

  stats_df <- data.frame(
    sample_id = names(mean_cor),
    mean_correlation = mean_cor,
    z_score = z_score,
    is_outlier = z_score < -sd_threshold,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  list(
    method = "connectivity",
    stats = stats_df,
    flagged_samples = stats_df$sample_id[stats_df$is_outlier]
  )
}

qc_outliers_iqr <- function(expr_mat, k) {
  sample_means <- colMeans(expr_mat, na.rm = TRUE)
  q <- stats::quantile(sample_means, probs = c(0.25, 0.75), na.rm = TRUE)
  iqr <- q[[2L]] - q[[1L]]
  lower <- q[[1L]] - k * iqr
  upper <- q[[2L]] + k * iqr
  is_out <- !is.na(sample_means) & (sample_means < lower | sample_means > upper)

  stats_df <- data.frame(
    sample_id = colnames(expr_mat),
    mean_signal = sample_means,
    lower_fence = lower,
    upper_fence = upper,
    is_outlier = is_out,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  list(
    method = "iqr",
    stats = stats_df,
    flagged_samples = stats_df$sample_id[stats_df$is_outlier]
  )
}

# ---- internal shared helpers -------------------------------------------

safe_z <- function(x) {
  sd_x <- stats::sd(x, na.rm = TRUE)
  if (is.na(sd_x) || sd_x == 0) {
    return(rep(0, length(x)))
  }
  (x - mean(x, na.rm = TRUE)) / sd_x
}

mean_impute_rows <- function(mat) {
  if (!anyNA(mat)) return(mat)
  row_means <- rowMeans(mat, na.rm = TRUE)
  row_means[is.nan(row_means)] <- 0
  idx <- which(is.na(mat), arr.ind = TRUE)
  mat[idx] <- row_means[idx[, 1L]]
  mat
}
