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
#' They run on a log scale whatever the layer holds: raw counts become
#' log2-CPM and linear intensities, TPM or FPKM become log2(x + 1), the
#' same conversion [run_diff()] makes for limma. On the raw scale a
#' correlation or a PCA is driven by the handful of most abundant
#' features, and the deepest library looks like the odd one out.
#'
#' A z-score over `n` samples cannot exceed `(n - 1) / sqrt(n)`, so with
#' ten or fewer samples a threshold of 3 can never be reached. The result
#' then carries a `note` saying so; look at the PCA plot instead.
#'
#' @param input An `omics_input`.
#' @param method One of `"pca"`, `"connectivity"`, `"iqr"`. Multiple methods
#'   may be supplied to run them in parallel and take the union of the
#'   flagged sets.
#' @param sd_threshold Z-score or IQR multiplier used to flag outliers
#'   (default `3`).
#'
#' @return When `method` is length 1, a list with `method`, `stats`,
#'   `flagged_samples`, and, when the data were transformed or the
#'   threshold cannot be reached, `note`. When length > 1, the same shape
#'   with `stats` row-bound across methods and a `by_method` field holding
#'   per-method results.
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
  scaled <- qc_log_scale(input)
  expr_mat <- scaled$mat

  per_method <- lapply(method, function(m) {
    switch(m,
      pca          = qc_outliers_pca(expr_mat, sd_threshold),
      connectivity = qc_outliers_connectivity(expr_mat, sd_threshold),
      iqr          = qc_outliers_iqr(expr_mat, sd_threshold)
    )
  })
  names(per_method) <- method

  # The z-score methods cannot flag anything at this sample size.
  n <- ncol(expr_mat)
  unreachable <- any(method %in% c("pca", "connectivity")) && n >= 2L &&
    (n - 1) / sqrt(n) <= sd_threshold
  note <- c(scaled$note, if (unreachable) sprintf(
    "With %d samples a z-score cannot exceed %.2f, so a threshold of %g flags nothing; inspect the PCA plot instead.",
    n, (n - 1) / sqrt(n), sd_threshold))

  if (length(per_method) == 1L) {
    out <- per_method[[1L]]
    if (length(note)) out$note <- note
    return(out)
  }

  all_stats <- do.call(rbind, lapply(per_method, function(r) {
    df <- r$stats[, c("sample_id", "is_outlier")]
    df$method <- r$method
    df
  }))
  rownames(all_stats) <- NULL

  out <- list(
    method = method,
    stats = all_stats,
    flagged_samples = unique(unlist(lapply(per_method, `[[`, "flagged_samples"))),
    by_method = per_method
  )
  if (length(note)) out$note <- note
  out
}

# The matrix the outlier methods and the QC PCA look at: on a log scale,
# by the same rule run_diff() uses for limma.
qc_log_scale <- function(input) {
  sc <- prepare_diff_scale(input, "limma")
  list(mat = sc$input$expr_mat,
       note = if (!is.null(sc$note)) sprintf(
         "`%s` values were put on a log2 scale%s for outlier detection.",
         input$assay_type,
         if (identical(input$assay_type, "raw_count")) " (log2-CPM)" else ""))
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
  # Pairwise-complete correlation is one pass per pair of samples and
  # took most of a minute on 60k x 300; without missing values the plain
  # matrix product gives the same answer.
  cor_mat <- if (anyNA(expr_mat)) pairwise_cor(expr_mat) else stats::cor(expr_mat)
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
