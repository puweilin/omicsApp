#' Detect sample-level outliers
#'
#' Outlier detection on an [omics_input()]. Four methods are supported:
#'
#' * `"pca"` — flags samples whose absolute z-score on PC1 or PC2 exceeds
#'   `sd_threshold`.
#' * `"connectivity"` — flags samples whose mean inter-sample correlation
#'   is more than `sd_threshold` standard deviations below the cohort mean
#'   (i.e. poorly-connected samples).
#' * `"iqr"` — flags samples whose mean log-intensity falls outside
#'   `[Q1 - k*IQR, Q3 + k*IQR]`; `sd_threshold` is reused as `k`.
#' * `"loo"` — a leave-one-out test built for small studies. Each
#'   sample's distance to its most similar other sample (correlation
#'   distance, `sqrt(1 - r)`, on a log scale) is compared with the same
#'   distance computed among the remaining samples without it: their
#'   median, and their MAD as the spread. A sample is flagged when it sits
#'   more than `sd_threshold` spreads beyond that median. The spread is
#'   never taken as less than `0.17` (on the natural-log scale, a factor
#'   of about 1.2): with four samples the MAD of the other three can be
#'   near zero by chance, and a z-score over nothing is not evidence. At
#'   the default threshold of 3 a flagged sample is therefore at least
#'   1.67 times further from its nearest neighbour than the others are
#'   from theirs. The nearest neighbour, rather than the mean over all
#'   samples, is what keeps group structure from looking like an
#'   outlier: in a 3 vs 3 design every sample has a close partner in its
#'   own group, while a failed run or a swapped-in sample has none. It
#'   needs at least four samples.
#'
#' All four methods are pure base R; no Suggests packages required.
#'
#' They run on a log scale whatever the layer holds: raw counts become
#' log2-CPM and linear intensities, TPM or FPKM become log2(x + 1), the
#' same conversion [run_diff()] makes for limma. On the raw scale a
#' correlation or a PCA is driven by the handful of most abundant
#' features, and the deepest library looks like the odd one out.
#'
#' A z-score over `n` samples cannot exceed `(n - 1) / sqrt(n)`, so with
#' ten or fewer samples a threshold of 3 can never be reached by `"pca"`
#' or `"connectivity"`. The result then carries a `note` saying so; the
#' `"loo"` test does not have this limit, since the sample under test is
#' not part of the reference it is scored against.
#'
#' @param input An `omics_input`.
#' @param method One or more of `"pca"`, `"connectivity"`, `"iqr"`,
#'   `"loo"`. Multiple methods run in parallel and the union of the
#'   flagged sets is taken; the default runs all four.
#' @param sd_threshold Z-score or IQR multiplier used to flag outliers
#'   (default `3`); for `"loo"`, the number of robust spreads.
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
  method = c("pca", "connectivity", "iqr", "loo"),
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
      iqr          = qc_outliers_iqr(expr_mat, sd_threshold),
      loo          = qc_outliers_loo(expr_mat, sd_threshold)
    )
  })
  names(per_method) <- method

  # The z-score methods cannot flag anything at this sample size.
  n <- ncol(expr_mat)
  unreachable <- any(method %in% c("pca", "connectivity")) && n >= 2L &&
    (n - 1) / sqrt(n) <= sd_threshold
  note <- c(scaled$note, if (unreachable) {
    if ("loo" %in% method && n >= LOO_MIN_SAMPLES) sprintf(
      "With %d samples a z-score cannot exceed %.2f, so the PCA and connectivity tests flag nothing at a threshold of %g; the leave-one-out test still applies.",
      n, (n - 1) / sqrt(n), sd_threshold)
    else sprintf(
      "With %d samples a z-score cannot exceed %.2f, so a threshold of %g flags nothing; inspect the PCA plot instead.",
      n, (n - 1) / sqrt(n), sd_threshold)
  }, if ("loo" %in% method && n < LOO_MIN_SAMPLES) sprintf(
    "The leave-one-out test needs at least %d samples; with %d it flags nothing.",
    LOO_MIN_SAMPLES, n))

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
  cor_mat <- sample_cor(expr_mat)
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

# Leave-one-out against the nearest neighbour; see the "loo" entry in
# the documentation above for why each choice was made. The floor and
# the minimum sample count were set by simulation (clean log2 intensity
# matrices of 4-10 samples in one to three groups, with left-censored
# missingness and unequal per-sample noise): at threshold 3 a clean
# study flags a sample in at most ~2% of seeds, while a sample whose
# intensities were scrambled by +/-1.5 log2 is flagged in nearly all.
LOO_MIN_SPREAD <- 0.17
LOO_MIN_SAMPLES <- 4L

qc_outliers_loo <- function(expr_mat, sd_threshold) {
  n <- ncol(expr_mat)
  ids <- colnames(expr_mat)
  empty <- data.frame(sample_id = ids, nearest_sample = NA_character_,
                      nearest_correlation = NA_real_,
                      reference_correlation = NA_real_, z_score = NA_real_,
                      is_outlier = rep(FALSE, n),
                      stringsAsFactors = FALSE, row.names = NULL)
  if (n < LOO_MIN_SAMPLES) {
    return(list(method = "loo", stats = empty, flagged_samples = character(0)))
  }
  r <- sample_cor(expr_mat)
  # Log of the correlation distance sqrt(1 - r). On the log scale a
  # sample twice as noisy as the rest is the same step away whatever the
  # overall noise level, and the spread of the reference is symmetric
  # rather than a long right tail that a MAD reads as outliers.
  d <- 0.5 * log(pmax(1 - r, 1e-12))
  # A pair with too few values in common to correlate is not a near
  # neighbour of anything.
  d[is.na(d)] <- Inf
  diag(d) <- Inf
  # The nearest and second-nearest neighbour of every sample: leaving
  # sample i out changes another sample's nearest neighbour only when
  # that neighbour was i, and then the second-nearest takes its place.
  ord <- apply(d, 1L, order)
  idx <- seq_len(n)
  first <- ord[1L, ]
  nn1 <- d[cbind(idx, first)]
  nn2 <- d[cbind(idx, ord[2L, ])]

  z <- ref <- rep(NA_real_, n)
  for (i in idx) {
    js <- idx[-i]
    others <- ifelse(first[js] == i, nn2[js], nn1[js])
    others <- others[is.finite(others)]
    if (length(others) < LOO_MIN_SAMPLES - 1L || !is.finite(nn1[i])) next
    ref[i] <- stats::median(others)
    spread <- max(stats::mad(others), LOO_MIN_SPREAD)
    z[i] <- (nn1[i] - ref[i]) / spread
  }
  is_out <- !is.na(z) & z > sd_threshold

  # Reported as correlations: "0.91 with its closest sample, where the
  # others manage 0.98" is the statement a reader can check.
  to_cor <- function(x) 1 - exp(2 * x)
  stats_df <- data.frame(
    sample_id = ids,
    nearest_sample = ifelse(is.finite(nn1), ids[first], NA_character_),
    nearest_correlation = ifelse(is.finite(nn1), to_cor(nn1), NA_real_),
    reference_correlation = to_cor(ref),
    z_score = z,
    is_outlier = is_out,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  list(method = "loo", stats = stats_df,
       flagged_samples = stats_df$sample_id[is_out])
}

# ---- internal shared helpers -------------------------------------------

# Sample-by-sample correlation, shared by the connectivity and
# leave-one-out tests: on 60k x 300 it is six seconds, and "all" methods
# would otherwise pay it twice per QC run. Remembered by the matrix's
# content (the last one), as pca_over_samples() remembers its PCA.
sample_cor <- function(expr_mat) {
  key <- rlang::hash(expr_mat)
  hit <- .cor_cache[[key]]
  if (!is.null(hit)) return(hit)
  r <- if (anyNA(expr_mat)) pairwise_cor(expr_mat) else stats::cor(expr_mat)
  rm(list = ls(.cor_cache), envir = .cor_cache)
  assign(key, r, envir = .cor_cache)
  r
}
.cor_cache <- new.env(parent = emptyenv())

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
