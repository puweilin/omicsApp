#' Winsorize per-gene outliers in a count matrix
#'
#' For each expressed gene (row), counts above `Q3 + k * IQR` are clipped
#' to that bound. This removes extreme leverage points that can drive
#' spurious differential expression results while preserving the overall
#' distribution.
#'
#' The quartiles are taken over the gene's non-zero counts, and only genes
#' with at least `min_expressed` non-zero counts are winsorized. A gene that
#' is switched on in some samples and off in the rest is the case this
#' protects: over all its counts its quartiles are both 0, the bound is 0,
#' and every sample where it is expressed used to be clipped to zero --
#' turning the strongest kind of difference into no difference at all.
#' Zeros are never changed.
#'
#' A gene whose non-zero counts do not spread at all (an IQR of 0) is left
#' alone too: there is no scale to measure an outlier against, and clipping
#' at its upper quartile would flatten every count above it.
#'
#' Clipped values are rounded down to a whole number, so the matrix stays
#' a count matrix that DESeq2 and edgeR accept; an integer matrix stays
#' integer.
#'
#' `legacy = TRUE` reproduces the legacy pipeline exactly: quartiles over
#' all counts (zeros included), every gene eligible, and clipping to the
#' unrounded bound. It is kept so results produced with that pipeline can
#' be reproduced; it has the on/off problem described above.
#'
#' Pure base R; no `matrixStats` dependency.
#'
#' @param count_mat Integer or numeric count matrix (genes × samples).
#' @param k IQR multiplier for the upper fence. Larger values are more
#'   conservative (fewer values clipped). Default `20`.
#' @param min_expressed Minimum number of samples with a non-zero count for
#'   a gene to be winsorized. Default `3`. Ignored when `legacy = TRUE`.
#' @param legacy Use the legacy pipeline's rule (see Details). Default
#'   `FALSE`.
#'
#' @return A list with components:
#'   \item{`count_mat`}{Winsorized count matrix (same dimensions / names).}
#'   \item{`stats`}{`data.frame` with per-gene outlier statistics:
#'     `feature_id`, `q1`, `q3`, `iqr`, `threshold` (`NA` for a gene that
#'     was not winsorized), `max_original`, `n_clipped`, `n_nonzero` and
#'     `winsorized`.}
#'   \item{`n_clipped`}{Total number of values clipped.}
#'   \item{`n_genes_affected`}{Number of genes with at least one clipped value.}
#'   \item{`k`}{The multiplier used.}
#'   \item{`min_expressed`}{The expression floor used (`NA` for legacy).}
#'   \item{`legacy`}{Whether the legacy rule was used.}
#' @export
#' @family qc
winsorize_counts <- function(count_mat, k = 20, min_expressed = 3L,
                             legacy = FALSE) {
  count_mat <- assert_numeric_matrix(count_mat, "count_mat")
  assert_number(k, "k", lower = 0)
  assert_count(min_expressed, "min_expressed", lower = 1L)
  assert_flag(legacy, "legacy")
  n_genes <- nrow(count_mat)

  observed <- !is.na(count_mat)
  n_nonzero <- as.integer(rowSums(observed & count_mat > 0))

  # Quartiles per gene: over all counts for the legacy rule, over the
  # non-zero ones otherwise.
  quartiles <- function(x) {
    if (!legacy) x <- x[!is.na(x) & x > 0]
    if (!length(x)) return(c(NA_real_, NA_real_))
    stats::quantile(x, probs = c(0.25, 0.75), na.rm = TRUE, names = FALSE)
  }
  qmat <- matrix(vapply(seq_len(n_genes), function(i) quartiles(count_mat[i, ]),
                        numeric(2L)),
                 nrow = 2L)
  q1_vec <- qmat[1L, ]
  q3_vec <- qmat[2L, ]
  iqr_vec <- q3_vec - q1_vec
  threshold_vec <- q3_vec + k * iqr_vec

  eligible <- if (legacy) {
    rep(TRUE, n_genes)
  } else {
    n_nonzero >= min_expressed & !is.na(iqr_vec) & iqr_vec > 0
  }
  if (!legacy) threshold_vec[!eligible] <- NA_real_
  # What a clipped count becomes: the bound itself for the legacy rule,
  # the largest whole number not above it otherwise.
  cap_vec <- if (legacy) threshold_vec else floor(threshold_vec)
  if (!legacy && is.integer(count_mat)) cap_vec <- as.integer(cap_vec)

  max_orig <- suppressWarnings(apply(count_mat, 1L, max, na.rm = TRUE))

  clipped_mat <- count_mat
  n_clipped_per_gene <- integer(n_genes)
  for (i in which(eligible)) {
    above <- observed[i, ] & count_mat[i, ] > threshold_vec[i]
    if (any(above)) {
      n_clipped_per_gene[i] <- sum(above)
      clipped_mat[i, above] <- cap_vec[i]
    }
  }

  stats_df <- data.frame(
    feature_id = rownames(count_mat) %||% rep(NA_character_, n_genes),
    q1 = round(q1_vec, 2),
    q3 = round(q3_vec, 2),
    iqr = round(iqr_vec, 2),
    threshold = round(threshold_vec, 2),
    max_original = max_orig,
    n_clipped = n_clipped_per_gene,
    n_nonzero = n_nonzero,
    winsorized = eligible,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  list(
    count_mat = clipped_mat,
    stats = stats_df,
    n_clipped = sum(n_clipped_per_gene),
    n_genes_affected = sum(n_clipped_per_gene > 0L),
    k = k,
    min_expressed = if (legacy) NA_integer_ else as.integer(min_expressed),
    legacy = legacy
  )
}
