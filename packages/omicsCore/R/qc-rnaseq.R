#' Winsorize per-gene outliers in a count matrix
#'
#' For each gene (row), values exceeding `Q3 + k * IQR` are clipped to that
#' threshold. This removes extreme leverage points that can drive spurious
#' differential expression results while preserving the overall distribution.
#'
#' Pure base R; no `matrixStats` dependency.
#'
#' @param count_mat Integer or numeric count matrix (genes × samples).
#' @param k IQR multiplier for the upper fence. Larger values are more
#'   conservative (fewer values clipped). Default `20`.
#'
#' @return A list with components:
#'   \item{`count_mat`}{Winsorized count matrix (same dimensions / names).}
#'   \item{`stats`}{`data.frame` with per-gene outlier statistics.}
#'   \item{`n_clipped`}{Total number of values clipped.}
#'   \item{`n_genes_affected`}{Number of genes with at least one clipped value.}
#'   \item{`k`}{The multiplier used.}
#' @export
#' @family qc
winsorize_counts <- function(count_mat, k = 20) {
  count_mat <- assert_numeric_matrix(count_mat, "count_mat")
  assert_number(k, "k", lower = 0)
  n_genes <- nrow(count_mat)

  qmat <- apply(
    count_mat, 1L,
    function(x) stats::quantile(x, probs = c(0.25, 0.75), na.rm = TRUE),
    simplify = TRUE
  )
  # qmat is 2 x n_genes
  q1_vec <- qmat[1L, ]
  q3_vec <- qmat[2L, ]
  iqr_vec <- q3_vec - q1_vec
  threshold_vec <- q3_vec + k * iqr_vec

  max_orig <- apply(count_mat, 1L, max, na.rm = TRUE)

  clipped_mat <- count_mat
  n_clipped_per_gene <- integer(n_genes)
  for (i in seq_len(n_genes)) {
    above <- !is.na(count_mat[i, ]) & count_mat[i, ] > threshold_vec[i]
    if (any(above)) {
      n_clipped_per_gene[i] <- sum(above)
      clipped_mat[i, above] <- threshold_vec[i]
    }
  }

  stats_df <- data.frame(
    feature_id = rownames(count_mat),
    q1 = round(q1_vec, 2),
    q3 = round(q3_vec, 2),
    iqr = round(iqr_vec, 2),
    threshold = round(threshold_vec, 2),
    max_original = max_orig,
    n_clipped = n_clipped_per_gene,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  list(
    count_mat = clipped_mat,
    stats = stats_df,
    n_clipped = sum(n_clipped_per_gene),
    n_genes_affected = sum(n_clipped_per_gene > 0L),
    k = k
  )
}
