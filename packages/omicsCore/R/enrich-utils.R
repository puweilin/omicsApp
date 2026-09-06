#' Filter standardized enrichment results
#'
#' Subsets a standardized enrichment `data.frame` (the schema returned by
#' [run_enrichment()]) to pathways that pass a significance cutoff on either
#' the adjusted, raw, or q-value column, with an optional minimum gene-set
#' size and direction filter.
#'
#' @param enrich_df Standardized enrichment `data.frame`.
#' @param p_cutoff Significance cutoff.
#' @param p_preference One of `"adjusted"` (default), `"raw"`, or `"qvalue"`.
#' @param min_genes Optional minimum number of overlapping / leading genes.
#' @param direction Optional direction filter (`"up"` or `"down"`). Useful
#'   for GSEA where pathways carry a sign.
#'
#' @return Filtered standardized enrichment `data.frame`.
#' @export
#' @family enrich
filter_enrich_results <- function(
  enrich_df,
  p_cutoff = 0.05,
  p_preference = c("adjusted", "raw", "qvalue"),
  min_genes = NULL,
  direction = NULL
) {
  p_preference <- match.arg(p_preference)
  assert_number(p_cutoff, "p_cutoff", lower = 0, upper = 1)
  assert_count(min_genes, "min_genes", allow_null = TRUE)

  if (!is.data.frame(enrich_df)) {
    if (methods::is(enrich_df, "enrichResult") || methods::is(enrich_df, "gseaResult")) {
      enrich_df <- as.data.frame(enrich_df)
    } else {
      stop("`enrich_df` must be a data.frame, enrichResult, or gseaResult.")
    }
  }

  p_col <- resolve_enrich_p_col(enrich_df, p_preference)
  out <- enrich_df[!is.na(enrich_df[[p_col]]) & enrich_df[[p_col]] < p_cutoff, , drop = FALSE]

  if (!is.null(min_genes)) {
    candidate_cols <- intersect(
      c("overlap_size", "gene_set_size", "Count", "setSize"),
      colnames(out)
    )
    if (length(candidate_cols) > 0L) {
      gc_col <- candidate_cols[[1L]]
      out <- out[!is.na(out[[gc_col]]) & out[[gc_col]] >= min_genes, , drop = FALSE]
    }
  }

  if (!is.null(direction)) {
    direction <- match.arg(direction, choices = c("up", "down"))
    if ("direction" %in% colnames(out)) {
      out <- out[!is.na(out$direction) & out$direction == direction, , drop = FALSE]
    }
  }

  rownames(out) <- NULL
  out
}

# ---- internal helpers --------------------------------------------------

resolve_enrich_p_col <- function(enrich_df, p_preference) {
  candidates <- switch(
    p_preference,
    adjusted = c("adj_p_value", "p.adjust"),
    raw      = c("p_value", "pvalue"),
    qvalue   = c("q_value", "qvalue")
  )
  hit <- intersect(candidates, colnames(enrich_df))
  if (length(hit) == 0L) {
    stop("No column found for p_preference = '", p_preference,
         "'. Looked for: ", paste(candidates, collapse = ", "))
  }
  hit[[1L]]
}
