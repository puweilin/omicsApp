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

# Which gene list an ORA tested, in words, for a plot subtitle: "" for
# GSEA, or for a bundle that does not record it.
ora_list_caption <- function(params) {
  if (!identical(params$type, "ora") || is.null(params$direction)) return("")
  switch(params$direction,
         separate = " · up- and down-regulated genes tested separately",
         up = " · up-regulated genes",
         down = " · down-regulated genes",
         both = " · up and down pooled",
         "")
}

# ---- gene-symbol case -------------------------------------------------
#
# Gene sets are matched to the data by symbol, exactly. Symbol case is a
# species convention -- TP53 in human, Trp53 in mouse and rat, tp53 in
# zebrafish -- and data often break it: a mouse table exported in upper
# case, a fish table run through a tool that capitalised everything.
# Matched as written, such a table finds a handful of genes and the
# result is an empty or near-empty enrichment with nothing to say why.
#
# So when the exact match is poor and ignoring case would match clearly
# more genes, the data's symbols are rewritten to the gene sets' spelling
# before enrichment, and the run says so. Only unambiguous matches are
# used: a set symbol that differs from another only by case (rare, but
# it happens in some fly and yeast tables) is never a target.

# How the data's `symbols` would match the gene sets' `reference`.
# `apply` is TRUE when the case-insensitive match gains at least
# `min_gain` genes and the exact match finds under half of what the
# case-insensitive one does.
symbol_case_plan <- function(symbols, reference, min_gain = 10L) {
  u <- unique(as.character(symbols[!is.na(symbols) & nzchar(symbols)]))
  ref <- unique(as.character(reference[!is.na(reference)]))
  up_ref <- toupper(ref)
  ambiguous <- unique(up_ref[duplicated(up_ref)])
  keep <- !up_ref %in% ambiguous
  map <- stats::setNames(ref[keep], up_ref[keep])
  exact <- u %in% ref
  ci <- !exact & toupper(u) %in% names(map)
  n_exact <- sum(exact)
  n_ci <- n_exact + sum(ci)
  list(
    apply = sum(ci) >= min_gain && n_exact < 0.5 * n_ci,
    map = map,
    reference = ref,
    n_symbols = length(u),
    n_exact = n_exact,
    n_ci = n_ci
  )
}

# Rewrite `x` to the gene sets' spelling where only case differs.
apply_symbol_case <- function(x, plan) {
  x <- as.character(x)
  idx <- !is.na(x) & !(x %in% plan$reference) & toupper(x) %in% names(plan$map)
  x[idx] <- unname(plan$map[toupper(x[idx])])
  x
}

# Put the data's own spelling back into "/"-joined gene lists, so the
# genes a result names are the ones the user's table has.
restore_symbol_case <- function(gene_lists, back) {
  if (!length(back) || !length(gene_lists)) return(gene_lists)
  vapply(gene_lists, function(s) {
    if (is.na(s) || !nzchar(s)) return(s)
    g <- strsplit(s, "/", fixed = TRUE)[[1L]]
    hit <- g %in% names(back)
    g[hit] <- unname(back[g[hit]])
    paste(g, collapse = "/")
  }, character(1), USE.NAMES = FALSE)
}

# ---- multiple testing across databases --------------------------------

# Re-adjust p-values across every database of one result instead of
# within each. Each gene list is its own family: for ORA the up-, down-
# or pooled list (`direction`, NA for pooled), for GSEA the one ranking.
# Rows are then bounded at `cutoff` on raw and adjusted p alike, as
# clusterProfiler bounds its own tables. q-values are per-database
# quantities this does not recompute, so they are set to NA rather than
# left contradicting the new adjusted p.
adjust_enrich_across_databases <- function(df, type, p_adjust_method, cutoff) {
  if (!nrow(df)) return(df)
  family <- if (identical(type, "ora")) {
    ifelse(is.na(df$direction), "pooled", df$direction)
  } else {
    rep("all", nrow(df))
  }
  adj <- df$adj_p_value
  for (f in unique(family)) {
    i <- family == f
    adj[i] <- stats::p.adjust(df$p_value[i], method = p_adjust_method)
  }
  df$adj_p_value <- adj
  df$q_value <- NA_real_
  keep <- !is.na(df$p_value) & df$p_value <= cutoff &
    !is.na(df$adj_p_value) & df$adj_p_value <= cutoff
  out <- df[keep, , drop = FALSE]
  rownames(out) <- NULL
  out
}
