# Standardize enrichment results into the package-wide schema defined by
# ENRICH_RESULT_REQUIRED_COLS. Every backend (ORA, GSEA, GSVA association)
# must funnel through these helpers so downstream filtering, plotting, and
# integration code can rely on a single column set.

parse_overlap_size <- function(x) {
  if (is.null(x) || length(x) == 0L || is.na(x)) return(NA_real_)
  pieces <- strsplit(as.character(x), "/", fixed = TRUE)[[1L]]
  if (length(pieces) < 1L) return(NA_real_)
  suppressWarnings(as.numeric(pieces[[1L]]))
}

# Coerce an `enrichResult` / `gseaResult` / `data.frame` into a flat data
# frame, then re-pack into the standard schema. NULL or empty inputs return
# an empty template.
standardize_enrich_result <- function(
  enrich_obj,
  database,
  result_type,
  comparison,
  source_label = NULL,
  default_direction = NA_character_
) {
  template <- new_enrich_result_template()
  if (is.null(enrich_obj)) return(template)

  enrich_df <- if (methods::is(enrich_obj, "enrichResult") ||
                   methods::is(enrich_obj, "gseaResult")) {
    as.data.frame(enrich_obj)
  } else if (is.data.frame(enrich_obj)) {
    enrich_obj
  } else {
    stop("`enrich_obj` must be NULL, a data.frame, an enrichResult, or a gseaResult.")
  }

  if (nrow(enrich_df) == 0L) return(template)

  n <- nrow(enrich_df)
  out <- as.data.frame(
    lapply(template, function(col) vector(mode = typeof(col), length = n)),
    stringsAsFactors = FALSE
  )

  has_nes <- "NES" %in% colnames(enrich_df)

  out$database <- rep(database, n)
  out$result_type <- rep(result_type, n)
  out$comparison <- rep(comparison, n)
  out$pathway_id <- as.character(
    enrich_df$ID %||% enrich_df$gs_id %||% enrich_df$Description
  )
  out$pathway_name <- as.character(
    enrich_df$Description %||% enrich_df$ID %||% enrich_df$gs_name
  )
  out$effect <- if (has_nes) {
    as.numeric(enrich_df$NES)
  } else if ("logFC" %in% colnames(enrich_df)) {
    as.numeric(enrich_df$logFC)
  } else {
    rep(NA_real_, n)
  }
  out$effect_type <- if (has_nes) {
    rep("nes", n)
  } else if ("logFC" %in% colnames(enrich_df)) {
    rep("logFC", n)
  } else {
    rep(NA_character_, n)
  }
  out$direction <- if (has_nes) {
    ifelse(enrich_df$NES >= 0, "up", "down")
  } else {
    rep(default_direction, n)
  }
  out$p_value <- if ("pvalue" %in% colnames(enrich_df)) {
    as.numeric(enrich_df$pvalue)
  } else {
    rep(NA_real_, n)
  }
  out$adj_p_value <- if ("p.adjust" %in% colnames(enrich_df)) {
    as.numeric(enrich_df$p.adjust)
  } else {
    rep(NA_real_, n)
  }
  out$q_value <- if ("qvalue" %in% colnames(enrich_df)) {
    as.numeric(enrich_df$qvalue)
  } else {
    rep(NA_real_, n)
  }
  out$gene_set_size <- if ("setSize" %in% colnames(enrich_df)) {
    as.numeric(enrich_df$setSize)
  } else if ("Count" %in% colnames(enrich_df)) {
    as.numeric(enrich_df$Count)
  } else {
    rep(NA_real_, n)
  }
  out$overlap_size <- if ("Count" %in% colnames(enrich_df)) {
    as.numeric(enrich_df$Count)
  } else if ("GeneRatio" %in% colnames(enrich_df)) {
    vapply(enrich_df$GeneRatio, parse_overlap_size, numeric(1))
  } else {
    rep(NA_real_, n)
  }
  out$overlap_features <- if ("geneID" %in% colnames(enrich_df)) {
    as.character(enrich_df$geneID)
  } else {
    rep(NA_character_, n)
  }
  out$leading_features <- if ("core_enrichment" %in% colnames(enrich_df)) {
    as.character(enrich_df$core_enrichment)
  } else {
    out$overlap_features
  }
  out$source_label <- rep(source_label %||% result_type, n)
  out
}

# Concatenate per-database results into one standardized data frame.
standardize_enrich_result_list <- function(
  enrich_list,
  result_type,
  comparison,
  default_direction = NA_character_,
  source_prefix = NULL
) {
  if (length(enrich_list) == 0L) return(new_enrich_result_template())
  rows <- lapply(names(enrich_list), function(db) {
    standardize_enrich_result(
      enrich_obj = enrich_list[[db]],
      database = db,
      result_type = result_type,
      comparison = comparison,
      source_label = if (is.null(source_prefix)) result_type else paste0(source_prefix, "_", db),
      default_direction = default_direction
    )
  })
  dplyr::bind_rows(rows)
}

`%||%` <- function(x, y) if (is.null(x)) y else x
