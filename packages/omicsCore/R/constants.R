#' Package constants
#'
#' Internal vocabulary used across `omicsCore`. Not exported.
#' @keywords internal
#' @name omicsCore-constants
NULL

#' @rdname omicsCore-constants
SUPPORTED_OMICS_TYPES <- c("proteomics", "rnaseq")

#' @rdname omicsCore-constants
SUPPORTED_DIFF_ANALYSIS_TYPES <- c("group", "continuous", "anova")

#' @rdname omicsCore-constants
SUPPORTED_PREFERENCE <- c("adjusted", "raw")

#' @rdname omicsCore-constants
SUPPORTED_ENRICH_PREFERENCE <- c("adjusted", "raw", "qvalue")

#' @rdname omicsCore-constants
DEFAULT_CACHE_MAX_AGE_HOURS <- 24 * 30

#' @rdname omicsCore-constants
CANONICAL_DATABASES <- c("GO_BP", "GO_MF", "GO_CC", "KEGG", "Reactome", "Wiki")
