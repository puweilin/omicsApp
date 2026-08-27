#' Package constants
#'
#' Internal vocabulary used across `omicsCore`. Not exported.
#' @keywords internal
#' @name omicsCore-constants
NULL

#' @rdname omicsCore-constants
SUPPORTED_OMICS_TYPES <- c("proteomics", "rnaseq")

#' Assay types recognised per omics modality
#'
#' A named list with one entry per [SUPPORTED_OMICS_TYPES] value. The names are
#' deliberately explicit about scale, because nothing downstream re-derives it:
#' `raw_intensity` is linear instrument output, `normalized_intensity` has been
#' variance-stabilised (vsn) and is log-like, `raw_count` is untransformed
#' RNA-seq counts that DESeq2/edgeR model directly.
#'
#' Values outside this vocabulary are allowed -- `validate_omics_input()` only
#' warns -- so callers can drive omicsCore with modalities it does not ship
#' dispatchers for (CyTOF `arcsinh_intensity`, for example).
#'
#' @export
#' @family omics_input
SUPPORTED_ASSAY_TYPES <- list(
  proteomics = c(
    "raw_intensity",
    "normalized_intensity",
    "imputed_intensity",
    "filtered_intensity"
  ),
  rnaseq = c(
    "raw_count",
    "tpm",
    "fpkm",
    "vst",
    "logcpm"
  )
)

#' Assay types that already sit on a log-like scale
#'
#' Used by [check_assay_scale()] to decide which direction a scale mismatch
#' points in, and by callers deciding whether a layer still needs
#' [normalize_omics()]. Everything else in [SUPPORTED_ASSAY_TYPES] is linear.
#'
#' @export
#' @family omics_input
LOG_SCALE_ASSAY_TYPES <- c(
  "normalized_intensity",
  "imputed_intensity",
  "filtered_intensity",
  "vst",
  "logcpm"
)

#' Superseded assay-type spellings
#'
#' Names are the old values, elements the replacement. `omics_input()` rewrites
#' them; `validate_omics_input()` accepts them with a warning so projects saved
#' before the vocabulary existed still load.
#'
#' `"intensity"` was what the import view stamped on every proteomics upload,
#' and `"raw_counts"` appeared in the documentation while the code only ever
#' matched `"raw_count"`.
#'
#' @rdname omicsCore-constants
DEPRECATED_ASSAY_TYPE_ALIASES <- c(
  "intensity"  = "raw_intensity",
  "raw_counts" = "raw_count"
)

#' Largest value still plausible on a log scale
#'
#' log2 proteomics intensities top out near 35; linear intensities run into the
#' millions. [check_assay_scale()] uses this to tell them apart.
#'
#' @rdname omicsCore-constants
MAX_PLAUSIBLE_LOG_SCALE_VALUE <- 100

#' Assay types whose scale can be inferred from magnitude
#'
#' [check_assay_scale()] separates linear from log values by how large they
#' get, which only works for intensity-like assays. Counts, TPM and FPKM are
#' legitimately small -- a gene with five reads is ordinary -- so applying the
#' same rule to them would flag correct data.
#'
#' @rdname omicsCore-constants
SCALE_CHECKED_ASSAY_TYPES <- c(
  "raw_intensity",
  "normalized_intensity",
  "imputed_intensity",
  "filtered_intensity"
)

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
