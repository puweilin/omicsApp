# Proteomics normalization. RNA-seq does not go through here: DESeq2 and edgeR
# model raw counts directly, and the generic t-test / lm backends apply
# log2(x + 1) themselves when `assay_type` is "raw_count".

ensure_vsn <- function() {
  if (!requireNamespace("vsn", quietly = TRUE)) {
    stop(
      "Package 'vsn' is required for `method = \"vsn\"`. ",
      "Install with: omicsCore::install_optional('proteomics').",
      call. = FALSE
    )
  }
}

#' Normalize a proteomics `omics_input`
#'
#' Turns linear intensities into analysis-ready values and records that in
#' `assay_type`. Every differential backend reads `expr_mat` as-is -- limma
#' applies no transform at all -- so this step is what makes the numbers
#' comparable across samples. Without it, limma runs on raw instrument output:
#' heavily right-skewed, variance scaling with the mean, and fold changes that
#' are differences of linear intensities rather than log ratios.
#'
#' `"vsn"` is the method the legacy proteomics pipeline has always used. It
#' fits a variance-stabilising transform on the linear matrix and returns
#' glog2-scale values, which behave like log2 intensities at the high end while
#' staying finite near zero.
#'
#' The original matrix is kept in `raw_mat`, so the un-normalized values remain
#' available for QC views.
#'
#' @param input An `omics_input` with `omics_type = "proteomics"` carrying
#'   linear intensities.
#' @param method `"vsn"` (default) or `"log2"`.
#' @param offset Added before the log for `method = "log2"`, so that zeros
#'   survive. Ignored by `"vsn"`.
#'
#' @return The `omics_input` with `expr_mat` normalized, `raw_mat` holding the
#'   input matrix, and `assay_type` set to `"normalized_intensity"`.
#' @export
#' @family omics_input
#' @examples
#' expr <- matrix(2^rnorm(200, 20, 2), nrow = 50, dimnames = list(
#'   paste0("P", 1:50), paste0("s", 1:4)
#' ))
#' meta <- data.frame(group = c("A", "A", "B", "B"), row.names = paste0("s", 1:4))
#' feat <- data.frame(feature_id = paste0("P", 1:50), row.names = paste0("P", 1:50))
#' input <- omics_input(expr, meta, feat, omics_type = "proteomics",
#'                      assay_type = "raw_intensity")
#' normalized <- normalize_omics(input, method = "log2")
#' normalized$assay_type
normalize_omics <- function(
  input,
  method = c("vsn", "log2"),
  offset = 1
) {
  validate_omics_input(input)
  method <- match.arg(method)

  if (!identical(input$omics_type, "proteomics")) {
    stop(
      "`normalize_omics()` covers proteomics only; got omics_type = '",
      input$omics_type, "'. RNA-seq counts are handled by the analysis ",
      "backends: DESeq2 and edgeR model raw counts directly, and the ",
      "t-test / lm backends log2-transform `raw_count` themselves.",
      call. = FALSE
    )
  }

  # Refusing here is the point of the whole contract: normalizing twice
  # compresses the dynamic range and shrinks every fold change, without
  # erroring anywhere downstream.
  if (input$assay_type %in% LOG_SCALE_ASSAY_TYPES) {
    stop(
      "`assay_type` is already '", input$assay_type, "', so these values have ",
      "been normalized. Normalizing again would compress the dynamic range ",
      "and shrink every fold change. If the label is wrong, rebuild the input ",
      "with the correct `assay_type`.",
      call. = FALSE
    )
  }

  mat <- input$expr_mat
  storage.mode(mat) <- "double"

  n_zero <- sum(mat == 0, na.rm = TRUE)
  n_negative <- sum(mat < 0, na.rm = TRUE)
  if (n_negative > 0) {
    warning(
      n_negative, " negative value(s) treated as missing: an intensity below ",
      "zero has no physical meaning and neither transform is defined there.",
      call. = FALSE
    )
  }
  # A zero intensity means "not detected", not "zero abundance"
  mat[mat == 0] <- NA
  mat[mat < 0] <- NA

  normalized <- switch(
    method,
    vsn = {
      ensure_vsn()
      fit <- vsn::vsnMatrix(mat)
      vsn::predict(fit, mat)
    },
    log2 = log2(mat + offset)
  )
  dimnames(normalized) <- dimnames(input$expr_mat)

  message(
    "Normalized ", nrow(normalized), " features x ", ncol(normalized),
    " samples with ", method,
    if (n_zero > 0) paste0(" (", n_zero, " zero(s) set to NA)") else ""
  )

  out <- input
  out$expr_mat <- normalized
  # Fill raw_mat only if the caller was not already carrying one
  out$raw_mat <- input$raw_mat %||% input$expr_mat
  out$normalized_mat <- normalized
  out$assay_type <- "normalized_intensity"

  validate_omics_input(out)
  out
}
