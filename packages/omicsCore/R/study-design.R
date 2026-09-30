# Which metadata column holds the groups, and which group is the control.
#
# Decided once, when a layer is imported -- the person importing knows
# their study -- and read by every view after that, so QC colours by the
# same column the Differential view compares on, with the same reference.
# Before this each view guessed on its own: QC coloured by a column
# literally called `group`, Differential took the alphabetically first
# level as the control.

#' The study design recorded on a layer
#'
#' @param input An [omics_input()].
#' @return `list(group_col, reference)`, or `NULL` when none is recorded
#'   (or the recorded column is no longer in the metadata).
#' @export
#' @family omics_input
study_design <- function(input) {
  if (!is_omics_input(input)) stop("`input` must be an `omics_input`.", call. = FALSE)
  d <- input$design
  if (is.null(d) || is.null(d$group_col) || !d$group_col %in% names(input$meta_df)) {
    return(NULL)
  }
  levels <- unique(as.character(stats::na.omit(input$meta_df[[d$group_col]])))
  ref <- d$reference
  if (!is.null(ref) && !ref %in% levels) ref <- NULL
  list(group_col = d$group_col, reference = ref)
}

#' Record the study design on a layer
#'
#' @param input An [omics_input()].
#' @param group_col The sample-metadata column holding the groups, or
#'   `NULL` to clear the record.
#' @param reference The control / reference level of `group_col`.
#'   Optional.
#' @return `input`, with the design recorded.
#' @export
#' @family omics_input
#' @examples
#' expr <- matrix(1:12, 3, dimnames = list(paste0("g", 1:3), paste0("s", 1:4)))
#' meta <- data.frame(treatment = c("DMSO", "DMSO", "Drug", "Drug"),
#'                    row.names = paste0("s", 1:4))
#' feat <- data.frame(feature_id = paste0("g", 1:3), row.names = paste0("g", 1:3))
#' x <- omics_input(expr, meta, feat, omics_type = "proteomics",
#'                  assay_type = "normalized_intensity")
#' study_design(set_study_design(x, "treatment", "DMSO"))
set_study_design <- function(input, group_col, reference = NULL) {
  if (!is_omics_input(input)) stop("`input` must be an `omics_input`.", call. = FALSE)
  if (is.null(group_col)) {
    input$design <- NULL
    return(input)
  }
  assert_string(group_col, "group_col")
  assert_label(reference, "reference", allow_null = TRUE)
  if (!group_col %in% names(input$meta_df)) {
    stop("`group_col` '", group_col, "' is not a column of the sample metadata.",
         call. = FALSE)
  }
  levels <- unique(as.character(stats::na.omit(input$meta_df[[group_col]])))
  if (!is.null(reference) && !as.character(reference) %in% levels) {
    stop("`reference` '", reference, "' is not a level of '", group_col,
         "'. Levels: ", paste(levels, collapse = ", "), ".", call. = FALSE)
  }
  input$design <- list(group_col = group_col,
                       reference = if (!is.null(reference)) as.character(reference))
  input
}
