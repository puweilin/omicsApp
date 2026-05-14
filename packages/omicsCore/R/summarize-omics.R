#' Summarize an omics_input or omics_project
#'
#' Returns a one-row-per-input `tibble` reporting sample / feature
#' counts, missingness, and modality labels. This is the inspection
#' counterpart to [run_qc()] and is intended for quick pre-flight checks
#' in the Shiny app's import wizard.
#'
#' For an `omics_input`, returns a single-row `tibble`. For an
#' `omics_project`, returns one row per experiment with an extra `tag`
#' column.
#'
#' Columns: `tag` (project only), `omics_type`, `assay_type`,
#' `n_samples`, `n_features`, `n_missing`, `missing_pct`,
#' `n_meta_cols`, `feature_columns`.
#'
#' @param x An [`omics_input`][omics_input()] or
#'   [`omics_project`][is_omics_project()].
#'
#' @return A `tibble::tibble`.
#' @export
#' @family omics_input
summarize_omics <- function(x) {
  if (inherits(x, "omics_input")) {
    return(tibble::as_tibble(summarize_omics_input(x)))
  }
  if (is_omics_project(x)) {
    tags <- experiment_tags(x)
    if (length(tags) == 0L) {
      return(tibble::tibble(
        tag = character(0),
        omics_type = character(0),
        assay_type = character(0),
        n_samples = integer(0),
        n_features = integer(0),
        n_missing = integer(0),
        missing_pct = numeric(0),
        n_meta_cols = integer(0),
        feature_columns = character(0)
      ))
    }
    rows <- lapply(tags, function(t) {
      row <- summarize_omics_input(x$experiments[[t]])
      cbind(data.frame(tag = t, stringsAsFactors = FALSE), row)
    })
    return(tibble::as_tibble(do.call(rbind, rows)))
  }
  stop("`x` must be an `omics_input` or an `omics_project`.")
}

# ---- internal helpers --------------------------------------------------

summarize_omics_input <- function(input) {
  mat <- input$expr_mat
  if (is.null(mat)) {
    n_total <- 0L
    n_missing <- 0L
    n_samples <- 0L
    n_features <- 0L
  } else {
    n_samples <- ncol(mat)
    n_features <- nrow(mat)
    n_total <- length(mat)
    n_missing <- sum(is.na(mat))
  }
  missing_pct <- if (n_total > 0L) 100 * n_missing / n_total else 0
  meta_cols <- if (is.null(input$meta_df)) 0L else ncol(input$meta_df)
  feat_cols <- if (is.null(input$feature_df)) character(0L) else colnames(input$feature_df)

  data.frame(
    omics_type = input$omics_type %||% NA_character_,
    assay_type = input$assay_type %||% NA_character_,
    n_samples = n_samples,
    n_features = n_features,
    n_missing = n_missing,
    missing_pct = round(missing_pct, 3),
    n_meta_cols = as.integer(meta_cols),
    feature_columns = paste(feat_cols, collapse = ","),
    stringsAsFactors = FALSE
  )
}
