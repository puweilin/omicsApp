#' Per-sample / per-feature missingness summary
#'
#' Computes per-sample and per-feature missing rates for an
#' [omics_input()] and flags samples / features whose missing rate exceeds
#' the supplied thresholds.
#'
#' The feature filter can look at the whole study or group by group. A
#' protein measured in every treated sample and in no control is missing
#' in half the samples overall, so the global rule puts it right at the
#' usual 50% cutoff, although it is the clearest on/off result the study
#' has. The group rules, standard in proteomics, ask instead whether the
#' feature was measured well enough *within* a group:
#'
#' * `"global"` (default) — missing rate over all samples.
#' * `"any_group"` — kept when at least one group is at or below the
#'   cutoff, i.e. observed in at least `1 - feature_missing_cutoff` of
#'   the samples of some group.
#' * `"all_groups"` — kept only when every group is at or below the
#'   cutoff.
#'
#' Samples with no value in `group_col` belong to no group and do not
#' count towards any group's rate.
#'
#' @param input An `omics_input`.
#' @param sample_missing_cutoff Optional sample missing-rate threshold in
#'   `[0, 1]`. `NULL` (default) leaves all samples unflagged.
#' @param feature_missing_cutoff Feature missing-rate threshold in `[0, 1]`.
#' @param missing_filter How the feature cutoff is applied: `"global"`,
#'   `"any_group"` or `"all_groups"` (see above).
#' @param group_col Sample-metadata column holding the groups, for the
#'   group rules. `NULL` uses the layer's recorded study design
#'   ([study_design()]).
#'
#' @return A list with:
#'   \describe{
#'     \item{`sample_metrics`}{`data.frame` with `sample_id`, `missing_rate`.}
#'     \item{`feature_metrics`}{`data.frame` with `feature_id`, `missing_rate`
#'       and, for the group rules, `filter_missing_rate`: the group rate the
#'       cutoff was compared with (the lowest for `"any_group"`, the
#'       highest for `"all_groups"`).}
#'     \item{`flagged_samples`}{Character vector of sample IDs exceeding
#'       `sample_missing_cutoff`.}
#'     \item{`flagged_features`}{Character vector of feature IDs exceeding
#'       `feature_missing_cutoff`.}
#'     \item{`settings`}{Echo of the thresholds and the filter used (and the
#'       group column, for the group rules).}
#'   }
#' @export
#' @family qc
qc_missingness <- function(
  input,
  sample_missing_cutoff = NULL,
  feature_missing_cutoff = 0.5,
  missing_filter = c("global", "any_group", "all_groups"),
  group_col = NULL
) {
  assert_number(sample_missing_cutoff, "sample_missing_cutoff",
                lower = 0, upper = 1, allow_null = TRUE)
  assert_number(feature_missing_cutoff, "feature_missing_cutoff", lower = 0, upper = 1)
  missing_filter <- match.arg(missing_filter)
  validate_omics_input(input)
  expr_mat <- input$expr_mat
  if (!identical(missing_filter, "global")) {
    group_col <- resolve_missing_group_col(input, group_col, missing_filter)
  }

  sample_metrics <- data.frame(
    sample_id = colnames(expr_mat),
    missing_rate = colMeans(is.na(expr_mat)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  feature_metrics <- data.frame(
    feature_id = rownames(expr_mat),
    missing_rate = rowMeans(is.na(expr_mat)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  flagged_samples <- if (is.null(sample_missing_cutoff)) {
    character(0)
  } else {
    sample_metrics$sample_id[sample_metrics$missing_rate > sample_missing_cutoff]
  }

  settings <- list(
    sample_missing_cutoff = sample_missing_cutoff,
    feature_missing_cutoff = feature_missing_cutoff,
    missing_filter = missing_filter
  )
  if (identical(missing_filter, "global")) {
    filter_rate <- feature_metrics$missing_rate
  } else {
    grp <- as.character(input$meta_df[colnames(expr_mat), group_col])
    levels <- unique(grp[!is.na(grp)])
    by_group <- matrix(
      vapply(levels, function(g) {
        rowMeans(is.na(expr_mat[, which(grp == g), drop = FALSE]))
      }, numeric(nrow(expr_mat))),
      nrow = nrow(expr_mat))
    filter_rate <- if (identical(missing_filter, "any_group")) {
      apply(by_group, 1L, min)
    } else {
      apply(by_group, 1L, max)
    }
    feature_metrics$filter_missing_rate <- unname(filter_rate)
    settings$group_col <- group_col
  }

  flagged_features <- feature_metrics$feature_id[filter_rate > feature_missing_cutoff]

  list(
    sample_metrics = sample_metrics,
    feature_metrics = feature_metrics,
    flagged_samples = flagged_samples,
    flagged_features = flagged_features,
    settings = settings
  )
}

# The group column for a group-wise missing filter: as given, else the
# layer's recorded design. Never guessed from the column names: filtering
# by the wrong column would quietly keep or drop the wrong features.
resolve_missing_group_col <- function(input, group_col, missing_filter) {
  assert_string(group_col, "group_col", allow_null = TRUE)
  if (is.null(group_col)) group_col <- study_design(input)$group_col
  if (is.null(group_col)) {
    stop(sprintf(paste(
      "missing_filter = '%s' needs `group_col`, or a study design recorded",
      "on the layer with set_study_design()."), missing_filter), call. = FALSE)
  }
  if (!group_col %in% names(input$meta_df)) {
    stop("`group_col` '", group_col, "' is not a column of the sample metadata.",
         call. = FALSE)
  }
  if (all(is.na(input$meta_df[[group_col]]))) {
    stop("`group_col` '", group_col, "' has no values.", call. = FALSE)
  }
  group_col
}
