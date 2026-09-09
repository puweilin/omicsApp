# A grouping the sample names already carry.
#
# Some vendors ship a matrix and a sample list and nothing else: the
# design is in the names -- `Drug1-1`, `Drug1-2`, `Drug1-3` -- and the
# metadata sheet has one column of identifiers. Differential analysis then
# has nothing to compare, and the view says "pick a group column with at
# least two levels" to a user who has none to pick.
#
# The rule is the one Integration already uses to pair donors: everything
# before the first separator. Deliberately blunt, and deliberately only a
# suggestion -- deriving an experimental design from a naming convention
# and acting on it unasked is asserting what the experiment was.

#' A grouping suggested by the sample names
#'
#' Applies `sample_stem()` -- the same rule that pairs donors across
#' layers, so a project cannot pair by one convention and group by
#' another -- to the sample ids and reports whether the
#' result could serve as a grouping: at least two levels, and at least
#' `min_per_level` samples in each. A design with singleton groups cannot
#' be tested, so proposing it would only move the dead end.
#'
#' @param input An `omics_input`.
#' @param min_per_level Samples a level needs to be usable.
#'
#' @return `NULL` when the names carry no usable grouping, otherwise a
#'   list with `groups` (named by sample id), `levels`, `n_levels` and
#'   `sizes`.
#' @export
#' @family qc
suggest_grouping_from_names <- function(input, min_per_level = 2L) {
  validate_omics_input(input)
  assert_count(min_per_level, "min_per_level", lower = 1L)
  ids <- colnames(input$expr_mat)
  if (length(ids) < 2L) return(NULL)

  stems <- sample_stem(ids)
  # A stem equal to the whole id means the name had no separator, so the
  # rule found nothing -- reporting that as a grouping of one sample each
  # would be a suggestion nobody can use.
  if (all(stems == ids)) return(NULL)
  sizes <- table(stems)
  if (length(sizes) < 2L || min(sizes) < min_per_level) return(NULL)
  # Every sample in one group is not a grouping either.
  if (length(sizes) == length(ids)) return(NULL)

  list(groups = stats::setNames(stems, ids),
       levels = names(sizes),
       n_levels = length(sizes),
       sizes = as.integer(sizes))
}

#' Write a suggested grouping into a layer's sample metadata
#'
#' Adds `column` to `meta_df` from [suggest_grouping_from_names()] and
#' records the step, so the exported script rebuilds the same column
#' rather than reading a file that never had it.
#'
#' Recorded rather than merely applied, for the reason
#' [relabel_assay_type()] is: a change made in the interface and kept only
#' there is one the script cannot replay, and nothing about the result
#' looks wrong while that is true.
#'
#' @param input An `omics_input`.
#' @param column Name for the new column.
#' @param min_per_level Passed to [suggest_grouping_from_names()].
#'
#' @return The input, with the column added.
#' @export
#' @family qc
add_grouping_from_names <- function(input, column = "group",
                                    min_per_level = 2L) {
  validate_omics_input(input)
  assert_string(column, "column")
  if (column %in% colnames(input$meta_df)) {
    stop("`meta_df` already has a column called '", column,
         "'. Pick another name, or use the column that is there.",
         call. = FALSE)
  }
  suggestion <- suggest_grouping_from_names(input, min_per_level)
  if (is.null(suggestion)) {
    stop("The sample names do not carry a usable grouping: every level ",
         "needs at least ", min_per_level, " samples.", call. = FALSE)
  }
  input$meta_df[[column]] <- unname(
    suggestion$groups[rownames(input$meta_df)])
  record_step(input, "add_grouping_from_names",
              list(column = column, min_per_level = min_per_level))
}
