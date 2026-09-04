# Working out which sample in one layer is the same person as which
# sample in another.
#
# Two layers of the same study rarely share sample ids, and should not:
# `RD001-C` is a cheek biopsy and `RD001_Folli` is a hair follicle. They
# are different samples. What they share is the donor.
#
# build_sample_pairs() already knows how to use a `sample_link`; what was
# missing was any way to produce one. Two routes, in order of confidence:
#
#   1. A `donor` column in each layer's metadata. Stated by whoever knows
#      the study, so it is taken as given.
#   2. A guess from the ids themselves, when they share a leading stem.
#      Offered as a suggestion for someone to confirm, never applied
#      silently -- a wrong pairing produces correlations between
#      different people, which look like results.
#
# Route 2 exists because the alternative is a trap: a user who imported
# without a donor column would otherwise have to re-import and re-run
# everything to get integration, and nothing else in the analysis depends
# on donor at all.

DONOR_COLUMN_NAMES <- c("donor", "donor_id", "subject", "subject_id",
                        "patient", "patient_id", "individual")

#' Find the column naming the person a sample came from
#'
#' @param meta_df Sample metadata, one row per sample.
#' @return Column name, or `NULL`.
#' @keywords internal
#' @noRd
donor_column <- function(meta_df) {
  if (is.null(meta_df) || ncol(meta_df) == 0L) return(NULL)
  hit <- intersect(tolower(colnames(meta_df)), DONOR_COLUMN_NAMES)
  if (length(hit) == 0L) return(NULL)
  colnames(meta_df)[match(hit[[1L]], tolower(colnames(meta_df)))]
}

#' Everything before the first separator
#'
#' `RD001-C` and `RD001_Folli` both reduce to `RD001`. Deliberately blunt:
#' this only ever produces a suggestion someone confirms, and a rule
#' simple enough to explain is one a user can predict.
#'
#' @keywords internal
#' @noRd
sample_stem <- function(x) sub("[-_. ].*$", "", as.character(x))

#' Build a sample link from donor columns
#'
#' @param project An [omics_project].
#' @return A data frame with `tag`, `sample_id`, `donor_id`, or `NULL`
#'   when fewer than two layers name a donor.
#' @export
#' @family integrate
derive_sample_link <- function(project) {
  exps <- project$experiments
  if (length(exps) < 2L) return(NULL)

  rows <- list()
  for (tag in names(exps)) {
    meta <- exps[[tag]]$meta_df
    col <- donor_column(meta)
    if (is.null(col)) next
    ids <- rownames(meta)
    donors <- as.character(meta[[col]])
    keep <- !is.na(donors) & nzchar(donors)
    if (!any(keep)) next
    rows[[tag]] <- data.frame(tag = tag, sample_id = ids[keep],
                              donor_id = donors[keep],
                              stringsAsFactors = FALSE)
  }
  if (length(rows) < 2L) return(NULL)
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Suggest a sample link from the ids
#'
#' Only when the stems pair one-to-one. `A-1` and `A-2` both stem to `A`,
#' which would pair two samples of one layer to the same donor and pair
#' the wrong people in the other -- a suggestion worse than none, because
#' the result of acting on it still looks like a result.
#'
#' @param project An [omics_project].
#' @param tag_a,tag_b Experiment tags.
#' @return A data frame shaped like [derive_sample_link()], or `NULL`.
#' @export
#' @family integrate
suggest_sample_link <- function(project, tag_a, tag_b) {
  exps <- project$experiments
  if (!all(c(tag_a, tag_b) %in% names(exps))) return(NULL)

  a <- colnames(exps[[tag_a]]$expr_mat)
  b <- colnames(exps[[tag_b]]$expr_mat)
  if (is.null(a) || is.null(b)) return(NULL)

  stem_a <- sample_stem(a)
  stem_b <- sample_stem(b)
  if (anyDuplicated(stem_a) || anyDuplicated(stem_b)) return(NULL)

  shared <- intersect(stem_a, stem_b)
  if (length(shared) == 0L) return(NULL)

  out <- rbind(
    data.frame(tag = tag_a, sample_id = a[match(shared, stem_a)],
               donor_id = shared, stringsAsFactors = FALSE),
    data.frame(tag = tag_b, sample_id = b[match(shared, stem_b)],
               donor_id = shared, stringsAsFactors = FALSE)
  )
  rownames(out) <- NULL
  out
}

#' The pairing two layers would be integrated on
#'
#' What the Integration view shows before it runs anything. One row per
#' donor, one column per layer, and a `source` saying where the pairing
#' came from -- because "these two samples are the same person" is an
#' assertion the reader should be able to check rather than infer.
#'
#' @param project An [omics_project].
#' @param tag_a,tag_b Experiment tags.
#' @return A list with `pairs` (data frame, possibly zero rows) and
#'   `source`, one of `"linked"` (a pairing saved on the project),
#'   `"donor"` (a donor column in both layers), `"sample_id"` (ids that
#'   match outright), `"suggested"` (guessed from the ids) or `"none"`.
#' @export
#' @family integrate
sample_pairing_preview <- function(project, tag_a, tag_b) {
  exps <- project$experiments
  empty <- function(src) {
    list(pairs = data.frame(donor_id = character(0),
                            a = character(0), b = character(0),
                            stringsAsFactors = FALSE),
         source = src)
  }
  if (!all(c(tag_a, tag_b) %in% names(exps))) return(empty("none"))

  a <- colnames(exps[[tag_a]]$expr_mat)
  b <- colnames(exps[[tag_b]]$expr_mat)

  # Most confident first. A stated link beats a derived one, a derived
  # one beats ids that happen to match, and matching ids beat a guess
  # from their shape -- each step down is a weaker claim about two
  # samples being the same person, and the view labels which one was
  # used so a reader can judge it.
  link <- project$sample_link
  source <- "linked"
  if (is.null(link) || nrow(link) == 0L) {
    link <- derive_sample_link(project)
    source <- "donor"
  }

  if (is.null(link) || nrow(link) == 0L) {
    shared <- intersect(a, b)
    if (length(shared) > 0L) {
      return(list(
        pairs = data.frame(donor_id = shared, a = shared, b = shared,
                           stringsAsFactors = FALSE),
        source = "sample_id"
      ))
    }
    link <- suggest_sample_link(project, tag_a, tag_b)
    source <- "suggested"
  }

  if (is.null(link) || nrow(link) == 0L) return(empty("none"))

  la <- link[link$tag == tag_a, c("sample_id", "donor_id"), drop = FALSE]
  lb <- link[link$tag == tag_b, c("sample_id", "donor_id"), drop = FALSE]
  merged <- merge(la, lb, by = "donor_id", suffixes = c("_a", "_b"))
  merged <- merged[merged$sample_id_a %in% a & merged$sample_id_b %in% b, ,
                   drop = FALSE]
  if (nrow(merged) == 0L) return(empty(source))

  list(
    pairs = data.frame(donor_id = merged$donor_id,
                       a = merged$sample_id_a,
                       b = merged$sample_id_b,
                       stringsAsFactors = FALSE),
    source = source
  )
}
