# Internal helpers shared by the integration backends. Two cross-omics
# join strategies are supported:
#
#   * feature-level: join on `feature_symbol` (or `feature_id`) across the
#     two layers. Used by correlation + concordance, since the unit of
#     analysis is a feature pair.
#   * sample-level: align samples across layers either via the project
#     `sample_link$donor_id` map, or by direct sample-ID matching if the
#     two layers share IDs. Used by correlation.

# Resolve the two experiment tags. If `experiments` is NULL we expect the
# project to contain exactly two layers; if more, the caller must pick.
resolve_experiment_pair <- function(project, experiments) {
  if (!is_omics_project(project)) {
    stop("`project` must be an `omics_project`.")
  }
  tags <- experiment_tags(project)
  if (length(tags) < 2L) {
    stop("Integration requires a project with at least two experiments.")
  }
  if (is.null(experiments)) {
    if (length(tags) != 2L) {
      stop("`experiments` must name the two layers to integrate when the ",
           "project has more than two: ", paste(tags, collapse = ", "), ".")
    }
    return(tags)
  }
  if (!is.character(experiments) || length(experiments) != 2L) {
    stop("`experiments` must be a length-2 character vector of experiment tags.")
  }
  missing <- setdiff(experiments, tags)
  if (length(missing) > 0L) {
    stop("Experiments not found in project: ", paste(missing, collapse = ", "))
  }
  experiments
}

# Build a sample mapping data.frame between two experiments. Returns a
# `data.frame` with columns `donor_id`, `<tag_a>`, `<tag_b>`.
#
# Resolved exactly as the Integration view previews it
# (`sample_pairing_preview()`): a saved `sample_link`, then a donor column
# in both layers, then sample ids that match outright. A pairing that is
# only *guessed* from the shape of the ids is never used here -- it has
# to be accepted (saved as a `sample_link`) first. The two used to
# disagree: the view announced "12 pairs, from the donor column" and the
# run then failed with "No shared sample IDs", because this function only
# knew about a saved link and identical ids.
#
# A donor with more than one sample in a layer (technical replicates,
# several time points) cannot be paired one-to-one. Such donors keep
# their first sample in each layer; the number dropped is carried on the
# result as attribute `n_ambiguous` so the caller can report it rather
# than silently correlating one person's samples against themselves.
build_sample_pairs <- function(project, tag_a, tag_b) {
  prev <- sample_pairing_preview(project, tag_a, tag_b)
  if (identical(prev$source, "suggested")) {
    stop("The samples of '", tag_a, "' and '", tag_b, "' are only paired by ",
         "a guess from their ids. Accept the pairing (save it as the ",
         "project's `sample_link`) or add a donor column to both layers.",
         call. = FALSE)
  }
  pairs <- prev$pairs
  if (nrow(pairs) == 0L) {
    stop("No sample pairing between '", tag_a, "' and '", tag_b, "': ",
         "no `sample_link` covers both layers, no donor column is shared, ",
         "and no sample ids match.", call. = FALSE)
  }
  ambiguous <- duplicated(pairs$donor_id) | duplicated(pairs$a) |
    duplicated(pairs$b)
  n_ambiguous <- sum(ambiguous)
  pairs <- pairs[!ambiguous, , drop = FALSE]
  out <- data.frame(
    donor_id = pairs$donor_id,
    a = pairs$a,
    b = pairs$b,
    stringsAsFactors = FALSE
  )
  names(out)[2:3] <- c(tag_a, tag_b)
  rownames(out) <- NULL
  attr(out, "source") <- prev$source
  attr(out, "n_ambiguous") <- n_ambiguous
  out
}

# One row per join key. Where several features share a key (protein
# isoforms, several probes or protein groups naming one gene) the one kept
# is the most abundant -- a choice that does not look at the test result,
# unlike "the most significant", which would inflate agreement between
# layers by picking each gene's luckiest row. Ties, and rows without an
# abundance, fall back to their original order.
dedupe_by_key <- function(key, abundance = NULL) {
  n <- length(key)
  if (n == 0L) return(integer(0))
  if (is.null(abundance)) abundance <- rep(NA_real_, n)
  abundance <- suppressWarnings(as.numeric(abundance))
  ord <- order(is.na(abundance), -abundance, seq_len(n), na.last = TRUE)
  ord[!duplicated(key[ord])] |> sort()
}

# Case- and whitespace-insensitive join key. RNA and protein tables of one
# study routinely disagree on nothing but case (`Tp53` / `TP53`,
# UniProt gene names vs Ensembl symbols) or carry a trailing space from a
# spreadsheet; neither is a different gene.
integration_join_key <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[!is.na(x) & !nzchar(x)] <- NA_character_
  x
}

# Differential backends label a group contrast `up` / `down` and a
# continuous one `positive` / `negative`. The concordance quadrants are
# about sign, so both vocabularies mean the same thing here; anything
# else (`ns`, a zero effect, NA) has no sign.
direction_sign <- function(direction, effect = NULL) {
  d <- tolower(as.character(direction))
  out <- rep(NA_character_, length(d))
  out[d %in% c("up", "positive")] <- "up"
  out[d %in% c("down", "negative")] <- "down"
  if (!is.null(effect)) {
    miss <- is.na(out) & !is.na(effect) & !(d %in% "ns")
    out[miss & effect > 0] <- "up"
    out[miss & effect < 0] <- "down"
  }
  out
}

# Build a feature mapping between two experiments. Returns a data.frame
# with columns `feature_a`, `feature_b`, `feature_symbol`, `feature_id`.
# Uses `feature_symbol` first (since cross-omics integration is typically
# gene-symbol space), falling back to `feature_id`.
build_feature_pairs <- function(project, tag_a, tag_b, by = "feature_symbol") {
  input_a <- project$experiments[[tag_a]]
  input_b <- project$experiments[[tag_b]]
  feat_a <- input_a$feature_df
  feat_b <- input_b$feature_df
  if (!by %in% colnames(feat_a) || !by %in% colnames(feat_b)) {
    stop("`", by, "` must be a column in both experiments' `feature_df`.")
  }
  a <- data.frame(
    feature_a = feat_a$feature_id,
    key = feat_a[[by]],
    stringsAsFactors = FALSE
  )
  b <- data.frame(
    feature_b = feat_b$feature_id,
    key = feat_b[[by]],
    stringsAsFactors = FALSE
  )
  a$symbol <- a$key
  a$key <- integration_join_key(a$key)
  b$key <- integration_join_key(b$key)
  abund <- function(input, ids) {
    m <- input$expr_mat
    if (is.null(m) || !length(ids)) return(NULL)
    suppressWarnings(rowMeans(m[intersect(ids, rownames(m)), , drop = FALSE],
                              na.rm = TRUE))[ids]
  }
  a <- a[!is.na(a$key), , drop = FALSE]
  b <- b[!is.na(b$key), , drop = FALSE]
  a <- a[dedupe_by_key(a$key, abund(input_a, a$feature_a)), , drop = FALSE]
  b <- b[dedupe_by_key(b$key, abund(input_b, b$feature_b)), , drop = FALSE]
  pairs <- merge(a, b, by = "key")
  if (nrow(pairs) == 0L) {
    stop("No shared `", by, "` features between '", tag_a, "' and '", tag_b, "'.")
  }
  data.frame(
    feature_id = pairs$symbol,
    feature_symbol = pairs$symbol,
    feature_a = pairs$feature_a,
    feature_b = pairs$feature_b,
    stringsAsFactors = FALSE
  )
}

# Validate and coerce a `diff_bundles` argument: must be a named list keyed
# by experiment tags, each being a run_diff analysis_bundle.
validate_diff_bundles <- function(diff_bundles, experiments) {
  if (!is.list(diff_bundles) || is.null(names(diff_bundles))) {
    stop("`diff_bundles` must be a named list of run_diff bundles keyed by experiment tag.")
  }
  missing <- setdiff(experiments, names(diff_bundles))
  if (length(missing) > 0L) {
    stop("Missing diff bundles for experiments: ", paste(missing, collapse = ", "))
  }
  for (tag in experiments) {
    b <- diff_bundles[[tag]]
    if (!is_analysis_bundle(b) || !identical(b$analysis_name, "run_diff")) {
      stop("`diff_bundles$", tag, "` must be an analysis_bundle from run_diff().")
    }
  }
  invisible(TRUE)
}

# Classify a (direction_a, direction_b) pair into a 4-quadrant label.
classify_concordance_quadrant <- function(dir_a, dir_b) {
  out <- rep(NA_character_, length(dir_a))
  out[dir_a == "up"   & dir_b == "up"]   <- "up_up"
  out[dir_a == "down" & dir_b == "down"] <- "down_down"
  out[dir_a == "up"   & dir_b == "down"] <- "up_down"
  out[dir_a == "down" & dir_b == "up"]   <- "down_up"
  out
}
