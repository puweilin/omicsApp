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
# `data.frame` with columns `<tag_a>`, `<tag_b>`, `donor_id`. If a
# project-level `sample_link` is present it is used; otherwise direct
# matches on `sample_id` are emitted.
build_sample_pairs <- function(project, tag_a, tag_b) {
  input_a <- project$experiments[[tag_a]]
  input_b <- project$experiments[[tag_b]]
  samp_a <- colnames(input_a$expr_mat)
  samp_b <- colnames(input_b$expr_mat)

  if (!is.null(project$sample_link) && nrow(project$sample_link) > 0L) {
    sl <- project$sample_link
    a <- sl[sl$tag == tag_a, c("sample_id", "donor_id"), drop = FALSE]
    b <- sl[sl$tag == tag_b, c("sample_id", "donor_id"), drop = FALSE]
    pairs <- merge(a, b, by = "donor_id", suffixes = c("_a", "_b"))
    pairs <- pairs[pairs$sample_id_a %in% samp_a &
                     pairs$sample_id_b %in% samp_b, , drop = FALSE]
    out <- data.frame(
      donor_id = pairs$donor_id,
      a = pairs$sample_id_a,
      b = pairs$sample_id_b,
      stringsAsFactors = FALSE
    )
  } else {
    shared <- intersect(samp_a, samp_b)
    if (length(shared) == 0L) {
      stop("No shared sample IDs between '", tag_a, "' and '", tag_b,
           "', and no `sample_link` is set on the project.")
    }
    out <- data.frame(
      donor_id = shared,
      a = shared,
      b = shared,
      stringsAsFactors = FALSE
    )
  }
  names(out)[2:3] <- c(tag_a, tag_b)
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
  a <- a[!is.na(a$key) & nzchar(a$key), , drop = FALSE]
  b <- b[!is.na(b$key) & nzchar(b$key), , drop = FALSE]
  a <- a[!duplicated(a$key), , drop = FALSE]
  b <- b[!duplicated(b$key), , drop = FALSE]
  pairs <- merge(a, b, by = "key")
  if (nrow(pairs) == 0L) {
    stop("No shared `", by, "` features between '", tag_a, "' and '", tag_b, "'.")
  }
  data.frame(
    feature_id = pairs$key,
    feature_symbol = pairs$key,
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
