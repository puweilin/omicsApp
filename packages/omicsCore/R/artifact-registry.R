#' Create an empty artifact registry
#'
#' An artifact registry is a `data.frame` tracking files written to disk by
#' analysis functions. Each row records one artifact's type, label, file
#' path, and creation timestamp.
#'
#' @return An empty registry `data.frame`.
#' @keywords internal
new_artifact_registry <- function() {
  data.frame(
    artifact_type = character(0),
    label = character(0),
    path = character(0),
    created_at = character(0),
    stringsAsFactors = FALSE
  )
}

#' Register an analysis artifact
#'
#' @param registry Existing registry or `NULL`.
#' @param artifact_type Artifact category such as `"plot"` or `"table"`.
#' @param label Logical artifact label.
#' @param path Optional file path.
#'
#' @return Updated registry.
#' @keywords internal
#' @importFrom dplyr bind_rows
register_artifact <- function(registry, artifact_type, label, path = NA_character_) {
  if (is.null(registry)) {
    registry <- new_artifact_registry()
  }

  entry <- data.frame(
    artifact_type = artifact_type,
    label = label,
    path = path,
    created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    stringsAsFactors = FALSE
  )

  dplyr::bind_rows(registry, entry)
}

#' Merge artifact registries
#'
#' @param ... Artifact registry objects.
#' @param deduplicate Whether to remove duplicated rows after merging.
#'
#' @return Merged artifact registry.
#' @keywords internal
#' @importFrom dplyr bind_rows distinct
merge_artifact_registries <- function(..., deduplicate = TRUE) {
  registries <- list(...)
  registries <- registries[!vapply(registries, is.null, logical(1))]

  if (length(registries) == 0) {
    return(new_artifact_registry())
  }

  merged <- dplyr::bind_rows(registries)

  if (isTRUE(deduplicate) && nrow(merged) > 0) {
    merged <- dplyr::distinct(merged)
  }

  merged
}
