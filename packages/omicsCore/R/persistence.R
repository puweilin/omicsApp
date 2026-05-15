# Atomic save / load of an omics_project. We use qs2 because it round-trips
# arbitrary R objects (DESeqDataSet, MArrayLM, ComplexHeatmap, ggplot, etc.)
# losslessly and is faster than saveRDS on the kinds of objects analysis
# bundles tend to hold. qs2 is Suggests-gated so a base install of
# omicsCore stays slim.
#
# File format: a qs2 archive whose root is a list with two fields:
#   * `schema_version`  -- semantic-version-ish string for forward compat
#   * `payload`         -- the `omics_project`
# Future schema changes should bump the major component and add a
# migration in `load_project()`.

OMP_SCHEMA_VERSION <- "1.0.0"

ensure_qs2 <- function() {
  if (!requireNamespace("qs2", quietly = TRUE)) {
    stop(
      "Package 'qs2' is required for save_project() / load_project(). ",
      "Install with: install_optional('persistence').",
      call. = FALSE
    )
  }
}

#' Save an omics_project to disk
#'
#' Writes the full project (every experiment, sample_link, and any
#' attached `analysis_bundle`s) to a single `.omp` file using `qs2`. The
#' write is atomic: the payload first goes to `path.tmp`, then renames
#' over `path`, so an interrupted save will not corrupt an existing file.
#'
#' Requires the `qs2` package; install it with `install_optional("persistence")`.
#'
#' @param project An [`omics_project`][is_omics_project()].
#' @param path Output path. The file extension is up to the caller, but
#'   the convention is `.omp` for projects produced by `omicsCore`.
#' @param overwrite If `FALSE` (default) and `path` already exists,
#'   `save_project()` raises an error. Set to `TRUE` to replace.
#'
#' @return Invisibly returns `path`.
#' @export
#' @family persistence
#' @examples
#' \dontrun{
#'   p <- omics_project("demo", experiments = list(proteo = my_input))
#'   save_project(p, "demo.omp")
#'   q <- load_project("demo.omp")
#' }
save_project <- function(project, path, overwrite = FALSE) {
  if (!is_omics_project(project)) {
    stop("`project` must be an `omics_project`.")
  }
  if (!is.character(path) || length(path) != 1L || !nzchar(path)) {
    stop("`path` must be a non-empty single string.")
  }
  if (file.exists(path) && !isTRUE(overwrite)) {
    stop("Path already exists: ", path,
         " (pass `overwrite = TRUE` to replace).")
  }
  ensure_qs2()

  payload <- list(
    schema_version = OMP_SCHEMA_VERSION,
    payload = project
  )

  tmp <- paste0(path, ".tmp")
  on.exit(if (file.exists(tmp)) file.remove(tmp), add = TRUE)
  qs2::qs_save(payload, file = tmp)
  ok <- file.rename(tmp, path)
  if (!isTRUE(ok)) {
    stop("Failed to rename '", tmp, "' to '", path, "'.")
  }
  invisible(path)
}

#' Load an omics_project from disk
#'
#' Reads an `.omp` archive written by [save_project()]. Validates the
#' schema version and that the payload is an `omics_project`. Future
#' schema changes will add migrations here.
#'
#' Requires the `qs2` package; install it with `install_optional("persistence")`.
#'
#' @param path Path to the `.omp` file.
#'
#' @return The deserialized `omics_project`.
#' @export
#' @family persistence
load_project <- function(path) {
  if (!is.character(path) || length(path) != 1L) {
    stop("`path` must be a single string.")
  }
  if (!file.exists(path)) {
    stop("File does not exist: ", path)
  }
  ensure_qs2()
  payload <- qs2::qs_read(path)

  if (!is.list(payload) ||
      !all(c("schema_version", "payload") %in% names(payload))) {
    stop("`", path, "` does not appear to be an omicsCore project file.")
  }
  schema <- payload$schema_version
  if (!is.character(schema) || length(schema) != 1L) {
    stop("Invalid schema_version in '", path, "'.")
  }
  if (!schema_is_supported(schema)) {
    stop("Unsupported project schema version: '", schema,
         "'. This omicsCore supports up to '", OMP_SCHEMA_VERSION, "'.")
  }

  project <- payload$payload
  if (!is_omics_project(project)) {
    stop("Loaded payload is not an `omics_project`.")
  }
  project
}

# Accept anything with the same major version as the current schema.
schema_is_supported <- function(version) {
  current <- as.integer(strsplit(OMP_SCHEMA_VERSION, ".", fixed = TRUE)[[1L]][[1L]])
  parts <- strsplit(version, ".", fixed = TRUE)[[1L]]
  if (length(parts) == 0L) return(FALSE)
  major <- suppressWarnings(as.integer(parts[[1L]]))
  !is.na(major) && major <= current
}
