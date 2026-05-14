#' Multi-omics project container
#'
#' An `omics_project` groups multiple [omics_input()] experiments together
#' with optional cross-omics mapping tables, so that integration analyses
#' (RNA × Protein, ActivePathways, etc.) can operate on a single root object.
#'
#' Each entry in `experiments` is one omics layer (e.g. `"proteomics"`,
#' `"rnaseq"`), indexed by an arbitrary tag. `sample_link` maps sample IDs
#' across layers so paired samples can be aligned; `feature_link` maps
#' feature IDs (e.g. UniProt accession ↔ gene symbol) so cross-omics
#' integration can join on shared identifiers.
#'
#' @param name Human-readable project label.
#' @param experiments Named list of [omics_input()] objects. Names become
#'   experiment tags. May be empty; use [add_experiment()] to populate.
#' @param sample_link Optional `data.frame` with at least the columns
#'   `tag`, `sample_id`, `donor_id`. `tag` matches a name in `experiments`.
#'   Each row links a per-experiment sample to a shared donor.
#' @param feature_link Optional `data.frame` mapping feature IDs across
#'   omics (e.g. columns `uniprot`, `gene_symbol`).
#' @param metadata Optional named list of arbitrary project-level metadata.
#'
#' @return An object of class `omics_project`.
#' @export
#' @family omics_project
#' @examples
#' p <- omics_project(name = "demo")
#' p
omics_project <- function(
  name,
  experiments = list(),
  sample_link = NULL,
  feature_link = NULL,
  metadata = list()
) {
  if (!is.character(name) || length(name) != 1L || is.na(name)) {
    stop("`name` must be a single non-NA string.")
  }
  if (!is.list(experiments)) {
    stop("`experiments` must be a (possibly empty) named list of omics_input objects.")
  }
  if (length(experiments) > 0L) {
    if (is.null(names(experiments)) || any(!nzchar(names(experiments)))) {
      stop("All entries in `experiments` must be named.")
    }
    if (anyDuplicated(names(experiments))) {
      stop("Experiment tags must be unique.")
    }
    not_input <- !vapply(experiments, is_omics_input, logical(1))
    if (any(not_input)) {
      stop(
        "All `experiments` entries must be `omics_input` objects. ",
        "Offending tags: ", paste(names(experiments)[not_input], collapse = ", ")
      )
    }
  }

  if (!is.null(sample_link)) {
    validate_sample_link(sample_link, tags = names(experiments))
  }
  if (!is.null(feature_link)) {
    if (!is.data.frame(feature_link)) {
      stop("`feature_link` must be a data.frame or NULL.")
    }
  }

  structure(
    list(
      name = name,
      experiments = experiments,
      sample_link = sample_link,
      feature_link = feature_link,
      metadata = metadata,
      created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      schema_version = "1"
    ),
    class = "omics_project"
  )
}

#' Test whether an object is an `omics_project`
#'
#' @param x Object to test.
#'
#' @return Logical scalar.
#' @export
#' @family omics_project
is_omics_project <- function(x) {
  inherits(x, "omics_project")
}

#' Add an experiment to a project
#'
#' @param project An `omics_project`.
#' @param name Experiment tag (e.g. `"proteomics"`, `"rnaseq"`). Must be
#'   unique within the project.
#' @param input An `omics_input` to attach.
#'
#' @return The modified `omics_project`.
#' @export
#' @family omics_project
add_experiment <- function(project, name, input) {
  if (!is_omics_project(project)) {
    stop("`project` must be an `omics_project`.")
  }
  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
    stop("`name` must be a single non-empty string.")
  }
  if (name %in% names(project$experiments)) {
    stop("Experiment tag already exists in project: ", name)
  }
  if (!is_omics_input(input)) {
    stop("`input` must be an `omics_input`.")
  }
  validate_omics_input(input)

  project$experiments[[name]] <- input
  project
}

#' Remove an experiment from a project
#'
#' @param project An `omics_project`.
#' @param name Experiment tag to remove.
#'
#' @return The modified `omics_project`.
#' @export
#' @family omics_project
remove_experiment <- function(project, name) {
  if (!is_omics_project(project)) {
    stop("`project` must be an `omics_project`.")
  }
  if (!name %in% names(project$experiments)) {
    stop("Experiment tag not found in project: ", name)
  }
  project$experiments[[name]] <- NULL
  project
}

#' List experiment tags in a project
#'
#' @param project An `omics_project`.
#'
#' @return Character vector of experiment tags.
#' @export
#' @family omics_project
experiment_tags <- function(project) {
  if (!is_omics_project(project)) {
    stop("`project` must be an `omics_project`.")
  }
  names(project$experiments) %||% character(0)
}

#' Print method for `omics_project`
#'
#' @param x An `omics_project` object.
#' @param ... Unused.
#'
#' @return Invisibly returns `x`.
#' @export
print.omics_project <- function(x, ...) {
  cat("<omics_project>\n")
  cat("  name        :", x$name, "\n")
  cat("  experiments :", length(x$experiments), "\n")
  for (tag in names(x$experiments)) {
    e <- x$experiments[[tag]]
    cat(sprintf("    - %s  (%s, %d features x %d samples)\n",
                tag, e$omics_type %||% "<unset>",
                nrow(e$expr_mat), ncol(e$expr_mat)))
  }
  cat("  sample_link :", if (is.null(x$sample_link)) "<none>" else paste0(nrow(x$sample_link), " rows"), "\n")
  cat("  feature_link:", if (is.null(x$feature_link)) "<none>" else paste0(nrow(x$feature_link), " rows"), "\n")
  invisible(x)
}

# Internal: validate a sample_link frame against project experiment tags.
validate_sample_link <- function(sample_link, tags) {
  if (!is.data.frame(sample_link)) {
    stop("`sample_link` must be a data.frame or NULL.")
  }
  required <- c("tag", "sample_id", "donor_id")
  check_required_cols(sample_link, required, object_name = "sample_link")
  unknown_tags <- setdiff(unique(sample_link$tag), tags)
  if (length(tags) > 0L && length(unknown_tags) > 0L) {
    stop(
      "`sample_link$tag` references experiments not in the project: ",
      paste(unknown_tags, collapse = ", ")
    )
  }
  invisible(TRUE)
}
