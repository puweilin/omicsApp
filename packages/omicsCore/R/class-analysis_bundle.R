#' Create an analysis bundle
#'
#' An `analysis_bundle` is the standard return container for every analysis
#' function in `omicsCore` (QC, differential, enrichment, integration). It
#' bundles results, parameter provenance, artifact registry, and informational
#' messages so downstream code (and the Shiny UI) can consume them uniformly.
#'
#' @param analysis_name Bundle name (e.g. `"run_diff"`).
#' @param input_info Named list describing inputs (omics_type, sample counts,
#'   etc.) for provenance.
#' @param params Named list of analysis parameters.
#' @param results Named list of result objects (data frames, matrices, plots).
#' @param artifacts Artifact registry; defaults to an empty registry.
#' @param messages Character vector of informational messages.
#' @param warnings Character vector of warning messages.
#'
#' @return An object of class `analysis_bundle`.
#' @keywords internal
new_analysis_bundle <- function(
  analysis_name,
  input_info = list(),
  params = list(),
  results = list(),
  artifacts = NULL,
  messages = character(0),
  warnings = character(0)
) {
  if (is.null(artifacts)) {
    artifacts <- new_artifact_registry()
  }

  structure(
    list(
      analysis_name = analysis_name,
      input_info = input_info,
      params = params,
      results = results,
      artifacts = artifacts,
      messages = messages,
      warnings = warnings
    ),
    class = "analysis_bundle"
  )
}

#' Test whether an object is an analysis bundle
#'
#' @param x Object to test.
#'
#' @return Logical scalar.
#' @export
is_analysis_bundle <- function(x) {
  inherits(x, "analysis_bundle")
}

#' Print method for `analysis_bundle`
#'
#' @param x An `analysis_bundle` object.
#' @param ... Unused.
#'
#' @return Invisibly returns `x`.
#' @export
print.analysis_bundle <- function(x, ...) {
  cat("<analysis_bundle>\n")
  cat("  analysis :", x$analysis_name, "\n")
  cat("  results  :", paste(names(x$results), collapse = ", "), "\n")
  cat("  params   :", length(x$params), "entries\n")
  cat("  artifacts:", nrow(x$artifacts), "rows\n")
  if (length(x$messages) > 0) {
    cat("  messages :", length(x$messages), "\n")
  }
  if (length(x$warnings) > 0) {
    cat("  warnings :", length(x$warnings), "\n")
  }
  invisible(x)
}
