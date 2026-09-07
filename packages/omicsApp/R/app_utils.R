
#' Stamp a result with what it was computed from
#'
#' Called by the view that produced the result, at the moment it arrives,
#' with the layers and versions the view resolved when it launched the run.
#'
#' It belongs here rather than in the observer that files results into the
#' project, for two reasons. A downstream view reads its upstream's
#' `analysis_id` off the bundle it was handed, so a binding attached later --
#' to a copy -- is a binding the dependent never sees. And that observer
#' re-runs on every flush while the view goes on holding the same bundle, so
#' anything it minted there would be minted again each pass.
#'
#' The id is derived rather than drawn: the same analysis, on the same data,
#' with the same settings, is the same analysis. That makes re-filing it
#' idempotent, which is what keeps the observer from chasing its own writes.
#'
#' @param bundle The `analysis_bundle` to stamp.
#' @param project The live `omics_project`.
#' @param view Which view produced it.
#' @param tags Layer tags the result was computed on.
#' @param versions Input versions for those tags, in the same order.
#' @param upstream Analysis ids this result was derived from.
#' @param stage `"effective"`, or `"baseline"` for QC.
#'
#' @return The bundle, carrying a `provenance_binding`.
#' @keywords internal
#' @noRd
bind_result <- function(bundle, project, view, tags, versions,
                        upstream = character(0), stage = "effective") {
  if (is.null(project) || is.null(bundle)) return(bundle)
  if (length(tags) == 0L || length(tags) != length(versions)) return(bundle)
  inputs <- stats::setNames(as.list(as.character(versions)), tags)
  upstream <- as.character(upstream %||% character(0))
  pid <- project$metadata$project_id
  bundle$input_info$layer_tag <- tags
  bundle$input_info$input_version <- as.character(versions)
  bundle$input_info$upstream_ids <- upstream
  bundle$input_info$input_stage <- stage
  bundle$provenance_binding <- omicsCore::new_provenance_binding(
    analysis_id = paste0("ana_", omicsCore::params_signature(list(
      view = view, inputs = inputs, upstream = upstream,
      params = bundle$params, project = pid))),
    project_id = pid,
    inputs = inputs,
    upstream_ids = upstream,
    params_signature = omicsCore::params_signature(bundle$params),
    input_stage = stage
  )
  bundle
}
