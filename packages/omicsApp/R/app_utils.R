
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


#' Everything a result will be stamped with, captured before it is computed
#'
#' The whole record, at launch. The first attempt took the input versions
#' from launch and the project identity from arrival, which is worse than
#' taking both from arrival: it produced a record that looks internally
#' consistent -- versions that match, a project that exists -- while
#' describing a computation that happened somewhere else. A late result
#' from project A landed in project B, on identical data, and passed the
#' ownership check because the identity had been rewritten to B's on the
#' way in.
#'
#' @param project The live `omics_project` at launch.
#' @param view Which view is running.
#' @param tags Layer tags the run will read.
#' @param upstream Analysis ids it derives from.
#' @param params The resolved parameters.
#' @param stage `"effective"`, or `"baseline"` for QC.
#'
#' @return A list carrying `project_id`, `tags`, `versions`, `upstream`,
#'   `params_signature`, `stage`, and the `activation` it was launched under.
#' @keywords internal
#' @noRd
capture_launch <- function(project, view, tags, upstream = character(0),
                           params = list(), stage = "effective",
                           activation = NA_integer_) {
  if (is.null(project)) return(NULL)
  resolve <- if (identical(stage, "baseline")) {
    omicsCore::baseline_input
  } else {
    omicsCore::effective_input
  }
  versions <- vapply(tags, function(tg) resolve(project, tg)$version,
                     character(1), USE.NAMES = FALSE)
  list(project_id = project$metadata$project_id,
       view = view, tags = tags, versions = versions,
       upstream = as.character(upstream %||% character(0)),
       params_signature = omicsCore::params_signature(params),
       stage = stage, activation = activation)
}

#' Stamp a result from what was captured at launch
#'
#' Nothing here is read from the live application. If the project has moved
#' on, that is for the admission check to refuse -- not for this to paper
#' over by describing the result in terms of wherever it landed.
#'
#' @param bundle The `analysis_bundle` to stamp.
#' @param launch A `capture_launch()` record.
#'
#' @return The bundle, carrying a `provenance_binding`.
#' @keywords internal
#' @noRd
stamp_from_launch <- function(bundle, launch) {
  if (is.null(bundle) || is.null(launch)) return(bundle)
  inputs <- stats::setNames(as.list(launch$versions), launch$tags)
  bundle$input_info$layer_tag <- launch$tags
  bundle$input_info$input_version <- launch$versions
  bundle$input_info$upstream_ids <- launch$upstream
  bundle$input_info$input_stage <- launch$stage
  bundle$provenance_binding <- omicsCore::new_provenance_binding(
    # Derived, so re-filing the same result is idempotent -- the observer
    # that files results re-runs on every flush while the view goes on
    # holding the same bundle.
    analysis_id = paste0("ana_", omicsCore::params_signature(list(
      view = launch$view, inputs = inputs, upstream = launch$upstream,
      params = bundle$params, project = launch$project_id))),
    project_id = launch$project_id,
    inputs = inputs,
    upstream_ids = launch$upstream,
    params_signature = omicsCore::params_signature(bundle$params),
    input_stage = launch$stage
  )
  bundle
}

#' The gate a finished task has to pass before it may write
#'
#' Both halves of the contract, on both outcomes: the session token says the
#' user has not moved on, the captured sources say the data has not. Either
#' alone lets a stale result land, and a request number alone -- which is
#' all this used to check -- sees neither a change of project nor a change
#' of data.
#'
#' A context it cannot evaluate is a refusal. The earlier version defaulted
#' to admitting when it could not read the current request number, which is
#' the wrong way round for a check whose whole purpose is to catch the case
#' where things have moved.
#'
#' @param launch A `capture_launch()` record.
#' @param ctx A function returning the live
#'   `list(project, activation, request_seq)`.
#' @param request_seq The request number this task was launched with.
#'
#' @return A function suitable for `run_async(admissible = )`.
#' @keywords internal
#' @noRd
launch_admissible <- function(launch, ctx, request_seq) {
  function(outcome) {
    now <- tryCatch(ctx(), error = function(e) NULL)
    if (is.null(now)) {
      return(list(ok = FALSE, reason = "the session context cannot be read"))
    }
    # No project: the built-in demo, which runs the real analyses on
    # fixtures. There is no ownership question because there is nowhere to
    # attribute the result -- and stamp_from_launch() leaves it unstamped,
    # so the write-back observer will not file it either. What still
    # applies is the request number: pressing Run twice on the demo should
    # not let the first answer win.
    if (is.null(launch)) {
      if (!identical(now$request_seq, request_seq)) {
        return(list(ok = FALSE, reason = "a newer run was started"))
      }
      return(list(ok = TRUE, reason = ""))
    }
    if (is.null(now$project)) {
      return(list(ok = FALSE, reason = "the project is no longer open"))
    }
    if (!identical(now$activation, launch$activation)) {
      return(list(ok = FALSE, reason = "the project was reopened"))
    }
    if (!identical(now$request_seq, request_seq)) {
      return(list(ok = FALSE, reason = "a newer run was started"))
    }
    if (!identical(now$project$metadata$project_id, launch$project_id)) {
      return(list(ok = FALSE, reason = "the project changed"))
    }
    resolve <- if (identical(launch$stage, "baseline")) {
      omicsCore::baseline_input
    } else {
      omicsCore::effective_input
    }
    for (i in seq_along(launch$tags)) {
      tg <- launch$tags[[i]]
      if (!tg %in% names(now$project$experiments)) {
        return(list(ok = FALSE, reason = sprintf("layer '%s' is gone", tg)))
      }
      if (!identical(resolve(now$project, tg)$version, launch$versions[[i]])) {
        return(list(ok = FALSE, reason = "the data changed while this ran"))
      }
    }
    if (length(launch$upstream) > 0L) {
      alive <- omicsCore::valid_result_ids(now$project)
      if (!all(launch$upstream %in% alive)) {
        return(list(ok = FALSE, reason = "an upstream result is no longer valid"))
      }
    }
    list(ok = TRUE, reason = "")
  }
}
