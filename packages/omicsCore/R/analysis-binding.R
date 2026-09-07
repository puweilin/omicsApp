# Which data a result is about, and whether a running task may still write.
#
# Those are two questions, and folding them into one field is what made the
# first attempt at this contract self-defeating: a session counter that
# increments on every open was written into the record that gets saved to
# disk, so every restored project judged its own saved analyses invalid.
#
# So they are kept apart:
#
#   provenance_binding   persisted. Which project, which input version per
#                        layer, which upstream analyses, which parameters.
#                        Answers "is this saved result still about this
#                        project's current data".
#
#   session_token        not persisted. Session, activation, per-view request
#                        number, parameter signature. Answers "may this
#                        in-flight callback still write".
#
# A write has to satisfy both. The token alone says the user has not moved
# on; the binding alone says the data has not moved on; the data can change
# under a task without the user starting a new one.

# ---- identity ---------------------------------------------------------

# Deliberately not sample(): drawing from the analysis RNG stream both
# advances a sequence the analysis is about to use and repeats itself under
# a fixed seed, so two results in one seeded script would share an id.
new_opaque_id <- function(prefix) {
  paste0(prefix, "_",
         format(as.numeric(Sys.time()) * 1000, scientific = FALSE, digits = 15),
         "_", basename(tempfile("")))
}

#' Generate an analysis identifier
#' @keywords internal
#' @noRd
new_analysis_id <- function() new_opaque_id("ana")

#' Generate a project identifier
#' @keywords internal
#' @noRd
new_project_id <- function() new_opaque_id("prj")

#' Content version of an analysis input
#'
#' Everything a backend can see. `feature_df` is in because `feature_symbol`
#' is the enrichment gene mapping and the integration join key -- editing it
#' changes the analysis without touching a single measurement. Shape is not
#' enough on its own either: two unrelated 20x10 matrices have the same one.
#'
#' @param input An `omics_input`.
#' @return A single string.
#' @keywords internal
#' @noRd
input_version <- function(input) {
  rlang::hash(list(
    matrix     = input$expr_mat,
    samples    = colnames(input$expr_mat),
    features   = rownames(input$expr_mat),
    meta       = input$meta_df,
    feature_df = input$feature_df,
    omics_type = input$omics_type,
    assay_type = input$assay_type,
    provenance = input$provenance
  ))
}

#' Stable signature of a resolved parameter set
#'
#' Resolved, not as typed: `run_qc(impute_method = NULL)` and
#' `run_qc(impute_method = "MinProb")` are the same analysis on proteomics.
#'
#' @keywords internal
#' @noRd
params_signature <- function(params) {
  if (is.null(params)) params <- list()
  params <- params[order(names(params))]
  rlang::hash(params)
}

# ---- persisted provenance ---------------------------------------------

#' Record of what a result was computed from
#'
#' @param analysis_id Immutable id of this result.
#' @param project_id Project the result belongs to.
#' @param inputs Named list, tag -> input version. Multi-layer analyses name
#'   every layer they read; one entry per analysis is not assumed.
#' @param upstream_ids Analysis ids this result was derived from.
#' @param params_signature Signature of the resolved parameters.
#'
#' @keywords internal
#' @noRd
new_provenance_binding <- function(analysis_id, project_id, inputs,
                                   upstream_ids = character(0),
                                   params_signature = NA_character_) {
  list(analysis_id = analysis_id, project_id = project_id,
       inputs = inputs, upstream_ids = as.character(upstream_ids),
       params_signature = params_signature)
}

#' Check that a binding is well formed before believing what it says
#'
#' Structure first, currency second. `all(logical(0))` is `TRUE`, so a
#' binding that names no inputs answered "all my inputs are current" -- which
#' is how an empty record passed for a verified one.
#'
#' @return `list(ok = , reason = )`.
#' @keywords internal
#' @noRd
validate_provenance_binding <- function(pb) {
  bad <- function(reason) list(ok = FALSE, reason = reason)
  if (!is.list(pb)) return(bad("binding is not a list"))

  chr1 <- function(x) is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
  if (!chr1(pb$analysis_id)) return(bad("analysis_id is missing or blank"))
  if (!chr1(pb$project_id))  return(bad("project_id is missing or blank"))

  if (!is.list(pb$inputs) || length(pb$inputs) == 0L) {
    return(bad("binding names no inputs"))
  }
  nms <- names(pb$inputs)
  if (is.null(nms) || any(!nzchar(nms)) || anyDuplicated(nms) > 0L) {
    return(bad("input versions are not uniquely named by layer"))
  }
  if (!all(vapply(pb$inputs, chr1, logical(1)))) {
    return(bad("an input version is missing or blank"))
  }
  if (!is.character(pb$upstream_ids)) {
    return(bad("upstream_ids is not a character vector"))
  }
  if (pb$analysis_id %in% pb$upstream_ids) {
    return(bad("binding depends on itself"))
  }
  if (!chr1(pb$params_signature)) return(bad("params_signature is missing"))
  list(ok = TRUE, reason = "")
}

#' Are the inputs this result names still the project's current ones?
#'
#' Ownership is part of the question. Two projects built on the same file
#' have identical input versions, so comparing versions alone accepted one
#' project's results inside the other.
#'
#' @keywords internal
#' @noRd
result_inputs_current <- function(pb, project) {
  if (!isTRUE(validate_provenance_binding(pb)$ok)) return(FALSE)
  if (!identical(pb$project_id, project$metadata$project_id)) return(FALSE)
  tags <- names(pb$inputs)
  if (!all(tags %in% names(project$experiments))) return(FALSE)
  all(vapply(tags, function(tg) {
    identical(pb$inputs[[tg]], effective_input(project, tg)$version)
  }, logical(1)))
}

# ---- the analysis registry --------------------------------------------

#' Record an analysis under its own id
#'
#' The view slots in `project$bundles` say what is on screen; this says what
#' exists. They are different questions, and conflating them meant an
#' analysis that nothing displayed -- the second layer's differential inside
#' an integration -- could not be found by the result that depended on it.
#'
#' @keywords internal
#' @noRd
register_analysis <- function(project, bundle, id = NULL) {
  id <- id %||% bundle$provenance_binding$analysis_id %||% new_analysis_id()
  if (is.null(project$analyses)) project$analyses <- list()
  project$analyses[[id]] <- bundle
  project
}

#' The analyses that are still about this project's current data
#'
#' Iterated to a fixed point. Validity is recursive -- an integration is only
#' as valid as the two differentials under it -- and a single pass keeps
#' whichever dependant happened to be visited before its upstream was
#' removed, so the answer depended on list order.
#'
#' @keywords internal
#' @noRd
valid_result_ids <- function(project) {
  analyses <- project$analyses %||% list()
  ids <- names(analyses)
  if (length(ids) == 0L) return(character(0))
  repeat {
    keep <- Filter(function(id) {
      pb <- analyses[[id]]$provenance_binding
      # No binding: an .omp written before the contract existed. It loads,
      # it just cannot claim to be about the current data on its own.
      if (is.null(pb)) return(FALSE)
      result_inputs_current(pb, project) && all(pb$upstream_ids %in% ids)
    }, ids)
    if (identical(keep, ids)) return(ids)
    ids <- keep
    if (length(ids) == 0L) return(character(0))
  }
}

# Fields worth keeping when a result becomes history, versus the fitted
# objects that make an autosave grow without bound. A model lives under
# several names depending on the backend, so this drops by allow-list.
HISTORY_RESULT_FIELDS <- c("diff_result_df", "enrich_result_df",
                           "integration_df", "qc_summary")

compress_for_history <- function(bundle, reason) {
  keep <- intersect(names(bundle$results %||% list()), HISTORY_RESULT_FIELDS)
  list(
    analysis_name = bundle$analysis_name,
    provenance_binding = bundle$provenance_binding,
    params = bundle$params,
    input_info = bundle$input_info,
    results = bundle$results[keep],
    retired_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    retired_because = reason
  )
}

#' Move results that are no longer current into history
#'
#' `max_history` bounds what is kept for the user to look at. It deliberately
#' does not touch `project$analyses`: the cap is a storage and display limit,
#' and applying it to the dependency graph would delete the upstream node a
#' still-valid result points at.
#'
#' @keywords internal
#' @noRd
retire_invalid <- function(project, max_history = 20L,
                           reason = "inputs changed") {
  analyses <- project$analyses %||% list()
  if (length(analyses) == 0L) return(project)
  keep <- valid_result_ids(project)
  gone <- setdiff(names(analyses), keep)
  if (length(gone) > 0L) {
    history <- lapply(analyses[gone], compress_for_history, reason = reason)
    project$retired <- c(project$retired %||% list(), history)
    project$analyses <- analyses[keep]
    # Views hold results by slot name, not by id; drop any that just left.
    if (!is.null(project$bundles)) {
      still <- vapply(project$bundles, function(b) {
        (b$provenance_binding$analysis_id %||% NA_character_) %in% keep
      }, logical(1))
      project$bundles <- project$bundles[still]
    }
  }
  if (length(project$retired %||% list()) > max_history) {
    n <- length(project$retired)
    project$retired <- project$retired[seq.int(n - max_history + 1L, n)]
  }
  project
}

# ---- per-layer preprocessing ------------------------------------------

#' Version of a layer, computed once and carried with it
#'
#' Cached on the object rather than recomputed per query: `effective_input()`
#' is read from reactive expressions, and hashing a 63k x 258 matrix on every
#' read is not a cache.
#'
#' @keywords internal
#' @noRd
stamp_version <- function(input) {
  input$version <- input_version(input)
  input
}

#' Record the preprocessing applied to one layer
#'
#' Per layer, because a project can hold several and each is filtered and
#' imputed on its own terms. A single slot let the second layer's QC replace
#' the first layer's.
#'
#' `base_input_version` is what makes this safe to serve: the record proves
#' which import it was built from, so a replaced layer cannot silently keep
#' handing out preprocessing derived from the layer before it.
#'
#' @keywords internal
#' @noRd
set_preprocessing <- function(project, tag, cleaned, params,
                              base_input_version = NULL) {
  raw <- project$experiments[[tag]]
  if (is.null(raw)) stop("No such experiment layer: ", tag, call. = FALSE)
  if (is.null(raw$version)) {
    project$experiments[[tag]] <- stamp_version(raw)
    raw <- project$experiments[[tag]]
  }
  if (is.null(project$preprocessing)) project$preprocessing <- list()
  project$preprocessing[[tag]] <- list(
    base_input_version = base_input_version %||% raw$version,
    input_version = input_version(cleaned),
    cleaned = cleaned,
    params = params
  )
  project
}

#' Drop a layer's preprocessing
#' @keywords internal
#' @noRd
clear_preprocessing <- function(project, tag) {
  project$preprocessing[[tag]] <- NULL
  project
}

#' The imported layer, before any preprocessing
#'
#' QC computes from here rather than from `effective_input()`. Re-running QC
#' over already-filtered values cannot return the features the previous run
#' removed, so loosening a threshold would appear to do nothing, and changing
#' the imputation method would impute values that are no longer missing.
#'
#' @keywords internal
#' @noRd
baseline_input <- function(project, tag) {
  raw <- project$experiments[[tag]]
  if (is.null(raw)) stop("No such experiment layer: ", tag, call. = FALSE)
  if (is.null(raw$version)) raw <- stamp_version(raw)
  list(input = raw, version = raw$version, provenance = raw$provenance)
}

#' The input every analysis of this layer must use
#'
#' The single reader for "which data does this layer mean right now". Modules
#' that went to `project$experiments[[tag]]` directly are how a project ended
#' up with one layer analysed after QC and another before it.
#'
#' @return `list(input, version, source, provenance, stale)`. `stale` marks
#'   preprocessing that was left behind by a layer replacement -- the caller
#'   gets the baseline and can say so, rather than being handed cleaned data
#'   derived from an import that is no longer here.
#' @keywords internal
#' @noRd
effective_input <- function(project, tag) {
  base <- baseline_input(project, tag)
  pp <- project$preprocessing[[tag]]
  if (!is.null(pp) && identical(pp$base_input_version, base$version)) {
    return(list(input = pp$cleaned, version = pp$input_version,
                source = "qc_cleaned",
                provenance = pp$cleaned$provenance, stale = FALSE))
  }
  list(input = base$input, version = base$version, source = "imported",
       provenance = base$provenance, stale = !is.null(pp))
}

# ---- session admission ------------------------------------------------

#' Token stamped on a task when it is launched
#'
#' Not persisted. `activation` counts opens within a session, so writing it
#' into a saved result would make that result invalid the next time the file
#' is opened.
#'
#' @keywords internal
#' @noRd
new_session_token <- function(session_id, project_id, activation, request_seq,
                              view, params_signature) {
  list(session_id = session_id, project_id = project_id,
       activation = activation, request_seq = request_seq, view = view,
       params_signature = params_signature)
}

#' May this callback write its result?
#'
#' Both halves, on both outcomes. A failed task writing its error into a view
#' the user has already moved on from is the same defect as a stale success,
#' so `outcome` is recorded rather than used to skip the check.
#'
#' Missing fields are refusals. `identical(NULL, NULL)` is `TRUE`, so a token
#' built without a parameter signature satisfied the comparison meant to
#' catch it.
#'
#' @param token The token captured at launch.
#' @param ctx Current session state: `session_id`, `project_id`,
#'   `activation`, `request_seq` (named by view).
#' @param pb The provenance binding the task would write.
#' @param project The project as it stands now.
#' @param outcome `"success"` or `"error"`.
#'
#' @return `list(ok = , reason = )`.
#' @keywords internal
#' @noRd
callback_admissible <- function(token, ctx, pb, project,
                                outcome = c("success", "error")) {
  outcome <- match.arg(outcome)
  bad <- function(reason) list(ok = FALSE, reason = reason, outcome = outcome)

  if (!is.list(token) || !is.list(ctx)) return(bad("no session context"))
  if (is.null(token$params_signature)) return(bad("token carries no params signature"))
  if (is.null(token$view) || !nzchar(token$view)) return(bad("token names no view"))

  if (!identical(token$session_id, ctx$session_id)) return(bad("session ended"))
  if (!identical(token$project_id, ctx$project_id)) return(bad("project changed"))
  if (!identical(token$activation, ctx$activation)) return(bad("project reopened"))

  current_seq <- ctx$request_seq[[token$view]]
  if (is.null(current_seq)) return(bad("no request outstanding for this view"))
  if (!identical(token$request_seq, current_seq)) return(bad("a newer request was made"))

  if (!identical(token$params_signature, pb$params_signature)) {
    return(bad("params signature does not match the task"))
  }

  structure_ok <- validate_provenance_binding(pb)
  if (!isTRUE(structure_ok$ok)) return(bad(structure_ok$reason))
  if (!result_inputs_current(pb, project)) return(bad("input versions moved on"))
  if (!all(pb$upstream_ids %in% valid_result_ids(project))) {
    return(bad("an upstream result is no longer valid"))
  }
  list(ok = TRUE, reason = "", outcome = outcome)
}
