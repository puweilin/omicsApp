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
#' @return A single opaque string, unique per call.
#' @keywords internal
#' @export
new_analysis_id <- function() new_opaque_id("ana")

#' Generate a project identifier
#' @return A single opaque string, unique per call.
#' @keywords internal
#' @export
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
#' @export
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
#' @param params A named list of resolved parameters.
#' @return A single string.
#' @keywords internal
#' @export
params_signature <- function(params) {
  if (is.null(params)) params <- list()
  # Sorted so that the same settings hash the same however they were
  # assembled -- but only when there is something to sort by. An analysis
  # with no parameters is ordinary (QC on defaults, an integration with a
  # fixed method), and `order(NULL)` is an error, not an empty ordering.
  nms <- names(params)
  if (length(params) > 0L && !is.null(nms)) {
    params <- params[order(nms)]
  }
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
#' @param input_stage Which resolver the recorded versions refer to.
#'   `"effective"` for analyses, which read whatever preprocessing is applied;
#'   `"baseline"` for QC, which reads the import. Without the distinction,
#'   applying QC would retire the QC result that produced the preprocessing,
#'   because the effective version had just stopped matching the one QC ran on.
#'
#' @return A list describing what the result was computed from.
#' @keywords internal
#' @export
new_provenance_binding <- function(analysis_id, project_id, inputs,
                                   upstream_ids = character(0),
                                   params_signature = NA_character_,
                                   input_stage = c("effective", "baseline")) {
  input_stage <- match.arg(input_stage)
  list(analysis_id = analysis_id, project_id = project_id,
       inputs = inputs, upstream_ids = as.character(upstream_ids),
       params_signature = params_signature, input_stage = input_stage)
}

#' Check that a binding is well formed before believing what it says
#'
#' Structure first, currency second. `all(logical(0))` is `TRUE`, so a
#' binding that names no inputs answered "all my inputs are current" -- which
#' is how an empty record passed for a verified one.
#'
#' @return `list(ok = , reason = )`.
#' @param pb A provenance binding.
#' @keywords internal
#' @export
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
  stage <- pb$input_stage %||% "effective"
  if (!stage %in% c("effective", "baseline")) {
    return(bad("input_stage names no known resolver"))
  }
  list(ok = TRUE, reason = "")
}

#' Are the inputs this result names still the project's current ones?
#'
#' Ownership is part of the question. Two projects built on the same file
#' have identical input versions, so comparing versions alone accepted one
#' project's results inside the other.
#'
#' @param pb A provenance binding.
#' @param project An `omics_project`.
#' @return Logical scalar.
#' @keywords internal
#' @export
result_inputs_current <- function(pb, project) {
  if (!isTRUE(validate_provenance_binding(pb)$ok)) return(FALSE)
  if (!identical(pb$project_id, project$metadata$project_id)) return(FALSE)
  tags <- names(pb$inputs)
  if (!all(tags %in% names(project$experiments))) return(FALSE)
  resolve <- if (identical(pb$input_stage %||% "effective", "baseline")) {
    baseline_input
  } else {
    effective_input
  }
  all(vapply(tags, function(tg) {
    identical(pb$inputs[[tg]], resolve(project, tg)$version)
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
#' @param project An `omics_project`.
#' @param bundle An `analysis_bundle`.
#' @param id Identifier to file it under; defaults to the bundle's own.
#' @return The project, with the analysis registered.
#' @keywords internal
#' @export
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
#' @param project An `omics_project`.
#' @return Character vector of analysis ids.
#' @keywords internal
#' @export
valid_result_ids <- function(project) {
  analyses <- project$analyses %||% list()
  ids <- names(analyses)
  if (length(ids) == 0L) return(character(0))
  # Coherence before currency. A result filed under a key that is not the
  # id it claims resolves an upstream reference to a different analysis,
  # and no amount of version comparison notices.
  broken <- validate_analysis_registry(project)$analysis_id
  ids <- setdiff(ids, broken)
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
#' @param project An `omics_project`.
#' @param max_history How many retired results to keep.
#' @param reason Recorded against each retired result.
#' @return The project, with invalid results moved to history.
#' @keywords internal
#' @export
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

#' Check the registry's own bookkeeping
#'
#' Three things the fixed-point algorithm cannot see, because it reasons
#' about which results are current rather than about whether the record is
#' coherent in the first place.
#'
#' The key a result is filed under has to be the id it claims. They are set
#' in two places -- `register_analysis()` takes an `id`, the binding carries
#' one -- and if they diverge, an upstream reference resolves to a different
#' analysis than the one it names, silently.
#'
#' The parameter signature has to match the parameters actually saved.
#' Otherwise a result can be re-filed under a different id after its
#' settings changed, and both copies look valid.
#'
#' And the graph has to be acyclic. A cycle makes every member invalid,
#' which is the right answer, but arrived at for the wrong reason: they are
#' not stale, the record is malformed, and saying so is the difference
#' between a bug report and a shrug.
#'
#' @param project An `omics_project`.
#' @return A data frame of `analysis_id` and `problem`; zero rows when the
#'   registry is coherent.
#' @keywords internal
#' @export
validate_analysis_registry <- function(project) {
  analyses <- project$analyses %||% list()
  problems <- list()
  add <- function(id, what) {
    problems[[length(problems) + 1L]] <<-
      data.frame(analysis_id = id, problem = what, stringsAsFactors = FALSE)
  }

  for (key in names(analyses)) {
    b <- analyses[[key]]
    pb <- b$provenance_binding
    if (is.null(pb)) {
      add(key, "no provenance binding")
      next
    }
    structure_ok <- validate_provenance_binding(pb)
    if (!isTRUE(structure_ok$ok)) {
      add(key, structure_ok$reason)
      next
    }
    if (!identical(key, pb$analysis_id)) {
      add(key, sprintf("filed under '%s' but claims '%s'", key, pb$analysis_id))
    }
    if (!identical(pb$params_signature, params_signature(b$params))) {
      add(key, "params signature does not match the saved parameters")
    }
  }

  # Reachability, iterated. A self-reference is caught by
  # validate_provenance_binding(); this finds the longer loops.
  edges <- lapply(analyses, function(b) b$provenance_binding$upstream_ids %||% character(0))
  for (start in names(edges)) {
    seen <- character(0)
    frontier <- edges[[start]]
    while (length(frontier) > 0L) {
      if (start %in% frontier) {
        add(start, "the dependency graph contains a cycle")
        break
      }
      nxt <- setdiff(frontier, seen)
      seen <- c(seen, nxt)
      frontier <- unique(unlist(edges[intersect(nxt, names(edges))]))
      frontier <- setdiff(frontier, seen)
    }
  }

  if (length(problems) == 0L) {
    return(data.frame(analysis_id = character(0), problem = character(0),
                      stringsAsFactors = FALSE))
  }
  do.call(rbind, problems)
}

#' Fork a project into a copy with its own identity
#'
#' Save-as is not the same file under a different name. Two projects that
#' share an id are one project as far as every ownership check is
#' concerned, so a copy that kept it would accept the original's results as
#' its own and the original would accept the copy's.
#'
#' This is the one place a result's ownership is rewritten. Everywhere else
#' an unattributed or foreign result is refused rather than adopted --
#' rewriting identity is how a record becomes a claim about data it was
#' never computed from, so it happens only here, deliberately, for every
#' result at once, and `derived_from` says where they came from.
#'
#' @param project An `omics_project`.
#' @param from_slug Where the original was saved, for the record.
#' @return The project, with a new identity and every result re-attributed.
#' @keywords internal
#' @export
fork_project <- function(project, from_slug = NA_character_) {
  old_id <- project$metadata$project_id
  new_id <- new_project_id()
  project$metadata$project_id <- new_id
  project$metadata$derived_from <- list(
    project_id = old_id,
    saved_as = from_slug,
    at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )

  retag <- function(b) {
    if (!is.null(b$provenance_binding) &&
        identical(b$provenance_binding$project_id, old_id)) {
      b$provenance_binding$project_id <- new_id
    }
    b
  }
  project$bundles  <- lapply(project$bundles %||% list(), retag)
  project$analyses <- lapply(project$analyses %||% list(), retag)
  # History too: a retired result whose owner no longer exists cannot be
  # told apart from one that was never attributed.
  project$retired  <- lapply(project$retired %||% list(), retag)

  # The ids are derived from the project id among other things, so they
  # change with it -- and the upstream references have to follow, or every
  # dependant names a node the copy does not contain.
  remap <- character(0)
  for (key in names(project$analyses)) {
    b <- project$analyses[[key]]
    fresh <- paste0("ana_", params_signature(list(
      view = b$analysis_name, inputs = b$provenance_binding$inputs,
      upstream = b$provenance_binding$upstream_ids,
      params = b$params, project = new_id)))
    remap[[key]] <- fresh
  }
  relabel <- function(b) {
    pb <- b$provenance_binding
    if (is.null(pb)) return(b)
    if (!is.na(remap[pb$analysis_id] %||% NA_character_)) {
      pb$analysis_id <- unname(remap[[pb$analysis_id]])
    }
    known <- intersect(pb$upstream_ids, names(remap))
    if (length(known)) {
      pb$upstream_ids[match(known, pb$upstream_ids)] <- unname(remap[known])
    }
    b$provenance_binding <- pb
    b
  }
  project$analyses <- stats::setNames(
    lapply(project$analyses, relabel), unname(remap[names(project$analyses)]))
  project$bundles <- lapply(project$bundles %||% list(), relabel)
  project
}

#' Bring a project saved before the contract existed into it
#'
#' Such a project carries results in `bundles` with nothing saying which data
#' they were computed on. There are three things one could do with them and
#' only one is honest.
#'
#' Adopting them -- writing the current project and input versions into their
#' binding -- makes them valid by assertion. That is the failure the whole
#' contract exists to prevent, arrived at during the migration instead of
#' during the analysis. Deleting them silently loses work the user can see on
#' screen. So they are moved to history with a reason, where the project view
#' can show them and the report and export will not treat them as current.
#'
#' Called on load, before any validity question is asked, so the general
#' algorithm never has to reason about a half-specified record.
#'
#' @param project An `omics_project`, freshly loaded.
#' @return The project, with unattributed results moved to history.
#' @keywords internal
#' @export
migrate_legacy_project <- function(project) {
  if (is.null(project$metadata$project_id)) {
    project$metadata$project_id <- new_project_id()
    project$metadata$schema_note <-
      "Opened from a file saved before results recorded which data they came from."
  }
  # Layers written before versions were stamped.
  project$experiments <- lapply(project$experiments, function(e) {
    if (is.null(e$version)) stamp_version(e) else e
  })

  bundles <- project$bundles %||% list()
  if (length(bundles) == 0L) return(project)
  unbound <- vapply(bundles, function(b) is.null(b$provenance_binding),
                    logical(1))
  if (!any(unbound)) return(project)

  history <- lapply(bundles[unbound], compress_for_history,
                    reason = "saved before results recorded which data they came from")
  # Keyed by generated id, not by the view slot they came from: history is
  # indexed by analysis and two migrations would otherwise both write "diff".
  # The slot survives on the record so the view can still be named.
  for (i in seq_along(history)) {
    history[[i]]$view_slot <- names(history)[[i]]
  }
  names(history) <- vapply(history, function(h) new_analysis_id(), character(1))

  project$retired <- c(project$retired %||% list(), history)
  project$bundles <- bundles[!unbound]
  project
}

# ---- per-layer preprocessing ------------------------------------------

#' Version of a layer, computed once and carried with it
#'
#' Cached on the object rather than recomputed per query: `effective_input()`
#' is read from reactive expressions, and hashing a 63k x 258 matrix on every
#' read is not a cache.
#'
#' @param input An `omics_input`.
#' @return The input, carrying its content version.
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
#' @param project An `omics_project`.
#' @param tag Experiment layer tag.
#' @param cleaned The preprocessed `omics_input`.
#' @param params The resolved parameters that produced it.
#' @param base_input_version Version of the layer it was built from.
#' @return The project, carrying that layer's applied preprocessing.
#' @keywords internal
#' @export
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
#' @param project An `omics_project`.
#' @param tag Experiment layer tag.
#' @return The project, with that layer's preprocessing removed.
#' @keywords internal
#' @export
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
#' @param project An `omics_project`.
#' @param tag Experiment layer tag.
#' @return `list(input, version, provenance)` for the imported layer.
#' @keywords internal
#' @export
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
#' @param project An `omics_project`.
#' @param tag Experiment layer tag.
#' @keywords internal
#' @export
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
#' @param session_id,project_id,activation,request_seq,view Session state at
#'   the moment the task was launched.
#' @param params_signature Signature of the parameters it was launched with.
#' @return A list carrying the session state a callback is checked against.
#' @keywords internal
#' @export
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
#' @export
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

# ---- how a layer came to be -------------------------------------------

#' Append a step to a layer's provenance
#'
#' One entry per operation that changed the values, with the arguments as
#' they resolved rather than as they were typed: `run_qc(impute_method =
#' NULL)` on proteomics is `MinProb`, and a script that replayed the `NULL`
#' would be replaying a different analysis on a different modality.
#'
#' The list is what `export_script()` writes out. Reconstructing it from the
#' final state instead is what dropped the normalization: the layer knew it
#' was on a log scale and not that anything had put it there, so the script
#' re-read the raw file and told the analysis it was already normalized.
#'
#' @param input An `omics_input`.
#' @param op Name of the function that was applied.
#' @param args Its resolved arguments.
#'
#' @return The input, with the step appended.
#' @keywords internal
#' @export
record_step <- function(input, op, args = list()) {
  step <- list(op = op, args = args)
  input$provenance <- c(input$provenance %||% list(), list(step))
  input
}

#' Change a layer's assay type, and its record of it
#'
#' `assay_type` is the only note of what scale `expr_mat` is on, and the
#' analysis backends act on it: `log2(x + 1)` is applied for
#' `"raw_count"` and not otherwise, and limma refuses raw counts outright.
#' So the import offers the user a chance to correct a wrong guess.
#'
#' Correcting it used to be a single assignment to the field, which left
#' the provenance still recording the value the file was read with. The
#' analysis then ran on the corrected layer and the exported script
#' replayed the original -- so a layer relabelled from `raw_count` to
#' `logcpm` analysed fine in the application and the downloaded script
#' stopped at the scale check. The results were right and could not be
#' reproduced, which is the harder failure to notice.
#'
#' Both together, here, so they cannot drift apart again.
#'
#' @param input An `omics_input`.
#' @param assay_type The corrected assay type.
#'
#' @return The input, relabelled, with its `read_omics` step updated.
#' @export
#' @family omics_input
relabel_assay_type <- function(input, assay_type) {
  validate_omics_input(input)
  assert_string(assay_type, "assay_type")
  # Refused, not warned about. validate_assay_type() warns because its
  # other job is a magnitude heuristic that a small or unusual matrix can
  # legitimately trip. This is not that: an assay type outside the
  # modality's list is a category error, the picker never offers one, and
  # a layer silently relabelled `raw_intensity` on RNA-seq is mislabelled
  # in a way nothing downstream catches.
  allowed <- SUPPORTED_ASSAY_TYPES[[input$omics_type %||% ""]]
  if (!is.null(allowed) && !assay_type %in% allowed) {
    stop("`assay_type = \"", assay_type, "\"` is not one an ",
         input$omics_type, " layer can be on. Expected one of: ",
         paste(sprintf("'%s'", allowed), collapse = ", "), ".",
         call. = FALSE)
  }
  # Still runs, for the magnitude warning a legitimate relabel may earn.
  validate_assay_type(assay_type, input$omics_type)
  input$assay_type <- assay_type
  steps <- input$provenance %||% list()
  for (i in seq_along(steps)) {
    if (identical(steps[[i]]$op, "read_omics")) {
      steps[[i]]$args$assay_type <- assay_type
    }
  }
  input$provenance <- steps
  input
}

#' The steps recorded on a layer
#'
#' @param input An `omics_input`.
#' @return A list of `list(op, args)`, empty when nothing was recorded.
#' @keywords internal
#' @export
provenance_steps <- function(input) {
  input$provenance %||% list()
}
