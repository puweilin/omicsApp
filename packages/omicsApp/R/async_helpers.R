# Async helper: wrap a long-running expression in a future and route
# success / error back to the caller via callbacks. Used by the Diff,
# Enrich, and Integration views.
#
# What this does and does not buy, stated accurately because the comment
# here used to claim more: with plain promises, the session that started
# the task still blocks on it. Other sessions stay responsive, which is
# what matters on a shared server, but the analyst who pressed Run cannot
# do anything else meanwhile. Making that session responsive too needs
# ExtendedTask; this is not that.
#
# In testServer (or when `future::plan()` is `sequential`), the
# future runs synchronously, so existing testServer tests continue
# to pass without modification.

#' Run `func` in a background worker
#'
#' @param func A function of no arguments. **Its environment is
#'   serialised along with it.** Build it with [detached_call()] rather
#'   than defining it inline in a module server, or the whole reactive
#'   scope -- the previous result, the project, every fixture -- travels
#'   to the worker with it.
#' @param on_success,on_error Callbacks.
#' @param message Progress label.
#' @param .future The future constructor. Injected rather than called
#'   through `future::` so a test can make it throw: the failure worth
#'   covering here is future refusing *before* it submits anything, and
#'   provoking that for real needs a payload too large to put in a test.
#' @keywords internal
#' @noRd
#' @param admissible A function of one argument -- `"success"` or
#'   `"error"` -- returning `list(ok = , reason = )`. Called when the task
#'   comes back, before either callback runs. A task that started against
#'   data the project no longer holds must not write its result, and must
#'   not write its *error* either: an error message about a run the user has
#'   already moved past reads as a fault in what is on screen now.
#' @param on_discard Called with the reason when a result is refused, so the
#'   view can say something rather than appearing to have done nothing.
run_async <- function(func, on_success, on_error, message = "Running...",
                      admissible = function(outcome) list(ok = TRUE, reason = ""),
                      on_discard = function(reason) NULL,
                      progress_poll = NULL,
                      .future = future::future) {
  # Both outcomes go through the same gate. Wrapping only on_success would
  # leave the failure path able to overwrite state the success path is no
  # longer allowed to touch.
  guard <- function(outcome, act) {
    verdict <- tryCatch(admissible(outcome),
                        error = function(e) list(ok = FALSE,
                                                 reason = conditionMessage(e)))
    if (isTRUE(verdict$ok)) act() else on_discard(verdict$reason)
  }
  # testServer does not have a real event loop; run synchronously
  # so that existing testServer tests continue to pass.
  if (isTRUE(getOption("shiny.allowoutputreads", FALSE))) {
    result <- tryCatch(func(), error = function(e) e)
    if (inherits(result, "error")) {
      guard("error", function() on_error(conditionMessage(result)))
    } else {
      guard("success", function() on_success(result))
    }
    return()
  }

  progress <- shiny::Progress$new()
  progress$set(message = message, value = 0.3)

  # Closed at most once. Progress$close() warns rather than errors on a
  # second call, so a tryCatch(error=) around it does not stop the
  # warning reaching the log -- five of them arrived with the failure
  # this guard was added for.
  closed <- FALSE
  close_progress <- function() {
    if (closed) return(invisible(NULL))
    closed <<- TRUE
    tryCatch(progress$close(), error = function(e) NULL,
             warning = function(w) NULL)
  }
  session <- shiny::getDefaultReactiveDomain()
  if (!is.null(session)) {
    session$onSessionEnded(close_progress)
  }

  # A bar that moves, when the caller can say where the work has got to.
  #
  # The task runs in another process, so it cannot touch this Progress --
  # it reports through whatever channel `progress_poll` reads, and this
  # side asks every half second. Without it the bar sits at 30% for the
  # length of the run, which on the enrichment grid is a minute of a
  # screen that looks exactly like a screen that has stopped.
  if (is.function(progress_poll)) {
    poll_obs <- shiny::observe({
      shiny::invalidateLater(500)
      st <- tryCatch(progress_poll(), error = function(e) NULL)
      if (is.null(st)) return()
      progress$set(value = st$value, message = message, detail = st$detail)
    })
    # Destroyed with the bar: an observer left running after the task
    # finishes keeps invalidating for the life of the session.
    local({
      inner <- close_progress
      close_progress <<- function() {
        tryCatch(poll_obs$destroy(), error = function(e) NULL)
        inner()
      }
    })
  }

  # future() can throw before anything is submitted -- most reliably by
  # refusing to export globals over future.globals.maxSize. Thrown from
  # inside an observer that is what kills the session and greys the
  # page, so it is routed to on_error like any other failure.
  f <- tryCatch(
    .future({ tryCatch(func(), error = function(e) e) }, seed = TRUE),
    error = function(e) e
  )
  if (inherits(f, "error")) {
    close_progress()
    guard("error", function() on_error(conditionMessage(f)))
    return(invisible(NULL))
  }

  p <- promises::as.promise(f)
  promises::then(
    p,
    onFulfilled = function(result) {
      if (inherits(result, "error")) {
        guard("error", function() on_error(conditionMessage(result)))
      } else {
        guard("success", function() on_success(result))
      }
      close_progress()
    },
    onRejected = function(err) {
      guard("error", function() on_error(conditionMessage(err)))
      close_progress()
    }
  )
}

#' Report progress across a process boundary, through a file
#'
#' The worker cannot reach the session's `Progress`: it is another
#' process, and under `multisession` it does not even share its
#' `tempdir()`. What both ends do share is the filesystem, so the parent
#' picks a path, hands it over, and the worker appends a line per step.
#'
#' Append-per-step rather than rewrite-per-step, so a read that lands
#' mid-write sees one fewer step rather than a truncated number. The
#' reader takes the last complete line and ignores a partial one.
#'
#' @param path A file the worker can write and the session can read.
#' @return A `function(done, total, label)` suitable for
#'   `run_enrichment_grid(on_progress = )`.
#' @keywords internal
#' @noRd
progress_writer <- function(path) {
  force(path)
  function(done, total, label) {
    line <- sprintf("%d\t%d\t%s\n", as.integer(done), as.integer(total),
                    label %||% "")
    tryCatch(cat(line, file = path, append = TRUE), error = function(e) NULL)
  }
}

#' @rdname progress_writer
#' @return For `progress_reader()`, a `function()` giving
#'   `list(value = , detail = )` or `NULL` before the first step lands.
#' @keywords internal
#' @noRd
progress_reader <- function(path) {
  force(path)
  function() {
    if (!file.exists(path)) return(NULL)
    size <- file.size(path)
    if (is.na(size) || size == 0) return(NULL)
    raw <- tryCatch(readChar(path, size, useBytes = TRUE),
                    error = function(e) "")
    if (!nzchar(raw)) return(NULL)

    # Only lines a newline has closed. `readLines()` hands back an
    # unterminated final line without complaint, and a read that landed
    # halfway through an append would then be parsed as a complete step
    # -- "7\t1" from a torn "7\t14\tGSEA ..." reads as 7 of 1, a bar that
    # jumps to full and stays there. A step behind is the right answer.
    lines <- strsplit(raw, "\n", fixed = TRUE)[[1L]]
    if (!endsWith(raw, "\n")) lines <- utils::head(lines, -1L)

    parts <- NULL
    for (i in rev(seq_along(lines))) {
      p <- strsplit(lines[[i]], "\t", fixed = TRUE)[[1L]]
      if (length(p) >= 2L && !is.na(suppressWarnings(as.integer(p[[1L]]))) &&
          !is.na(suppressWarnings(as.integer(p[[2L]])))) {
        parts <- p
        break
      }
    }
    if (is.null(parts)) return(NULL)
    done <- as.integer(parts[[1L]])
    total <- max(as.integer(parts[[2L]]), 1L)
    list(value = min(done / total, 1),
         detail = sprintf("%d / %d  %s", done, total,
                          if (length(parts) >= 3L) parts[[3L]] else ""))
  }
}

#' Build a zero-argument function that carries only what it is given
#'
#' A closure defined in a module server keeps that server's environment
#' as its parent, and `future` serialises the whole chain: the previous
#' analysis bundle, the loaded project, the demo fixtures. A differential
#' run measured 527 MB of "globals" this way, of which 501 MB was the
#' function itself, and future refused to export it.
#'
#' Re-parenting to `baseenv()` makes the payload exactly the named
#' values. Anything the body needs must therefore be named -- a
#' forgotten one becomes "object not found" in the worker rather than a
#' silent capture, which is the failure worth having.
#'
#' @param fn A function of no arguments.
#' @param ... Named values the body refers to.
#' @keywords internal
#' @noRd
detached_call <- function(fn, ...) {
  environment(fn) <- list2env(list(...), parent = baseenv())
  fn
}


#' A diff bundle with the parts the worker will not read removed
#'
#' `detached_call()` bounds *which* values travel; this bounds how big one
#' of them is. Both are needed, and the second was missing: `run_diff()` now
#' defaults to dropping the fitted model, but a default cannot reach an
#' object that already exists. A project saved before that change, restored
#' from its `.omp`, still carries the model -- and so does any bundle from
#' `keep_model = TRUE`, which is a supported thing to ask for.
#'
#' Measured on a restored bundle of that shape, 63k features x 48 samples:
#' `future` refused the enrichment payload at 802 MiB. The same bundle
#' through here is 25 MB, and stays there whatever the model would have
#' weighed, because what is left scales with features rather than samples.
#'
#' Subtraction rather than reconstruction. Building a "lightweight bundle"
#' from the fields enrichment is known to read would have dropped
#' `provenance_binding`, which `run_integration()` reads through
#' `side_analysis_id()` and enrichment does not -- the sort of omission that
#' surfaces as a blank column in a table nobody checks. Taking away two
#' named elements cannot silently lose a third.
#'
#' `diff_result_df` is kept whole. Filtering it here would change the
#' answer, not the payload: ORA's universe is every tested feature, and
#' GSEA ranks the full list.
#'
#' @param bundle An `analysis_bundle` from `run_diff()`, or `NULL`.
#' @keywords internal
#' @noRd
bundle_for_worker <- function(bundle) {
  if (!inherits(bundle, "analysis_bundle")) return(bundle)
  # The fitted model: nothing downstream of run_diff() reads it.
  bundle$results$diff_object <- NULL
  # The backend-native table: neither run_enrichment() nor
  # run_integration() reads it, and it is the standardized one they use.
  bundle$results$diff_raw_df <- NULL
  bundle
}

#' A project carrying its data but not its analyses
#'
#' `run_integration()` wants the project for one thing --
#' `resolve_experiment_pair()`, which reads `project$experiments`. The
#' bundles it works from arrive separately in `diff_bundles`.
#'
#' `project$bundles` therefore travels for nothing, and it is the largest
#' thing in a project that has been used: every analysis run against every
#' layer, including the diff bundles that are being passed in properly
#' alongside. On a restored project that is where the fitted models are.
#'
#' `assays = FALSE` goes further and drops the layers' data as well,
#' leaving the project's shape and its layer names. That is legitimate for
#' exactly one caller: the concordance path reads the project through
#' `resolve_experiment_pair()`, which asks `experiment_tags()` for names
#' and touches nothing else. The two result tables it works from arrive in
#' `diff_bundles`.
#'
#' It is also the difference between running and not: dropping the bundles
#' alone still left `proj` at 277 MiB of a 600 MiB payload, because on a
#' 63k x 129 layer the data *is* the project.
#'
#' A stripped project is deliberately unusable for anything that needs the
#' matrices -- `run_integration(method = "correlation")` would get an empty
#' one and stop. That is the right failure: loud, immediate, and at the
#' call site that asked for it. `test-integration.R` holds the other end,
#' asserting concordance gives the same answer either way, so this stops
#' being safe the moment that stops being true.
#'
#' @param project An `omics_project`, or `NULL`.
#' @param assays Keep the layers' matrices? `TRUE` by default.
#' @keywords internal
#' @noRd
project_for_worker <- function(project, assays = TRUE) {
  if (!inherits(project, "omics_project")) return(project)
  project$bundles <- NULL
  if (!isTRUE(assays)) {
    project$experiments <- lapply(project$experiments, function(e) {
      # Emptied rather than removed: `experiments` still has to be a list
      # of layers with the right names for experiment_tags() to read.
      e$expr_mat   <- e$expr_mat[0L, 0L, drop = FALSE]
      e$raw_mat    <- NULL
      e$feature_df <- e$feature_df[0L, , drop = FALSE]
      e$meta_df    <- e$meta_df[0L, , drop = FALSE]
      e
    })
  }
  project
}


#' A gate that admits only the newest request for a view
#'
#' The narrow half of the contract in omicsCore: a request number per view,
#' compared when the task returns. `callback_admissible()` is the full
#' check, and where a project is in scope the views use it; this covers the
#' case that needs no project -- the user pressing Run again while the
#' previous run is still out.
#'
#' Per view, not global: enrichment moving on does not make a differential's
#' outstanding request stale.
#'
#' @param seq_val The request number captured when the task was launched.
#' @param current A function returning the view's current request number.
#'
#' @return A function suitable for `run_async(admissible = )`.
#' @keywords internal
#' @noRd
newest_request_only <- function(seq_val, current) {
  function(outcome) {
    now <- tryCatch(current(), error = function(e) NULL)
    if (is.null(now) || identical(now, seq_val)) {
      return(list(ok = TRUE, reason = ""))
    }
    list(ok = FALSE, reason = "a newer run was started")
  }
}


#' Tell the user a late result was thrown away
#'
#' The three views all need this and it reads the same in each: a run that
#' vanishes without a word is indistinguishable from one that never
#' started, and the user is left waiting for something that has already
#' been decided against.
#'
#' @return A function suitable for `run_async(on_discard = )`.
#' @keywords internal
#' @noRd
discard_notice <- function() {
  function(reason) {
    shiny::showNotification(
      sprintf("An earlier run finished after you started another; its result was discarded (%s).",
              reason),
      type = "message", duration = 6)
  }
}
