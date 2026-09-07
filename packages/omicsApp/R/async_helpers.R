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
