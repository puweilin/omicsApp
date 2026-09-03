# Async helper: wrap a long-running expression in a future and route
# success / error back to the caller via callbacks. Used by the Diff,
# Enrich, and Integration views to keep the Shiny session responsive
# during computation.
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
run_async <- function(func, on_success, on_error, message = "Running...",
                      .future = future::future) {
  # testServer does not have a real event loop; run synchronously
  # so that existing testServer tests continue to pass.
  if (isTRUE(getOption("shiny.allowoutputreads", FALSE))) {
    result <- tryCatch(func(), error = function(e) e)
    if (inherits(result, "error")) {
      on_error(conditionMessage(result))
    } else {
      on_success(result)
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
    on_error(conditionMessage(f))
    return(invisible(NULL))
  }

  p <- promises::as.promise(f)
  promises::then(
    p,
    onFulfilled = function(result) {
      if (inherits(result, "error")) {
        on_error(conditionMessage(result))
      } else {
        on_success(result)
      }
      close_progress()
    },
    onRejected = function(err) {
      on_error(conditionMessage(err))
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
