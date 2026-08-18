# Async helper: wrap a long-running expression in a future and route
# success / error back to the caller via callbacks. Used by the Diff,
# Enrich, and Integration views to keep the Shiny session responsive
# during computation.
#
# In testServer (or when `future::plan()` is `sequential`), the
# future runs synchronously, so existing testServer tests continue
# to pass without modification.

run_async <- function(func, on_success, on_error, message = "Running...") {
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

  # Belt-and-braces close: if the session terminates before the
  # promise settles, onSessionEnded fires and we release the
  # progress handle. Using `tryCatch` around close() makes a
  # second call (from the promise handler below) a no-op.
  close_progress <- function() {
    tryCatch(progress$close(), error = function(e) NULL)
  }
  session <- shiny::getDefaultReactiveDomain()
  if (!is.null(session)) {
    session$onSessionEnded(close_progress)
  }

  f <- future::future({
    tryCatch(func(), error = function(e) e)
  }, seed = TRUE)
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
