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
                      .future = future::future,
                      on_cancel = function() on_error(CANCELLED_MESSAGE)) {
  # testServer does not have a real event loop; run synchronously
  # so that existing testServer tests continue to pass.
  if (isTRUE(getOption("shiny.allowoutputreads", FALSE))) {
    result <- tryCatch(func(), error = function(e) e)
    if (inherits(result, "error")) {
      on_error(conditionMessage(result))
    } else {
      on_success(result)
    }
    return(invisible(NULL))
  }

  session <- shiny::getDefaultReactiveDomain()
  task <- async_task_new(session, message)
  func <- with_step_file(func, task$step_file)

  # future() can throw before anything is submitted -- most reliably by
  # refusing to export globals over future.globals.maxSize. Thrown from
  # inside an observer that is what kills the session and greys the
  # page, so it is routed to on_error like any other failure.
  f <- tryCatch(
    .future({ tryCatch(func(), error = function(e) e) }, seed = TRUE),
    error = function(e) e
  )
  if (inherits(f, "error")) {
    task$finish()
    on_error(conditionMessage(f))
    return(invisible(NULL))
  }
  task$future <- f
  task$on_cancel <- on_cancel

  p <- promises::as.promise(f)
  promises::then(
    p,
    onFulfilled = function(result) {
      # A cancelled run has already told its caller; whatever the worker
      # finished with afterwards is not wanted.
      if (task$cancelled) return(invisible(NULL))
      task$finish()
      if (inherits(result, "error")) {
        on_error(conditionMessage(result))
      } else {
        on_success(result)
      }
    },
    onRejected = function(err) {
      if (task$cancelled) return(invisible(NULL))
      task$finish()
      on_error(conditionMessage(err))
    }
  )
  invisible(task)
}

CANCELLED_MESSAGE <- "The run was cancelled."

# ---- progress and cancel -----------------------------------------------
# What a long run shows while it runs. It used to be a progress bar set
# to 30% and left there, with no way to stop the run: a DESeq2 fit on a
# big matrix looked exactly like a hung page. Now it is a notification
# that says what is running and for how long, with a bar that keeps
# moving (it cannot know the true fraction -- the engines do not report
# one -- so it approaches the end without claiming it), and a Cancel
# button that interrupts the worker where the future backend supports
# it and, either way, stops the result from being used.

async_task_new <- function(session, message) {
  task <- new.env(parent = emptyenv())
  task$cancelled <- FALSE
  task$done <- FALSE
  task$future <- NULL
  task$on_cancel <- NULL
  # Where the worker says which step it is on (see with_step_file()).
  task$step_file <- tempfile("omics-step-", fileext = ".txt")
  task$finish <- function() {
    if (task$done) return(invisible(NULL))
    task$done <- TRUE
    unlink(task$step_file)
    if (!is.null(task$ticker)) task$ticker$destroy()
    if (!is.null(session)) {
      tryCatch(shiny::removeNotification(task$id, session = session),
               error = function(e) NULL)
    }
    invisible(NULL)
  }
  if (is.null(session)) return(task)

  registry <- async_registry(session)
  registry$n <- registry$n + 1L
  task$id <- sprintf("async-task-%d", registry$n)
  registry$tasks[[task$id]] <- task
  started <- Sys.time()

  show <- function() {
    secs <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    shiny::showNotification(
      async_progress_ui(message, secs, task$id, step = read_step(task$step_file)),
      id = task$id, duration = NULL, closeButton = FALSE,
      type = "message", session = session)
  }
  show()
  task$ticker <- shiny::observe({
    shiny::invalidateLater(1000, session)
    if (!task$done) show()
  }, domain = session)
  session$onSessionEnded(function() {
    task$done <- TRUE
    if (!is.null(task$ticker)) task$ticker$destroy()
  })
  task
}

# One per session: the tasks it has running, and the observer that hears
# their Cancel buttons (a notification lives outside every module, so
# the click arrives as a top-level input).
async_registry <- function(session) {
  root <- session$rootScope()
  reg <- root$userData$async_registry
  if (!is.null(reg)) return(reg)
  reg <- new.env(parent = emptyenv())
  reg$n <- 0L
  reg$tasks <- list()
  root$userData$async_registry <- reg
  shiny::observeEvent(root$input$omics_async_cancel, {
    id <- root$input$omics_async_cancel
    task <- reg$tasks[[id]]
    if (is.null(task) || task$done) return()
    async_cancel(task)
  }, domain = root)
  reg
}

async_cancel <- function(task) {
  task$cancelled <- TRUE
  if (!is.null(task$future) && "cancel" %in% getNamespaceExports("future")) {
    tryCatch(future::cancel(task$future), error = function(e) NULL,
             warning = function(w) NULL)
  }
  task$finish()
  if (is.function(task$on_cancel)) task$on_cancel()
  invisible(TRUE)
}

async_progress_ui <- function(message, secs, id, step = NULL) {
  # 1 - exp(-t / 20): half-way at ~14 s, 90% at ~46 s, never 100%.
  frac <- 0.05 + 0.9 * (1 - exp(-secs / 20))
  htmltools::tags$div(
    class = "async-progress",
    htmltools::tags$div(
      class = "async-progress-head",
      htmltools::tags$strong(message),
      htmltools::tags$span(class = "muted", format_elapsed(secs))
    ),
    if (length(step) && nzchar(step)) {
      htmltools::tags$div(class = "async-progress-step muted", step)
    },
    htmltools::tags$div(
      class = "async-progress-track",
      htmltools::tags$div(class = "async-progress-bar",
                          style = sprintf("width:%.0f%%", 100 * frac))
    ),
    htmltools::tags$button(
      type = "button", class = "btn btn-sm btn-link async-cancel",
      onclick = sprintf(
        "Shiny.setInputValue('omics_async_cancel', '%s', {priority: 'event'});",
        id),
      "Cancel")
  )
}

# The engine reports the step it starts (omicsCore's report_progress(),
# read through the `omicsCore.progress` option). A worker is another
# process, so the step travels through a file: the worker writes it,
# and the notification, redrawn every second, reads it. The file is
# written whole and renamed into place, so a read never sees half a line.
with_step_file <- function(func, step_file) {
  detached_call(
    function() {
      options(omicsCore.progress = function(step) {
        tmp <- paste0(step_file, ".tmp")
        writeLines(step, tmp)
        file.rename(tmp, step_file)
      })
      func()
    },
    func = func, step_file = step_file
  )
}

read_step <- function(step_file) {
  if (is.null(step_file) || !file.exists(step_file)) return(NULL)
  tryCatch(readLines(step_file, n = 1L, warn = FALSE), error = function(e) NULL)
}

# The same steps for work done in the main process (QC), as a Shiny
# progress panel: its messages go out while the computation still runs,
# which an output or notification redrawn by the server would not.
with_step_progress <- function(message, expr,
                               session = shiny::getDefaultReactiveDomain()) {
  if (is.null(session) || isTRUE(getOption("shiny.allowoutputreads", FALSE))) {
    return(expr)
  }
  shiny::withProgress(message = message, value = 0.2, session = session, {
    old <- options(omicsCore.progress = function(step) {
      shiny::setProgress(detail = step, session = session)
    })
    on.exit(options(old), add = TRUE)
    expr
  })
}

format_elapsed <- function(secs) {
  secs <- max(0, round(secs))
  if (secs < 60) sprintf("%d s", secs) else sprintf("%d min %02d s", secs %/% 60, secs %% 60)
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

#' A counter that tells a current async result from a stale one
#'
#' `start()` takes a number for a new run; `bump()` takes one for a reset
#' (a layer or project change) with no run behind it. A callback asks
#' `is_current(id)` before writing its result, so a run that finishes
#' after the user moved on does not land on the new state.
#' `is_last_started(id)` says whether no newer run has begun, which is
#' when the busy state may be cleared.
#'
#' @return A list of closures.
#' @keywords internal
#' @noRd
run_epoch <- function() {
  current <- 0L
  last_started <- 0L
  list(
    start = function() {
      current <<- current + 1L
      last_started <<- current
      current
    },
    bump = function() {
      current <<- current + 1L
      invisible(current)
    },
    is_current = function(id) identical(id, current),
    is_last_started = function(id) identical(id, last_started)
  )
}
