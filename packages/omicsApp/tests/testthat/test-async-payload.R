# future serialises the function it is handed, and serialising a closure
# takes its environment with it. A closure written inline in a module
# server therefore carries that server's whole scope to the worker: the
# previous analysis bundle, the loaded project, every demo fixture. On a
# real workbook that reached 527 MB of "globals", future refused to
# export it, and the throw -- from inside an observer -- ended the
# session and greyed the page.

test_that("a detached call carries only what it is handed", {
  heavy <- matrix(rnorm(2e5), ncol = 10)
  small <- 3L
  inline   <- function() small + 1L          # keeps `heavy` in scope too
  detached <- detached_call(function() small + 1L, small = small)

  size <- function(f) length(serialize(f, NULL))
  expect_lt(size(detached), size(inline) / 10)
  expect_identical(detached(), inline())
})

test_that("a value the body needs but was not named is an error, not a capture", {
  # The point of re-parenting: a forgotten value fails loudly in the
  # worker instead of being smuggled along with everything around it.
  hidden <- 42L
  fn <- detached_call(function() hidden + 1L)   # not passed
  expect_error(fn(), "hidden")
})

test_that("the detached function still reaches package code", {
  # baseenv() as parent still resolves ::, which is how the callers
  # invoke omicsCore.
  fn <- detached_call(function() omicsCore::is_omics_input(x), x = 1)
  expect_false(fn())
})

# ---- failures must not reach the session ------------------------------

# What future does when the payload is too large, without needing a
# payload that large: it throws from getGlobalsAndPackages, before
# anything is submitted, so the tryCatch *inside* the future never runs.
refusing_future <- function(...) {
  stop("The total size of the 7 globals exported for future expression ",
       "is 527.50 MiB.. This exceeds the maximum allowed size of ",
       "500.00 MiB (option 'future.globals.maxSize')")
}

in_session <- function(code) {
  shiny::withReactiveDomain(shiny::MockShinySession$new(), code)
}

test_that("a refusal to export globals is reported, not thrown", {
  # Thrown from inside an observer this ended the session and greyed
  # the page; the user saw a dead tab and no message.
  got <- NULL
  expect_silent(in_session(
    run_async(function() 1L,
              on_success = function(x) got <<- "success",
              on_error   = function(msg) got <<- msg,
              .future    = refusing_future)
  ))
  expect_true(is.character(got))
  expect_match(got, "future.globals.maxSize", fixed = TRUE)
})

test_that("progress is closed once, however the run ends", {
  # Progress$close() warns rather than errors on a second call, so the
  # old tryCatch(error=) let five warnings through on the failing path.
  expect_no_warning(in_session(
    run_async(function() 1L,
              on_success = function(x) NULL,
              on_error   = function(msg) NULL,
              .future    = refusing_future)
  ))
})

test_that("a normal run still succeeds and still closes its progress", {
  skip_if_not_installed("later")
  old <- future::plan(future::sequential)
  on.exit(future::plan(old), add = TRUE)
  got <- NULL
  expect_no_warning(in_session({
    run_async(detached_call(function() n * 2L, n = 21L),
              on_success = function(x) got <<- x,
              on_error   = function(msg) got <<- msg)
    # then() schedules its callback on the event loop, which nothing
    # runs outside a live app. Drain it, or the assertion below reads
    # the value from before the promise settled.
    while (!later::loop_empty()) later::run_now()
  }))
  expect_identical(got, 42L)
})

test_that("a running task can be cancelled, and its late result is dropped", {
  skip_if_not_installed("later")
  old <- future::plan(future::sequential)
  on.exit(future::plan(old), add = TRUE)
  got <- character(0)
  session <- shiny::MockShinySession$new()
  shiny::withReactiveDomain(session, {
    task <- run_async(detached_call(function() 1L),
                      on_success = function(x) got <<- c(got, "success"),
                      on_error   = function(msg) got <<- c(got, msg))
    # Cancelled before the promise settles: the caller hears "cancelled"
    # once, and the finished result that arrives afterwards is ignored.
    async_cancel(task)
    while (!later::loop_empty()) later::run_now()
  })
  expect_identical(got, CANCELLED_MESSAGE)
  expect_true(task$done)
})

test_that("the progress notification says how long, and offers Cancel", {
  html <- paste(as.character(async_progress_ui("Running limma", 75, "async-task-3")),
                collapse = "")
  expect_match(html, "1 min 15 s", fixed = TRUE)
  expect_match(html, "Cancel", fixed = TRUE)
  expect_match(html, "omics_async_cancel", fixed = TRUE)
  expect_match(html, "async-task-3", fixed = TRUE)
  # Moves with time, never claims completion.
  w <- function(s) as.numeric(sub(".*width:([0-9]+)%.*", "\\1",
                                  paste(as.character(async_progress_ui("x", s, "i")), collapse = "")))
  expect_lt(w(1), w(30))
  expect_lt(w(3600), 100)
})
