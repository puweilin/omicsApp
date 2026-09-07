# Tests for async_helpers.R — run_async function.

# Enable synchronous test mode so run_async does not try to create a
# shiny::Progress outside of a real Shiny session.
op <- options(shiny.allowoutputreads = TRUE)
on.exit(options(op), add = TRUE)

test_that("run_async synchronous fallback calls on_success", {
  # In testServer context (shiny.allowoutputreads = TRUE), run_async
  # runs synchronously and calls on_success with the result.
  result <- NULL
  error_msg <- NULL
  run_async(
    function() 42L,
    on_success = function(x) { result <<- x },
    on_error = function(msg) { error_msg <<- msg },
    message = "Testing..."
  )
  expect_equal(result, 42L)
  expect_null(error_msg)
})

test_that("run_async synchronous fallback calls on_error on failure", {
  result <- NULL
  error_msg <- NULL
  run_async(
    function() stop("boom"),
    on_success = function(x) { result <<- x },
    on_error = function(msg) { error_msg <<- msg },
    message = "Testing..."
  )
  expect_null(result)
  expect_match(error_msg, "boom")
})

test_that("run_async passes with NULL result", {
  result <- NULL
  error_msg <- NULL
  run_async(
    function() NULL,
    on_success = function(x) { result <<- x },
    on_error = function(msg) { error_msg <<- msg },
    message = "Testing..."
  )
  expect_null(result)
  expect_null(error_msg)
})

test_that("run_async passes with character result", {
  result <- NULL
  run_async(
    function() "success",
    on_success = function(x) { result <<- x },
    on_error = function(msg) { },
    message = "Testing..."
  )
  expect_equal(result, "success")
})

test_that("run_async handles error with custom message", {
  msg <- NULL
  run_async(
    function() stop("custom failure"),
    on_success = function(x) { },
    on_error = function(m) { msg <<- m },
    message = "Custom message"
  )
  expect_equal(msg, "custom failure")
})

test_that("run_async passes list result through", {
  result <- NULL
  run_async(
    function() list(a = 1, b = 2, c = list(d = 3)),
    on_success = function(x) { result <<- x },
    on_error = function(m) { },
    message = "Testing..."
  )
  expect_equal(result$a, 1)
  expect_equal(result$c$d, 3)
})

test_that("run_async preserves closure over local variables", {
  local_var <- "hello"
  captured <- NULL
  run_async(
    function() paste(local_var, "world"),
    on_success = function(x) { captured <<- x },
    on_error = function(m) { },
    message = "Testing..."
  )
  expect_equal(captured, "hello world")
})

# --- a result that arrives after the user has moved on --------------------

# Both outcomes go through the gate. Guarding only success leaves the
# failure path able to overwrite state the success path is no longer
# allowed to touch: an error about a run the user has already replaced
# reads as a fault in what is on screen now.
test_that("a stale result is discarded, and said so", {
  wrote <- NULL
  discarded <- NULL
  run_async(
    func = function() "late result",
    on_success = function(x) wrote <<- x,
    on_error = function(m) wrote <<- m,
    admissible = function(outcome) list(ok = FALSE, reason = "a newer run was started"),
    on_discard = function(reason) discarded <<- reason
  )
  expect_null(wrote)
  expect_match(discarded, "newer run")
})

test_that("a stale failure is discarded too", {
  wrote <- NULL
  discarded <- NULL
  run_async(
    func = function() stop("the engine fell over"),
    on_success = function(x) wrote <<- x,
    on_error = function(m) wrote <<- m,
    admissible = function(outcome) list(ok = FALSE, reason = "stale"),
    on_discard = function(reason) discarded <<- reason
  )
  expect_null(wrote)
  expect_equal(discarded, "stale")
})

test_that("the gate sees which outcome it is deciding about", {
  seen <- character(0)
  run_async(function() 1, function(x) NULL, function(m) NULL,
            admissible = function(outcome) {
              seen <<- c(seen, outcome); list(ok = TRUE, reason = "")
            })
  run_async(function() stop("x"), function(x) NULL, function(m) NULL,
            admissible = function(outcome) {
              seen <<- c(seen, outcome); list(ok = TRUE, reason = "")
            })
  expect_equal(seen, c("success", "error"))
})

test_that("a gate that throws refuses rather than letting the result through", {
  wrote <- NULL
  run_async(function() "x", function(x) wrote <<- x, function(m) wrote <<- m,
            admissible = function(outcome) stop("cannot tell"),
            on_discard = function(reason) NULL)
  expect_null(wrote)
})

test_that("newest_request_only admits the current request and no other", {
  expect_true(newest_request_only(3L, function() 3L)("success")$ok)
  expect_false(newest_request_only(2L, function() 3L)("success")$ok)
  expect_match(newest_request_only(2L, function() 3L)("error")$reason, "newer run")
  # A counter that cannot be read is not evidence of staleness.
  expect_true(newest_request_only(2L, function() stop("gone"))("success")$ok)
})

test_that("by default everything is admitted, so existing callers are unchanged", {
  got <- NULL
  run_async(function() "ok", function(x) got <<- x, function(m) got <<- m)
  expect_equal(got, "ok")
})
