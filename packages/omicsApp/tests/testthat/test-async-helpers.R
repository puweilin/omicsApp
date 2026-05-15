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
