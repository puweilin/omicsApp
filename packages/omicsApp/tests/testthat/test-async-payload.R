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

test_that("bundle_for_worker drops the model and keeps everything else", {
  b <- omicsCore::new_analysis_bundle(
    analysis_name = "run_diff",
    input_info = list(omics_type = "rnaseq", n_features = 3L),
    params = list(comparison = "case_vs_ctrl", method = "deseq2"),
    results = list(
      diff_result_df = data.frame(feature_id = c("a", "b"),
                                  effect_type = "log2FC"),
      diff_raw_df = data.frame(x = 1:2),
      diff_object = matrix(rnorm(2e5), ncol = 10)))
  # Read by run_integration() through side_analysis_id(), and by nothing
  # in enrichment -- which is exactly how a rebuilt "lightweight bundle"
  # loses it, and why this one is built by subtraction.
  b$provenance_binding <- list(analysis_id = "abc123")

  out <- bundle_for_worker(b)
  expect_null(out$results$diff_object)
  expect_null(out$results$diff_raw_df)
  expect_identical(out$results$diff_result_df, b$results$diff_result_df)
  expect_identical(out$params, b$params)
  expect_identical(out$input_info, b$input_info)
  expect_identical(out$provenance_binding, b$provenance_binding)
  expect_s3_class(out, "analysis_bundle")
  expect_lt(length(serialize(out, NULL)),
            length(serialize(b$results$diff_object, NULL)))
})

test_that("bundle_for_worker passes anything that is not a bundle through", {
  expect_null(bundle_for_worker(NULL))
  expect_identical(bundle_for_worker(list(a = 1)), list(a = 1))
})

test_that("project_for_worker keeps the data and drops the analyses", {
  m <- matrix(rnorm(30), nrow = 3,
              dimnames = list(paste0("f", 1:3), paste0("s", 1:10)))
  inp <- omicsCore::omics_input(
    m,
    data.frame(group = rep(c("a", "b"), 5), row.names = colnames(m)),
    data.frame(feature_id = rownames(m), row.names = rownames(m)),
    omics_type = "proteomics", assay_type = "normalized_intensity")
  p <- omicsCore::omics_project("p", experiments = list(prot = inp))
  p$bundles <- list(diff = list(big = matrix(rnorm(2e5), ncol = 10)))

  out <- project_for_worker(p)
  expect_null(out$bundles)
  # resolve_experiment_pair() reads this, so it has to survive
  expect_identical(names(out$experiments), "prot")
  expect_identical(out$experiments$prot$expr_mat, m)
  expect_s3_class(out, "omics_project")
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
