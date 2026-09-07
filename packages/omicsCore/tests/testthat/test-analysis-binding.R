# The data contract behind "which data is this result about, and may this
# callback still write it".
#
# Two questions that look like one and are not:
#
#   * a result loaded from disk is valid if its inputs are still the current
#     ones and its upstream analyses are still valid. Nothing about the
#     session it was computed in enters that.
#   * an in-flight callback may write if its session token is still current
#     AND the sources it captured are still current. Either alone lets a
#     stale result land: the token says the user has not moved on, the
#     sources say the data has not moved on, and the data can move under a
#     task without the user starting a new one.
#
# Each block below names the scenario it came from.

mk_input <- function(seed = 1, n_features = 10, symbols = NULL) {
  set.seed(seed)
  e <- matrix(rnorm(n_features * 4), nrow = n_features,
              dimnames = list(paste0("g", seq_len(n_features)),
                              paste0("s", seq_len(4))))
  md <- data.frame(group = rep(c("A", "B"), each = 2), row.names = colnames(e))
  fd <- data.frame(feature_id = rownames(e),
                   feature_symbol = symbols %||% rownames(e),
                   row.names = rownames(e), stringsAsFactors = FALSE)
  omics_input(e, md, fd, omics_type = "proteomics",
              assay_type = "normalized_intensity")
}

mk_project <- function(...) {
  exps <- list(...)
  if (!length(exps)) exps <- list(prot = mk_input())
  omics_project("p", experiments = exps)
}

mk_bundle <- function(project, tag, id = new_analysis_id(), upstream = character(0),
                      params = list(a = 1), name = "run_diff", inputs = NULL) {
  b <- new_analysis_bundle(analysis_name = name, results = list(x = 1))
  b$provenance_binding <- new_provenance_binding(
    analysis_id = id,
    project_id  = project$metadata$project_id,
    inputs      = inputs %||% stats::setNames(
      list(effective_input(project, tag)$version), tag),
    upstream_ids = upstream,
    params_signature = params_signature(params)
  )
  b
}

# --- identity ------------------------------------------------------------

test_that("ids do not consume the analysis random stream", {
  # sample() would both advance the RNG an analysis is about to use and
  # repeat itself under a fixed seed.
  set.seed(99); a <- new_analysis_id(); x1 <- runif(1)
  set.seed(99); b <- new_analysis_id(); x2 <- runif(1)
  expect_false(identical(a, b))
  expect_equal(x1, x2)
})

test_that("input_version answers to the whole input, not just the matrix", {
  base <- mk_input()
  expect_identical(input_version(base), input_version(mk_input()))

  # feature_symbol is the enrichment gene mapping and the integration join
  # key, so changing it changes the analysis input.
  resym <- mk_input(symbols = paste0("SYM", 1:10))
  expect_false(identical(input_version(base), input_version(resym)))

  remeta <- base; remeta$meta_df$batch <- rep(1:2, 2)
  expect_false(identical(input_version(base), input_version(remeta)))

  retype <- base; retype$omics_type <- "rnaseq"
  expect_false(identical(input_version(base), input_version(retype)))
})

# --- C2: what makes a saved result valid ---------------------------------

test_that("a binding from another project is not current here", {
  # Two projects on identical data is an ordinary workflow. Comparing input
  # versions alone accepted the other project's result.
  p1 <- mk_project(); p2 <- mk_project()
  foreign <- mk_bundle(p2, "prot")
  expect_false(result_inputs_current(foreign$provenance_binding, p1))
  expect_true(result_inputs_current(mk_bundle(p1, "prot")$provenance_binding, p1))
})

test_that("an empty input map is not a verified binding", {
  # all(logical(0)) is TRUE, so a binding that names no inputs passed.
  p <- mk_project()
  pb <- mk_bundle(p, "prot")$provenance_binding
  pb$inputs <- list()
  expect_false(validate_provenance_binding(pb)$ok)
  expect_false(result_inputs_current(pb, p))
})

test_that("structurally broken bindings are reported, not trusted", {
  p <- mk_project()
  good <- mk_bundle(p, "prot")$provenance_binding
  cases <- list(
    no_id       = within(good, analysis_id <- NULL),
    blank_id    = within(good, analysis_id <- ""),
    no_project  = within(good, project_id <- NULL),
    unnamed     = within(good, inputs <- list("v1")),
    na_version  = within(good, inputs <- list(prot = NA_character_)),
    bad_upstream = within(good, upstream_ids <- list(1L))
  )
  for (nm in names(cases)) {
    res <- validate_provenance_binding(cases[[nm]])
    expect_false(isTRUE(res$ok), info = nm)
    expect_true(nzchar(res$reason), info = nm)
  }
  expect_true(validate_provenance_binding(good)$ok)
})

test_that("a bundle with no binding at all does not break the algorithm", {
  # An .omp written before the contract existed. It must load, and it must
  # not be adopted into the current project by default.
  p <- mk_project()
  legacy <- new_analysis_bundle(analysis_name = "run_diff", results = list(x = 1))
  p <- register_analysis(p, legacy, id = "legacy-1")
  expect_silent(ids <- valid_result_ids(p))
  expect_false("legacy-1" %in% ids)
})

test_that("a saved result whose inputs are unchanged stays valid", {
  p <- mk_project()
  p <- register_analysis(p, mk_bundle(p, "prot", id = "d1"))
  expect_true("d1" %in% valid_result_ids(p))

  # Round-tripping through .omp must not change the answer -- activation is
  # a session counter and has no business in this decision.
  skip_if_not_installed("qs2")
  f <- withr::local_tempfile(fileext = ".omp")
  save_project(p, f)
  expect_true("d1" %in% valid_result_ids(load_project(f)))
})

# --- C1: what makes an in-flight callback admissible ---------------------

ctx_for <- function(project, activation = 1L, seq = list(diff = 1L),
                    session = "sess-1") {
  list(session_id = session, project_id = project$metadata$project_id,
       activation = activation, request_seq = seq)
}

test_that("a callback is refused when the data moved under it", {
  # The scenario the review reproduced: same project, same session, same
  # request number, but the layer was replaced while the task ran.
  p <- mk_project()
  tok <- new_session_token(session_id = "sess-1",
                           project_id = p$metadata$project_id,
                           activation = 1L, request_seq = 1L, view = "diff",
                           params_signature = params_signature(list(a = 1)))
  pb <- mk_bundle(p, "prot")$provenance_binding
  expect_true(callback_admissible(tok, ctx_for(p), pb, p)$ok)

  moved <- p
  moved$experiments$prot <- mk_input(seed = 2)
  res <- callback_admissible(tok, ctx_for(moved), pb, moved)
  expect_false(res$ok)
  expect_match(res$reason, "input")
})

test_that("a callback is refused when a newer request has been made", {
  p <- mk_project()
  tok <- new_session_token("sess-1", p$metadata$project_id, 1L, 1L, "diff",
                           params_signature(list(a = 1)))
  pb <- mk_bundle(p, "prot")$provenance_binding
  expect_false(callback_admissible(tok, ctx_for(p, seq = list(diff = 2L)), pb, p)$ok)
})

test_that("request numbers do not cross between views", {
  p <- mk_project()
  tok <- new_session_token("sess-1", p$metadata$project_id, 1L, 1L, "diff",
                           params_signature(list(a = 1)))
  # Enrichment has moved on; the diff task is still the current diff task.
  ctx <- ctx_for(p, seq = list(diff = 1L, enrich = 7L))
  expect_true(callback_admissible(tok, ctx, mk_bundle(p, "prot")$provenance_binding, p)$ok)
})

test_that("a missing parameter signature is refused rather than matched", {
  # Two NULLs are identical(), so a token built without the field passed the
  # check that was supposed to compare it.
  p <- mk_project()
  tok <- new_session_token("sess-1", p$metadata$project_id, 1L, 1L, "diff",
                           params_signature(list(a = 1)))
  tok$params_signature <- NULL
  pb <- mk_bundle(p, "prot")$provenance_binding
  res <- callback_admissible(tok, ctx_for(p), pb, p)
  expect_false(res$ok)
  expect_match(res$reason, "signature")
})

test_that("both callback outcomes are held to the same test", {
  # A failed task writing its error into a view the user has already moved
  # on from is the same defect as a stale success.
  p <- mk_project()
  tok <- new_session_token("sess-1", p$metadata$project_id, 1L, 1L, "diff",
                           params_signature(list(a = 1)))
  pb <- mk_bundle(p, "prot")$provenance_binding
  stale_ctx <- ctx_for(p, activation = 2L)
  expect_false(callback_admissible(tok, stale_ctx, pb, p, outcome = "success")$ok)
  expect_false(callback_admissible(tok, stale_ctx, pb, p, outcome = "error")$ok)
})

# --- C3: the dependency graph must be able to find its nodes -------------

integration_project <- function(register_secondary = TRUE) {
  p <- mk_project(prot = mk_input(1), rna = mk_input(2))
  d_a <- mk_bundle(p, "prot", id = "diff-a")
  d_b <- mk_bundle(p, "rna",  id = "diff-b")
  integ <- mk_bundle(p, "prot", id = "integ", name = "run_integration",
                     upstream = c("diff-a", "diff-b"),
                     inputs = list(prot = effective_input(p, "prot")$version,
                                   rna  = effective_input(p, "rna")$version))
  p <- register_analysis(p, d_a)
  if (register_secondary) p <- register_analysis(p, d_b)
  p <- register_analysis(p, integ)
  p
}

test_that("an integration is valid the moment it completes", {
  # The secondary diff used to exist only inside the worker, so an
  # integration that honestly recorded both upstream ids was immediately
  # judged invalid for referencing one that could not be found.
  expect_true("integ" %in% valid_result_ids(integration_project(TRUE)))
  expect_false("integ" %in% valid_result_ids(integration_project(FALSE)))
})

test_that("invalidation is transitive and order-independent", {
  p <- integration_project()
  p <- register_analysis(p, mk_bundle(p, "rna", id = "enrich-b",
                                      upstream = "diff-b", name = "run_enrichment"))
  expect_setequal(valid_result_ids(p),
                  c("diff-a", "diff-b", "integ", "enrich-b"))

  # Only the second layer changes. Its diff, the enrichment on that diff and
  # the integration that reads it all stop being about the project; the
  # first layer's diff is untouched.
  moved <- p
  moved$experiments$rna <- mk_input(seed = 42)
  expect_setequal(valid_result_ids(moved), "diff-a")

  # A single pass over the list would keep whichever dependant happened to
  # be visited before its upstream was removed.
  shuffled <- moved
  shuffled$analyses <- shuffled$analyses[rev(names(shuffled$analyses))]
  expect_setequal(valid_result_ids(shuffled), "diff-a")
})

test_that("retiring keeps the graph consistent and bounded", {
  p <- integration_project()
  p$experiments$rna <- mk_input(seed = 42)
  p <- retire_invalid(p)
  expect_setequal(names(p$analyses), "diff-a")
  expect_true(all(c("diff-b", "integ") %in% names(p$retired)))
})

test_that("history is compressed and capped", {
  # 21 rounds of re-analysis must not put 21 fitted models into every
  # autosave.
  p <- mk_project()
  for (i in 1:21) {
    b <- mk_bundle(p, "prot", id = sprintf("d%02d", i))
    b$results$diff_object <- runif(5000)     # a stand-in for a fitted model
    b$model_object <- runif(5000)
    p <- register_analysis(p, b)
    p$experiments$prot <- mk_input(seed = 100 + i)   # invalidates the last one
    p <- retire_invalid(p, max_history = 20L)
  }
  expect_lte(length(p$retired), 20L)
  expect_equal(anyDuplicated(names(p$retired)), 0L)
  for (h in p$retired) {
    expect_null(h$results$diff_object)
    expect_null(h$model_object)
    expect_true(nzchar(h$provenance_binding$analysis_id))
    expect_false(is.null(h$retired_at))
  }
})

test_that("retiring never drops a node a valid result still references", {
  # The history cap is a display and storage limit. Applying it to the
  # dependency graph would invalidate the integration it was trying to keep.
  p <- integration_project()
  p <- retire_invalid(p, max_history = 0L)
  expect_true("integ" %in% valid_result_ids(p))
  expect_true(all(c("diff-a", "diff-b") %in% names(p$analyses)))
})

# --- per-layer preprocessing --------------------------------------------

test_that("each layer keeps its own applied preprocessing", {
  p <- mk_project(prot = mk_input(1), rna = mk_input(2))
  cleaned_p <- mk_input(1, n_features = 8)
  cleaned_r <- mk_input(2, n_features = 6)
  p <- set_preprocessing(p, "prot", cleaned_p, params = list(missing = 0.3))
  p <- set_preprocessing(p, "rna",  cleaned_r, params = list(missing = 0.5))

  expect_equal(nrow(effective_input(p, "prot")$input$expr_mat), 8L)
  expect_equal(nrow(effective_input(p, "rna")$input$expr_mat), 6L)
  expect_equal(effective_input(p, "prot")$source, "qc_cleaned")
})

test_that("preprocessing built on a replaced layer is not served", {
  p <- mk_project()
  p <- set_preprocessing(p, "prot", mk_input(1, n_features = 8),
                         params = list(missing = 0.3))
  expect_equal(effective_input(p, "prot")$source, "qc_cleaned")

  p$experiments$prot <- mk_input(seed = 5)
  eff <- effective_input(p, "prot")
  expect_equal(eff$source, "imported")
  expect_true(eff$stale)
})

test_that("QC reads the imported baseline, not the last cleaned layer", {
  # Re-running QC on already-filtered data cannot return features the
  # previous run removed, so loosening a threshold would do nothing.
  p <- mk_project()
  p <- set_preprocessing(p, "prot", mk_input(1, n_features = 4),
                         params = list(missing = 0.1))
  expect_equal(nrow(baseline_input(p, "prot")$input$expr_mat), 10L)
  expect_equal(nrow(effective_input(p, "prot")$input$expr_mat), 4L)
})

test_that("a layer's version is cached rather than recomputed per query", {
  p <- mk_project()
  calls <- 0L
  testthat::local_mocked_bindings(
    input_version = function(...) { calls <<- calls + 1L; "v" }
  )
  for (i in 1:3) effective_input(p, "prot")
  expect_equal(calls, 0L)
})

test_that("an analysis with no parameters still has a signature", {
  # order(NULL) is an error rather than an empty ordering, and an analysis
  # with nothing to configure is ordinary.
  expect_true(nzchar(params_signature(list())))
  expect_true(nzchar(params_signature(NULL)))
  expect_identical(params_signature(list()), params_signature(NULL))
  expect_identical(params_signature(list(b = 2, a = 1)),
                   params_signature(list(a = 1, b = 2)))
})
