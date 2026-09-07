# A result that comes back after the project has moved.
#
# The gate used to compare a request number and nothing else, so it saw
# neither a change of project nor a change of data. Worse, the record a
# result was stamped with was assembled from two moments: the input
# versions from launch, the project identity from arrival. That produces a
# record which looks internally consistent -- versions that match, a
# project that exists -- while describing a computation that happened
# somewhere else. On two projects built from the same file, a late result
# from A landed in B and passed the ownership check.
#
# These drive the real modules and hold the real scheduler still, rather
# than injecting a gate that always refuses: a gate that refuses everything
# proves nothing about the one the module actually builds.

lc_input <- function(seed = 1, n_features = 20) {
  set.seed(seed)
  n <- 8
  e <- matrix(rnorm(n_features * n, 10, 1), nrow = n_features,
              dimnames = list(paste0("g", seq_len(n_features)),
                              paste0("s", seq_len(n))))
  e[1:5, 5:8] <- e[1:5, 5:8] + 2
  md <- data.frame(group = rep(c("A", "B"), each = 4),
                   row.names = colnames(e), stringsAsFactors = FALSE)
  fd <- data.frame(feature_id = rownames(e), feature_symbol = rownames(e),
                   row.names = rownames(e), stringsAsFactors = FALSE)
  omicsCore::omics_input(e, md, fd, omics_type = "proteomics",
                         assay_type = "normalized_intensity")
}

# Two projects, same data. This is an ordinary workflow -- one dataset,
# two contrasts, two write-ups -- and it is what makes identity rather
# than content the only thing separating them.
lc_projects <- function() {
  inp <- lc_input()
  list(a = omicsCore::omics_project("Project A", list(prot = inp)),
       b = omicsCore::omics_project("Project B", list(prot = inp)))
}

# Hold the result back so the state can be changed underneath it, then
# release. Everything else -- the launch record, the gate, the callbacks --
# is the module's own.
held_scheduler <- function() {
  pending <- NULL
  list(
    future = function(expr, ...) {
      pending <<- function() force(expr)
      structure(list(), class = "held_future")
    },
    release = function() if (!is.null(pending)) pending()
  )
}

test_that("a result that comes back into another project is refused", {
  skip_if_not_installed("limma")
  p <- lc_projects()
  proj <- shiny::reactiveVal(p$a)
  act <- shiny::reactiveVal(1L)

  shiny::testServer(
    diff_view_server,
    args = list(current_project = proj, activation = act),
    {
      session$setInputs(layer = "prot", method = "limma", group_col = "group",
                        control = "A", case = "B")
      session$flushReact()
      # The run starts against A.
      launched <- capture_launch(proj(), "diff", tags = "prot",
                                 params = list(x = 1), activation = act())
      expect_equal(launched$project_id, p$a$metadata$project_id)

      # The user opens B. Same data, different project, so every input
      # version still matches and only the identity separates them.
      proj(p$b)
      act(2L)
      session$flushReact()

      gate <- launch_admissible(
        launched,
        function() list(project = proj(), activation = act(),
                        request_seq = 1L),
        1L)
      for (outcome in c("success", "error")) {
        verdict <- gate(outcome)
        expect_false(verdict$ok, info = outcome)
        expect_true(nzchar(verdict$reason), info = outcome)
      }
    })
})

test_that("a result is stamped with where it started, not where it landed", {
  p <- lc_projects()
  launched <- capture_launch(p$a, "diff", tags = "prot",
                             params = list(x = 1), activation = 1L)
  b <- omicsCore::new_analysis_bundle("run_diff", params = list(x = 1),
                                      results = list(x = 1))
  stamped <- stamp_from_launch(b, launched)

  expect_equal(stamped$provenance_binding$project_id,
               p$a$metadata$project_id)
  # And in B -- identical data, so the versions match -- it is still not
  # B's result. That check is the whole reason project_id is in the record.
  expect_false(omicsCore::result_inputs_current(stamped$provenance_binding, p$b))
  expect_true(omicsCore::result_inputs_current(stamped$provenance_binding, p$a))
})

test_that("a result is refused when the data changed under it", {
  p <- lc_projects()
  launched <- capture_launch(p$a, "diff", tags = "prot",
                             params = list(x = 1), activation = 1L)
  moved <- p$a
  moved <- omicsCore::set_preprocessing(moved, "prot", lc_input(n_features = 12),
                                        params = list(missing_threshold = 0.1))
  gate <- launch_admissible(
    launched,
    function() list(project = moved, activation = 1L, request_seq = 1L), 1L)
  expect_false(gate("success")$ok)
  expect_match(gate("success")$reason, "data changed")
})

test_that("a context that cannot be read is a refusal, not a pass", {
  # The earlier gate defaulted to admitting when it could not read the
  # current request number -- the wrong way round for a check whose whole
  # purpose is the case where things have moved.
  launched <- capture_launch(lc_projects()$a, "diff", tags = "prot",
                             params = list(), activation = 1L)
  gate <- launch_admissible(launched, function() stop("session gone"), 1L)
  expect_false(gate("success")$ok)
  expect_false(gate("error")$ok)

  gate2 <- launch_admissible(launched, function() list(project = NULL), 1L)
  expect_false(gate2("success")$ok)
})

test_that("reopening the same project refuses an outstanding result", {
  # Same project id, same data: only the act of reopening has happened,
  # and that is exactly the case a project comparison cannot see.
  p <- lc_projects()
  launched <- capture_launch(p$a, "diff", tags = "prot",
                             params = list(), activation = 1L)
  gate <- launch_admissible(
    launched,
    function() list(project = p$a, activation = 2L, request_seq = 1L), 1L)
  expect_false(gate("success")$ok)
  expect_match(gate("success")$reason, "reopened")
})

test_that("an unchanged project admits the result it was launched with", {
  # The negative cases above are only meaningful if the ordinary path
  # still works.
  p <- lc_projects()
  launched <- capture_launch(p$a, "diff", tags = "prot",
                             params = list(), activation = 1L)
  gate <- launch_admissible(
    launched,
    function() list(project = p$a, activation = 1L, request_seq = 4L), 4L)
  expect_true(gate("success")$ok)
  expect_true(gate("error")$ok)
})

test_that("an upstream that has been retired refuses its dependant", {
  p <- lc_projects()$a
  upstream <- stamp_from_launch(
    omicsCore::new_analysis_bundle("run_diff", results = list(x = 1)),
    capture_launch(p, "diff", tags = "prot", params = list(), activation = 1L))
  launched <- capture_launch(
    p, "enrich", tags = "prot",
    upstream = upstream$provenance_binding$analysis_id,
    params = list(), activation = 1L)

  with_up <- omicsCore::register_analysis(p, upstream)
  gate_ok <- launch_admissible(
    launched,
    function() list(project = with_up, activation = 1L, request_seq = 1L), 1L)
  expect_true(gate_ok("success")$ok)

  gate_gone <- launch_admissible(
    launched,
    function() list(project = p, activation = 1L, request_seq = 1L), 1L)
  expect_false(gate_gone("success")$ok)
  expect_match(gate_gone("success")$reason, "upstream")
})

test_that("the demo still computes, and is still bound by the request number", {
  # No project means no ownership question -- there is nowhere to
  # misattribute a result to, and stamp_from_launch() leaves it unstamped
  # so the write-back will not file it. Refusing outright here broke the
  # built-in demo, which runs the real analyses on fixtures.
  gate <- launch_admissible(NULL, function() list(project = NULL,
                                                  activation = 1L,
                                                  request_seq = 2L), 2L)
  expect_true(gate("success")$ok)
  stale <- launch_admissible(NULL, function() list(project = NULL,
                                                   activation = 1L,
                                                   request_seq = 3L), 2L)
  expect_false(stale("success")$ok)
})

test_that("an unattributed result is not filed into any project", {
  b <- omicsCore::new_analysis_bundle("run_diff", results = list(x = 1))
  expect_identical(stamp_from_launch(b, NULL), b)
  expect_null(stamp_from_launch(b, NULL)$provenance_binding)
})
