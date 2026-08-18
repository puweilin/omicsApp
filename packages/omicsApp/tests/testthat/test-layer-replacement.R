# Replacing an omics layer must take the analyses computed on the old
# data with it. Clearing `project$bundles` alone does not do that: the
# analysis views own what they computed, and the bundle-attach observer
# in app_server() puts it straight back. These tests cover both halves --
# the fingerprint that decides whether a re-import is a real change, and
# the signal that tells the views to let go.

fake_input <- function(fingerprint = NULL, n_features = 3L) {
  mat <- matrix(as.numeric(seq_len(n_features * 4L)), nrow = n_features,
                dimnames = list(paste0("g", seq_len(n_features)),
                                paste0("s", 1:4)))
  meta <- data.frame(group = c("A", "A", "B", "B"),
                     row.names = paste0("s", 1:4))
  feat <- data.frame(feature_id = paste0("g", seq_len(n_features)))
  omicsCore::omics_input(mat, meta, feat, omics_type = "proteomics",
                         source_fingerprint = fingerprint)
}

project_holding <- function(input, bundles = list()) {
  proj <- omicsCore::omics_project(name = "P",
                                   experiments = list(proteomics = input))
  proj$bundles <- bundles
  proj
}

# ---- fingerprints -----------------------------------------------------

test_that("the same file yields the same fingerprint", {
  path <- tempfile(); on.exit(unlink(path), add = TRUE)
  writeLines("some,content", path)
  expect_identical(input_fingerprint(path, "proteomics", "intensity"),
                   input_fingerprint(path, "proteomics", "intensity"))
})

test_that("different content yields a different fingerprint", {
  a <- tempfile(); b <- tempfile()
  on.exit(unlink(c(a, b)), add = TRUE)
  writeLines("content A", a)
  writeLines("content B", b)
  expect_false(identical(input_fingerprint(a, "proteomics", "intensity"),
                         input_fingerprint(b, "proteomics", "intensity")))
})

test_that("the parse settings are part of the fingerprint", {
  path <- tempfile(); on.exit(unlink(path), add = TRUE)
  writeLines("same,bytes", path)
  # The same workbook read as proteomics and as RNA-seq is two different
  # inputs, so re-importing across that switch is a real change.
  expect_false(identical(input_fingerprint(path, "proteomics", "intensity"),
                         input_fingerprint(path, "rnaseq", "raw_count")))
})

test_that("an unreadable path has no fingerprint", {
  expect_null(input_fingerprint(tempfile(), "proteomics", "intensity"))
})

test_that("a match needs a fingerprint on both sides", {
  expect_true(fingerprints_match(fake_input("fp-A"), fake_input("fp-A")))
  expect_false(fingerprints_match(fake_input("fp-A"), fake_input("fp-B")))
  # Absence is not evidence of sameness: an input built straight from
  # matrices carries no fingerprint, and assuming "unchanged" there
  # would silently keep results belonging to other data.
  expect_false(fingerprints_match(fake_input(NULL), fake_input(NULL)))
  expect_false(fingerprints_match(fake_input("fp-A"), fake_input(NULL)))
})

# ---- the Confirm gate -------------------------------------------------

test_that("a first import commits without asking", {
  shiny::testServer(
    import_view_server,
    args = list(current_project = shiny::reactiveVal(NULL)),
    {
      parsed(list(input = fake_input("fp-A"), report = NULL))
      session$setInputs(confirm = 1)
      expect_false(is.null(confirmed_input()))
    }
  )
})

test_that("re-picking the same file commits nothing", {
  proj <- project_holding(fake_input("fp-A"))
  shiny::testServer(
    import_view_server,
    args = list(current_project = shiny::reactiveVal(proj)),
    {
      parsed(list(input = fake_input("fp-A"), report = NULL))
      session$setInputs(confirm = 1)
      # A mis-click must not cost the user their analyses, so this is a
      # no-op rather than a replace.
      expect_null(confirmed_input())
    }
  )
})

test_that("new data does not commit until the replace is confirmed", {
  proj <- project_holding(fake_input("fp-A"),
                          bundles = list(diff = "OLD_DIFF"))
  shiny::testServer(
    import_view_server,
    args = list(current_project = shiny::reactiveVal(proj)),
    {
      parsed(list(input = fake_input("fp-B"), report = NULL))
      session$setInputs(confirm = 1)
      # The dialog is up; nothing has been handed to app_server yet.
      expect_null(confirmed_input())

      session$setInputs(confirm_replace = 1)
      committed <- confirmed_input()
      expect_false(is.null(committed))
      expect_identical(committed$source_fingerprint, "fp-B")
    }
  )
})

test_that("the replace dialog names the analyses that will be lost", {
  proj <- project_holding(fake_input("fp-A"),
                          bundles = list(qc = 1, diff = 2))
  html <- render_html(replace_layer_modal(function(x) x, "proteomics", proj))
  expect_match(html, "Quality control", fixed = TRUE)
  expect_match(html, "Differential analysis", fixed = TRUE)
  expect_match(html, "Replace and clear", fixed = TRUE)
})

# ---- the invalidation signal ------------------------------------------
#
# Each analysis view holds its own result. NULL is every module's
# start-up state, so a bump only rewinds it there.

expect_clears_on_invalidate <- function(server, args, bundle_name) {
  gen <- shiny::reactiveVal(0L)
  shiny::testServer(
    server,
    args = c(args, list(invalidate = gen)),
    {
      # `ignoreInit` skips the observer's *first* execution, and under
      # testServer that has not happened yet. A real session flushes at
      # start-up long before anyone imports a file; flushing here puts
      # the module in the same state, so the bump below is a genuine
      # change rather than the one ignoreInit swallows.
      session$flushReact()

      held <- get(bundle_name, envir = environment())
      held("PRETEND_BUNDLE")
      expect_equal(held(), "PRETEND_BUNDLE")
      gen(1L)
      session$flushReact()
      expect_null(held())
    }
  )
}

test_that("replacing the layer clears the QC result", {
  expect_clears_on_invalidate(
    qc_view_server, list(current_project = shiny::reactiveVal(NULL)),
    "last_bundle")
})

test_that("replacing the layer clears the differential result", {
  expect_clears_on_invalidate(
    diff_view_server, list(current_project = shiny::reactiveVal(NULL)),
    "diff_bundle")
})

test_that("replacing the layer clears the enrichment result", {
  expect_clears_on_invalidate(
    enrich_view_server, list(diff_bundle = shiny::reactiveVal(NULL)),
    "enrich_bundle")
})

test_that("replacing the layer clears the integration result", {
  expect_clears_on_invalidate(
    integration_view_server, list(current_project = shiny::reactiveVal(NULL)),
    "integration_bundle")
})

test_that("a view keeps its result when nothing was replaced", {
  gen <- shiny::reactiveVal(0L)
  shiny::testServer(
    diff_view_server,
    args = list(current_project = shiny::reactiveVal(NULL), invalidate = gen),
    {
      session$flushReact()
      diff_bundle("PRETEND_BUNDLE")
      session$flushReact()
      # No bump: importing a second, *different* omics type extends the
      # project rather than replacing a layer, and must leave existing
      # results alone.
      expect_equal(diff_bundle(), "PRETEND_BUNDLE")
    }
  )
})
