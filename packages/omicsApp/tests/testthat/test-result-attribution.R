# Which project a result belongs to, and where it stops being one.
#
# Two failures sit behind these. Results were dropped from the project only
# when a view had held them first, so a project restored from disk with a
# saved differential and an empty view kept that differential through every
# change of data -- while the report, the export and the autosave went on
# reading it. And "the same project" was decided by comparing data, so
# opening an earlier save of the same project changed nothing the views could
# see, and what they held was written over the file that had just been opened.

# Local copies: testthat gives each test file its own environment, so the
# helpers in the acceptance file are not in scope here.
ra_skip_no_xlsx <- function() {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("readxl")
  skip_if_not_installed("withr")
}

ra_store <- function(env = parent.frame()) {
  dir <- file.path(tempfile("attribution-"), "store")
  dir.create(dir, recursive = TRUE)
  withr::local_envvar(OMICSAPP_DATA_DIR = dir, .local_envir = env)
  withr::defer(unlink(dirname(dir), recursive = TRUE), envir = env)
  dir
}

ra_workbook <- function(env = parent.frame()) {
  path <- tempfile(fileext = ".xlsx")
  withr::defer(unlink(path), envir = env)
  write_tiny_omics_xlsx(path, n_features = 40L, n_samples = 8L, seed = 3)
  path
}

ra_upload <- function(session, path) {
  suppressWarnings(session$setInputs(
    `import-omics_type` = "proteomics",
    `import-file` = list(datapath = path, name = basename(path),
                         size = file.info(path)$size)
  ))
  session$setInputs(`import-confirm` = 1)
}

ra_input <- function(seed = 1, n = 8) {
  set.seed(seed)
  e <- matrix(rnorm(n * 4), nrow = n,
              dimnames = list(paste0("g", seq_len(n)), paste0("s", 1:4)))
  md <- data.frame(group = rep(c("A", "B"), each = 2), row.names = colnames(e))
  fd <- data.frame(feature_id = rownames(e), feature_symbol = rownames(e),
                   row.names = rownames(e), stringsAsFactors = FALSE)
  omicsCore::omics_input(e, md, fd, omics_type = "proteomics",
                         assay_type = "normalized_intensity")
}

ra_project <- function(name = "p", seed = 1) {
  omicsCore::omics_project(name, list(prot = ra_input(seed)))
}

# A differential that says which data it came from, as the views now stamp.
ra_result <- function(project, tag = "prot", params = list(comparison = "B_vs_A"),
                      view = "diff") {
  b <- omicsCore::new_analysis_bundle(
    "run_diff", params = params,
    results = list(diff_result_df = data.frame(feature_id = "g1", effect = 1)))
  bind_result(b, project, view, tags = tag,
              versions = omicsCore::effective_input(project, tag)$version)
}

# --- the project's own record decides ------------------------------------

test_that("a restored project drops a saved result when its data changes", {
  # The view is empty -- nothing has been opened -- so the old rule had
  # nothing to compare against and left the result in place.
  p <- ra_project()
  p <- omicsCore::register_analysis(p, ra_result(p))
  p$bundles <- list(diff = ra_result(p))
  expect_length(p$bundles, 1L)

  # QC applied: the effective input is no longer what the result ran on.
  cleaned <- ra_input(1, n = 5)
  p <- omicsCore::set_preprocessing(p, "prot", cleaned,
                                    params = list(missing_threshold = 0.1))
  p <- omicsCore::retire_invalid(p, reason = "QC applied")

  expect_length(p$bundles, 0L)
  expect_length(p$analyses, 0L)
  expect_length(p$retired, 1L)
  expect_match(p$retired[[1]]$retired_because, "QC applied")
})

test_that("a result whose data is unchanged survives retirement", {
  p <- ra_project()
  p <- omicsCore::register_analysis(p, ra_result(p))
  p$bundles <- list(diff = p$analyses[[1]])
  p <- omicsCore::retire_invalid(p)
  expect_length(p$bundles, 1L)
  expect_length(p$analyses, 1L)
})

test_that("retirement reaches the autosave, not just the page", {
  skip_if_not_installed("qs2")
  dir <- withr::local_tempdir()
  withr::local_envvar(OMICSAPP_DATA_DIR = dir)

  p <- ra_project()
  p <- omicsCore::register_analysis(p, ra_result(p))
  p$bundles <- list(diff = p$analyses[[1]])
  expect_true(store_autosave(p, dir = dir)$ok)
  expect_length(store_read_autosave(dir)$bundles, 1L)

  p <- omicsCore::set_preprocessing(p, "prot", ra_input(1, n = 5),
                                    params = list(missing_threshold = 0.1))
  p <- omicsCore::retire_invalid(p)
  expect_true(store_autosave(p, dir = dir)$ok)

  # A recovery copy holding results computed on data the project no longer
  # has is the stale-result bug, persisted.
  back <- store_read_autosave(dir)
  expect_length(back$bundles, 0L)
  expect_length(back$retired, 1L)
})

# --- attribution ---------------------------------------------------------

test_that("a result from another project is never adopted", {
  # Two projects on identical data is an ordinary workflow, so the versions
  # match and only the project id separates them.
  a <- ra_project("A")
  b <- ra_project("B")
  from_a <- ra_result(a)
  expect_false(omicsCore::result_inputs_current(from_a$provenance_binding, b))

  b <- omicsCore::register_analysis(b, from_a)
  expect_length(omicsCore::valid_result_ids(b), 0L)
})

test_that("a result with no attribution is not filled in with the current one", {
  # `%||% pid` looked like a fallback. It meant an unattributed result was
  # accepted as this project's, which is the mislabelling being prevented.
  p <- ra_project()
  bare <- omicsCore::new_analysis_bundle("run_diff", results = list(x = 1))
  expect_null(bare$provenance_binding)
  p <- omicsCore::register_analysis(p, bare, id = "bare-1")
  expect_false("bare-1" %in% omicsCore::valid_result_ids(p))
})

test_that("bind_result refuses to stamp what it cannot describe", {
  p <- ra_project()
  b <- omicsCore::new_analysis_bundle("run_diff", results = list(x = 1))
  expect_null(bind_result(b, NULL, "diff", "prot", "v")$provenance_binding)
  expect_null(bind_result(b, p, "diff", character(0), character(0))$provenance_binding)
  # Tags and versions have to line up, or the map would be a guess.
  expect_null(bind_result(b, p, "diff", c("a", "b"), "v1")$provenance_binding)
})

test_that("the same analysis on the same data keeps the same id", {
  # The id is derived, not drawn. Filing a result the view is still holding
  # has to be idempotent, or the observer that files it changes the project
  # on every flush and re-triggers itself.
  p <- ra_project()
  expect_identical(ra_result(p)$provenance_binding$analysis_id,
                   ra_result(p)$provenance_binding$analysis_id)
  # Different settings are a different analysis.
  expect_false(identical(
    ra_result(p, params = list(comparison = "B_vs_A"))$provenance_binding$analysis_id,
    ra_result(p, params = list(comparison = "A_vs_B"))$provenance_binding$analysis_id))
})

# --- transitive ----------------------------------------------------------

test_that("enrichment goes when the differential under it goes", {
  p <- ra_project()
  d <- ra_result(p)
  e <- omicsCore::new_analysis_bundle("run_enrichment",
                                      results = list(enrich_result_df = data.frame(x = 1)))
  e <- bind_result(e, p, "enrich", tags = "prot",
                   versions = omicsCore::effective_input(p, "prot")$version,
                   upstream = d$provenance_binding$analysis_id)
  p <- omicsCore::register_analysis(p, d)
  p <- omicsCore::register_analysis(p, e)
  expect_length(omicsCore::valid_result_ids(p), 2L)

  p$experiments$prot <- ra_input(seed = 99)
  p$experiments$prot$version <- omicsCore::input_version(p$experiments$prot)
  expect_length(omicsCore::valid_result_ids(p), 0L)
})

# --- opening a project -----------------------------------------------------

test_that("opening an earlier save of the same project is an activation", {
  ra_skip_no_xlsx()
  store <- ra_store()
  xlsx <- ra_workbook()
  shiny::testServer(app_server, {
    ra_upload(session, xlsx)
    first <- activation()
    session$setInputs(`project-save_name` = "run one",
                      `project-save_project` = 1)
    session$flushReact()

    # Same project id, same data. Only the act of opening has happened, and
    # comparing ids or fingerprints cannot see it.
    session$setInputs(`project-saved_pick` = "run_one",
                      `project-open_project` = 1)
    session$flushReact()
    expect_gt(activation(), first)
  })
})

test_that("restoring the autosave is an activation too", {
  ra_skip_no_xlsx()
  store <- ra_store()
  xlsx <- ra_workbook()
  shiny::testServer(app_server, {
    ra_upload(session, xlsx)
    before <- activation()
    session$setInputs(`project-restore_autosave` = 1)
    session$flushReact()
    expect_gt(activation(), before)
  })
})

# --- opening a project puts its results back on screen --------------------

# Clearing the views on activation is right -- what they held belongs to
# the last project. Stopping there left the page empty in front of a
# project whose differential was on disk and valid, so the only move was to
# re-run it and hope the settings matched. Re-running is not restoring.
test_that("opening a project shows the differential it was saved with", {
  skip_if_not_installed("qs2")
  p <- ra_project()
  saved <- ra_result(p, params = list(method = "limma", group_col = "group",
                                      control_group = "A", case_group = "B",
                                      comparison = "B_vs_A"))
  p <- omicsCore::register_analysis(p, saved)
  p$bundles <- list(diff = saved)

  proj <- shiny::reactiveVal(NULL)
  act <- shiny::reactiveVal(0L)
  shiny::testServer(diff_view_server,
                    args = list(current_project = proj, activation = act), {
    session$flushReact()
    # The project arrives, as Open does it.
    proj(p)
    act(1L)
    session$flushReact()
    got <- session$getReturned()$bundle()
    expect_false(is.null(got), info = "the saved result was not restored")
    expect_equal(got$params$comparison, "B_vs_A")
  })
})

test_that("a saved result the data has moved past is not restored", {
  # The restore is not "show whatever is in the file": a result whose
  # inputs no longer match is exactly what retirement is for.
  skip_if_not_installed("qs2")
  p <- ra_project()
  saved <- ra_result(p)
  p$bundles <- list(diff = saved)
  p <- omicsCore::set_preprocessing(p, "prot", ra_input(seed = 77),
                                    params = list(missing_threshold = 0.1))

  proj <- shiny::reactiveVal(NULL)
  act <- shiny::reactiveVal(0L)
  shiny::testServer(diff_view_server,
                    args = list(current_project = proj, activation = act), {
    proj(p); act(1L)
    session$flushReact()
    expect_null(session$getReturned()$bundle())
  })
})

test_that("an unattributed saved result is not restored either", {
  p <- ra_project()
  p$bundles <- list(diff = omicsCore::new_analysis_bundle(
    "run_diff", results = list(diff_result_df = data.frame(x = 1))))
  proj <- shiny::reactiveVal(NULL)
  act <- shiny::reactiveVal(0L)
  shiny::testServer(diff_view_server,
                    args = list(current_project = proj, activation = act), {
    proj(p); act(1L)
    session$flushReact()
    expect_null(session$getReturned()$bundle())
  })
})
