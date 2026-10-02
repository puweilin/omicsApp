# App defects found in the 2026-10 audit.

aa_store <- function(env = parent.frame()) {
  dir <- withr::local_tempdir(.local_envir = env)
  withr::local_envvar(OMICSAPP_DATA_DIR = dir, .local_envir = env)
  dir
}

test_that("a slug from the browser cannot reach outside the project store", {
  dir <- aa_store()
  victim <- file.path(dirname(dir), "victim.omp")
  writeLines("keep me", victim)
  on.exit(unlink(victim), add = TRUE)
  slug <- file.path("..", "victim")
  expect_false(store_load_project(slug, dir = dir)$ok)
  expect_false(store_delete_project(slug, dir = dir)$ok)
  expect_true(file.exists(victim))
  expect_false(is_stored_slug("../x"))
  expect_false(is_stored_slug("a/b"))
  expect_true(is_stored_slug(project_slug("My project 1")))
})

test_that("a stale async result is told apart from the current one", {
  ep <- run_epoch()
  first <- ep$start()
  ep$bump()                      # the layer changed while it ran
  expect_false(ep$is_current(first))
  expect_true(ep$is_last_started(first))
  second <- ep$start()
  expect_true(ep$is_current(second))
  expect_false(ep$is_last_started(first))
})

test_that("removing a layer takes the results computed on it", {
  p <- tutorial_project()
  diff_p <- omicsCore::run_diff(p$experiments$proteomics, method = "limma",
                                group_col = "group", control_group = "Control",
                                case_group = "TreatA")
  integ <- structure(list(analysis_name = "run_integration",
                          params = list(experiments = c("proteomics", "rnaseq")),
                          input_info = list()),
                     class = class(diff_p))
  kept <- drop_layer_bundles(list(diff = diff_p, integration = integ), "rnaseq")
  expect_identical(names(kept), "diff")
  expect_length(drop_layer_bundles(list(diff = diff_p), "proteomics"), 0L)
})

test_that("a sample-ID column is never offered as the group", {
  meta <- data.frame(sample = paste0("S", 1:6), batch = c("a", "a", "b", "c", "d", "e"),
                     stringsAsFactors = FALSE)
  expect_length(grouping_candidates(meta), 0L)
  expect_identical(grouping_candidates(meta, min_per_level = 1L, replicated = TRUE),
                   "batch")
})

test_that("the chosen method and covariates survive a run", {
  proj <- shiny::reactiveVal(tutorial_project())
  shiny::testServer(diff_view_server, args = list(current_project = proj), {
    session$setInputs(layer = "proteomics", group_col = "group",
                      control = "Control", case = "TreatA",
                      method = "limma", covariates = "batch", rerun = 1)
    expect_false(is.null(diff_bundle()))
    html <- paste(as.character(output$ui_method), collapse = "")
    expect_match(html, '<option value="limma" selected>', fixed = TRUE)
  })
})

test_that("a project with no layers does not error in the Differential view", {
  proj <- shiny::reactiveVal(omicsCore::omics_project("empty", list()))
  shiny::testServer(diff_view_server, args = list(current_project = proj), {
    session$flushReact()
    expect_error(output$ui_method, class = "shiny.silent.error")
  })
})

test_that("two sessions keep two autosaves, and the oldest beyond the limit are pruned", {
  dir <- aa_store()
  p <- tutorial_project()
  p$name <- "tab one"
  expect_true(store_autosave(p, dir = dir, id = "aaaa1111"))
  Sys.sleep(1.1)
  p$name <- "tab two"
  expect_true(store_autosave(p, dir = dir, id = "bbbb2222"))
  snaps <- list_autosaves(dir)
  expect_identical(nrow(snaps), 2L)
  expect_identical(snaps$name, c("tab two", "tab one"))       # newest first
  expect_identical(store_read_autosave(dir)$name, "tab two")
  expect_identical(store_read_autosave(dir, path = basename(snaps$path[[2L]]))$name, "tab one")
  # A path from the browser that is not one of the snapshots reads nothing.
  expect_null(store_read_autosave(dir, path = "../elsewhere.omp"))
  # Autosaves are not listed as saved projects.
  expect_identical(nrow(list_saved_projects(dir)), 0L)

  for (i in 1:5) store_autosave(p, dir = dir, id = paste0("cccc", i))
  expect_lte(nrow(list_autosaves(dir)), MAX_AUTOSAVES + 1L)
})

test_that("a session writes its autosave under its own id", {
  dir <- aa_store()
  shiny::testServer(app_server, {
    id <- session_autosave_id(session)
    expect_match(id, "^[A-Za-z0-9]+$")
    session$setInputs(`project-load_tutorial` = 1)
    session$flushReact()
    expect_true(file.exists(autosave_path(dir, id)))
  })
})

test_that("an import parse that lands after a newer one is dropped", {
  ep <- run_epoch()
  old <- ep$start()
  new <- ep$start()
  expect_false(ep$is_current(old))
  expect_true(ep$is_current(new))
})
