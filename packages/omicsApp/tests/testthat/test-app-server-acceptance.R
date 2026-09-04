# The acceptance table from deploy/README.md, run against app_server().
#
# The table is what a person walks through on a deployed instance with
# one real dataset. Every row is a seam between two modules -- import
# and project state, diff and enrichment, state and autosave -- and the
# seams are where every reactive bug so far has lived: a replaced layer
# whose old results came back, a slider that recoloured a figure it was
# not supposed to touch, an autosave that read the state from before
# the result that triggered it. The module tests each check a part; this
# checks the wiring, on the same server function the app mounts.

skip_if_no_xlsx <- function() {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("readxl")
  skip_if_not_installed("withr")
}

# A private store per test, so autosave and Save as land somewhere this
# test owns and the next test starts from nothing.
local_store <- function(env = parent.frame()) {
  dir <- file.path(tempfile("acceptance-"), "store")
  dir.create(dir, recursive = TRUE)
  withr::local_envvar(OMICSAPP_DATA_DIR = dir, .local_envir = env)
  withr::defer(unlink(dirname(dir), recursive = TRUE), envir = env)
  dir
}

workbook <- function(seed, env = parent.frame()) {
  path <- tempfile(fileext = ".xlsx")
  withr::defer(unlink(path), envir = env)
  write_tiny_omics_xlsx(path, n_features = 40L, n_samples = 8L, seed = seed)
  path
}

upload <- function(session, path, name = basename(path)) {
  # The import view parses with the modality's default label first and
  # relabels from the values once it has them; the scale check warns on
  # that first pass. That is the view's own noise, not this test's.
  suppressWarnings(session$setInputs(
    `import-omics_type` = "proteomics",
    `import-file` = list(datapath = path, name = name,
                         size = file.info(path)$size)
  ))
}

run_diff_in_app <- function(session, click) {
  session$setInputs(`diff-method` = "ttest", `diff-group_col` = "group",
                    `diff-control` = "G1", `diff-case` = "G2",
                    `diff-rerun` = click)
}

# plotly and DT stamp their JSON with ids drawn from tempfile(); two
# renders of the same figure differ only there.
strip_ids <- function(json) gsub("file[0-9a-f]{6,}", "ID", as.character(json))

autosaved <- function(dir) omicsCore::load_project(file.path(dir, "_autosave.omp"))

# ---- rows 1-3: import, re-import, replace --------------------------------

test_that("upload and confirm: the project appears, and is on disk before anything else happens", {
  skip_if_no_xlsx()
  store <- local_store()
  file_a <- workbook(1)
  shiny::testServer(app_server, {
    expect_null(current_project())
    expect_false(file.exists(file.path(store, "_autosave.omp")))

    upload(session, file_a, "cohort_a.xlsx")
    session$setInputs(`import-confirm` = 1)

    proj <- current_project()
    expect_true(omicsCore::is_omics_project(proj))
    expect_identical(names(proj$experiments), "proteomics")
    expect_identical(dim(proj$experiments$proteomics$expr_mat), c(40L, 8L))

    # Row 9 of the table: the store carries the snapshot and the upload
    expect_true(file.exists(file.path(store, "_autosave.omp")))
    expect_identical(names(autosaved(store)$experiments), "proteomics")
    raw <- list.files(file.path(store, "raw"))
    expect_length(raw, 1L)
    expect_match(raw, "^cohort_a__[0-9a-f]{12}\\.xlsx$")
  })
})

test_that("re-uploading the same file clears nothing", {
  skip_if_no_xlsx()
  store <- local_store()
  file_a <- workbook(1)
  shiny::testServer(app_server, {
    upload(session, file_a)
    session$setInputs(`import-confirm` = 1)
    run_diff_in_app(session, 1)
    expect_true("diff" %in% names(current_project()$bundles))
    before <- layer_generation()
    before_bundle <- diff_view$bundle()

    upload(session, file_a)
    session$setInputs(`import-confirm` = 2)

    expect_identical(layer_generation(), before)
    expect_identical(diff_view$bundle(), before_bundle)
    expect_true("diff" %in% names(current_project()$bundles))
    expect_true("diff" %in% names(autosaved(store)$bundles))
  })
})

test_that("uploading a different file asks, and confirming clears the analyses", {
  skip_if_no_xlsx()
  store <- local_store()
  file_a <- workbook(1)
  file_b <- workbook(2)
  shiny::testServer(app_server, {
    upload(session, file_a)
    session$setInputs(`import-confirm` = 1)
    run_diff_in_app(session, 1)
    expect_true("diff" %in% names(current_project()$bundles))
    matrix_a <- current_project()$experiments$proteomics$expr_mat
    before <- layer_generation()

    upload(session, file_b)
    # Confirm alone does not replace: the data changed, so a dialog
    # names what will be cleared and waits.
    session$setInputs(`import-confirm` = 2)
    expect_identical(current_project()$experiments$proteomics$expr_mat, matrix_a)
    expect_true("diff" %in% names(current_project()$bundles))

    session$setInputs(`import-confirm_replace` = 1)
    proj <- current_project()
    expect_false(identical(proj$experiments$proteomics$expr_mat, matrix_a))
    expect_gt(layer_generation(), before)
    # The results computed on the old data are gone from the views, the
    # project and the snapshot -- not merely from one of the three.
    expect_null(diff_view$bundle())
    expect_false("diff" %in% names(proj$bundles))
    expect_false("enrich" %in% names(proj$bundles))
    expect_false("diff" %in% names(autosaved(store)$bundles))
  })
})

# ---- rows 4-5: the differential view and its sliders ----------------------

test_that("the volcano is two-coloured at a stated cut, and the sliders do not touch it", {
  skip_if_no_xlsx()
  store <- local_store()
  file_a <- workbook(1)
  shiny::testServer(app_server, {
    upload(session, file_a)
    session$setInputs(`import-confirm` = 1)
    run_diff_in_app(session, 1)

    volcano_before <- strip_ids(output$`diff-volcano`)
    # The hit table is served server-side, so its JSON is only a handle;
    # the stat cards next to it carry the counts the sliders change.
    hits_before <- strip_ids(output$`diff-stats`)
    # Drawn at the default cut: the dashed line sits at -log10(0.05).
    # (ggplotly does not carry the caption over, so the line is the
    # evidence of the threshold the figure was drawn at.)
    expect_match(volcano_before, "yintercept: 1.30103", fixed = TRUE)
    expect_match(volcano_before, "-log10(adj_p_value)", fixed = TRUE)
    expect_identical(diff_view$thresholds()$p_cutoff, 0.05)

    # Open the thresholds right up
    session$setInputs(`diff-fdr_cut` = 0.9, `diff-fc_cut` = 0)
    session$elapse(400)   # the sliders are debounced
    expect_identical(diff_view$thresholds()$p_cutoff, 0.9)
    expect_identical(diff_view$thresholds()$effect_cutoff, 0)

    volcano_after <- strip_ids(output$`diff-volcano`)
    hits_after <- strip_ids(output$`diff-stats`)
    expect_false(identical(hits_before, hits_after))
    expect_identical(volcano_before, volcano_after)
  })
})

test_that("enrichment reads the thresholds the hit table is read at", {
  skip_if_no_xlsx()
  skip_if_not_installed("clusterProfiler")
  store <- local_store()
  file_a <- workbook(1)
  shiny::testServer(app_server, {
    upload(session, file_a)
    session$setInputs(`import-confirm` = 1)
    run_diff_in_app(session, 1)
    session$setInputs(`diff-fdr_cut` = 0.5, `diff-fc_cut` = 0.1,
                      `diff-p_kind` = "raw")
    session$elapse(400)
    th <- diff_view$thresholds()
    expect_identical(th$p_preference, "raw")
    session$setInputs(`enrich-type` = "ora", `enrich-database` = "hallmark",
                      `enrich-direction` = "both", `enrich-rerun` = 1)
    eb <- enrich_bundle()
    skip_if(is.null(eb), "enrichment did not produce a bundle on the tiny fixture")
    expect_identical(eb$params$p_cutoff, th$p_cutoff)
    expect_identical(eb$params$p_preference, th$p_preference)
    expect_identical(eb$params$effect_cutoff, th$effect_cutoff)
  })
})

# ---- row 6: the exported script ------------------------------------------

test_that("the report view's analysis code names the calls that were run", {
  skip_if_no_xlsx()
  store <- local_store()
  file_a <- workbook(1)
  shiny::testServer(app_server, {
    upload(session, file_a)
    session$setInputs(`import-confirm` = 1)
    run_diff_in_app(session, 1)
    script <- omicsCore::export_script(current_project())
    expect_no_error(parse(text = script))
    expect_true(any(grepl("run_diff(", script, fixed = TRUE)))
    expect_true(any(grepl("read_omics(", script, fixed = TRUE)))
    # It points at the archived upload, not at Shiny's temp copy
    archived <- current_project()$experiments$proteomics$source_path
    expect_true(file.exists(archived))
    expect_true(any(grepl(basename(archived), script, fixed = TRUE)))
  })
})

# ---- rows 7-8: save, and come back tomorrow ---------------------------------

test_that("Save as puts the project under My projects, and Open brings it back", {
  skip_if_no_xlsx()
  store <- local_store()
  file_a <- workbook(1)
  shiny::testServer(app_server, {
    upload(session, file_a)
    session$setInputs(`import-confirm` = 1)
    run_diff_in_app(session, 1)

    session$setInputs(`project-save_name` = "Cohort A, first pass",
                      `project-save_project` = 1)
    saved <- list_saved_projects(store)
    expect_identical(saved$slug, "Cohort_A,_first_pass")
    expect_true(file.exists(saved$path))

    # Lose the in-memory state, as a page reload would
    current_project(NULL)
    session$flushReact()
    expect_null(current_project())

    session$setInputs(`project-saved_pick` = saved$slug,
                      `project-open_project` = 1)
    proj <- current_project()
    expect_true(omicsCore::is_omics_project(proj))
    expect_identical(names(proj$experiments), "proteomics")
    expect_true("diff" %in% names(proj$bundles))
  })
})

test_that("a new session on the same store restores the last one without being asked", {
  skip_if_no_xlsx()
  store <- local_store()
  file_a <- workbook(1)
  shiny::testServer(app_server, {
    upload(session, file_a)
    session$setInputs(`import-confirm` = 1)
    run_diff_in_app(session, 1)
    expect_true("diff" %in% names(autosaved(store)$bundles))
  })

  # The container was recycled; the user logs back in.
  shiny::testServer(app_server, {
    session$flushReact()
    proj <- current_project()
    expect_true(omicsCore::is_omics_project(proj))
    expect_identical(names(proj$experiments), "proteomics")
    expect_true("diff" %in% names(proj$bundles))
    # And the views show the restored layer rather than the demo
    expect_identical(diff_view$layer(), "proteomics")
  })
})

# ---- the Project view's own links ------------------------------------------

test_that("the Experiments table's View link lands on QC with that layer selected", {
  skip_if_no_xlsx()
  store <- local_store()
  file_a <- workbook(1)
  shiny::testServer(app_server, {
    upload(session, file_a)
    session$setInputs(`import-confirm` = 1)
    expect_identical(current_view(), "project")
    session$setInputs(`project-view_layer_1` = 1)
    expect_identical(current_view(), "qc")
    expect_identical(requested_layer(), "proteomics")
  })
})

test_that("the sidebar switches views and never lands on a view that does not exist", {
  shiny::testServer(app_server, {
    session$setInputs(nav_enrich = 1)
    expect_identical(current_view(), "enrich")
    session$setInputs(nav_report = 1)
    expect_identical(current_view(), "report")
    set_view("no_such_view")
    expect_identical(current_view(), "report")
  })
})
