# Inputs the UI would never send, sent anyway.
#
# A browser sends what its controls allow; a script talking to the
# websocket, a stale tab after a deploy, or a numericInput mid-edit
# sends anything: a negative cutoff, NA, a level that is not in the
# column, a method that does not exist, a project name that is a
# script tag. In Shiny an error inside an observer ends the session
# for that user, so each of these has to be answered in the notices,
# not thrown. Every scenario here must leave the session alive and
# the last good result in place.

skip_if_no_xlsx <- function() {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("readxl")
  skip_if_not_installed("withr")
}

local_store <- function(env = parent.frame()) {
  dir <- file.path(tempfile("hostile-"), "store")
  dir.create(dir, recursive = TRUE)
  withr::local_envvar(OMICSAPP_DATA_DIR = dir, .local_envir = env)
  withr::defer(unlink(dirname(dir), recursive = TRUE), envir = env)
  dir
}

fixture_xlsx <- function(env = parent.frame()) {
  path <- tempfile(fileext = ".xlsx")
  withr::defer(unlink(path), envir = env)
  write_tiny_omics_xlsx(path, n_features = 40L, n_samples = 8L, seed = 1)
  path
}

upload_and_run <- function(session, path) {
  suppressWarnings(session$setInputs(
    `import-omics_type` = "proteomics",
    `import-file` = list(datapath = path, name = "a.xlsx", size = file.info(path)$size)))
  session$setInputs(`import-confirm` = 1)
  session$setInputs(`diff-method` = "ttest", `diff-group_col` = "group",
                    `diff-control` = "G1", `diff-case` = "G2", `diff-rerun` = 1)
}

# req() cancels a render with a condition whose message is empty;
# anything else is an error the user would see.
render_state <- function(output, id) {
  tryCatch({ force(output[[id]]); "rendered" },
           error = function(e) if (nzchar(conditionMessage(e))) conditionMessage(e) else "waiting")
}

watched <- c("diff-stats", "diff-notices", "diff-volcano", "diff-hits", "diff-header",
             "enrich-notices", "enrich-header", "enrich-dot", "qc-stats", "qc-pca",
             "qc-missing", "qc-header", "project-body", "project-storage", "project-stats",
             "import-schema_summary", "import-upload_status", "import-confirm_state",
             "report-bundle_cards", "report-notices", "integration-stats",
             "integration-notices", "project_picker")

scenarios <- list(
  `negative FDR cutoff` = list(`diff-fdr_cut` = -1),
  `NA FDR cutoff` = list(`diff-fdr_cut` = NA),
  `text FDR cutoff` = list(`diff-fdr_cut` = "abc"),
  `FDR cutoff above one` = list(`diff-fdr_cut` = 2),
  `NULL FDR cutoff` = list(`diff-fdr_cut` = NULL),
  `NA fold-change cutoff` = list(`diff-fc_cut` = NA),
  `negative fold-change cutoff` = list(`diff-fc_cut` = -5),
  `infinite fold-change cutoff` = list(`diff-fc_cut` = Inf),
  `negative label count` = list(`diff-label_top` = -1),
  `huge label count` = list(`diff-label_top` = 1e6),
  `NA label count` = list(`diff-label_top` = NA),
  `text label count` = list(`diff-label_top` = "ten"),
  `group column that does not exist` = list(`diff-group_col` = "nope", `diff-rerun` = 2),
  `control equal to case` = list(`diff-control` = "G1", `diff-case` = "G1", `diff-rerun` = 2),
  `level that is not in the column` = list(`diff-case` = "G9", `diff-rerun` = 2),
  `method that does not exist` = list(`diff-method` = "bogus", `diff-rerun` = 2),
  `layer that does not exist` = list(`diff-layer` = "nope", `diff-rerun` = 2),
  `p kind that does not exist` = list(`diff-p_kind` = "bogus"),
  `covariate that does not exist` = list(`diff-covariates` = "nope", `diff-rerun` = 2),
  `covariate that is the group` = list(`diff-covariates` = "group", `diff-rerun` = 2),
  `enrichment database that does not exist` = list(`enrich-database` = "bogus", `enrich-type` = "ora", `enrich-rerun` = 1),
  `enrichment type that does not exist` = list(`enrich-type` = "bogus", `enrich-rerun` = 1),
  `enrichment direction that does not exist` = list(`enrich-direction` = "sideways", `enrich-rerun` = 1),
  `enrichment display thresholds` = list(`enrich-show_p` = -1, `enrich-show_cutoff` = NA),
  `negative QC threshold` = list(`qc-missing_threshold` = -1),
  `QC threshold above one` = list(`qc-missing_threshold` = 2),
  `NA QC threshold` = list(`qc-missing_threshold` = NA),
  `outlier method that does not exist` = list(`qc-outlier_method` = "bogus"),
  `imputation method that does not exist` = list(`qc-impute_method` = "bogus"),
  `QC layer that does not exist` = list(`qc-layer` = "nope"),
  `QC view that does not exist` = list(`qc-quality_view` = "bogus"),
  `project name that is a script tag` = list(`project-save_name` = "<script>alert(1)</script>", `project-save_project` = 1),
  `project name that escapes the store` = list(`project-save_name` = "../../escape", `project-save_project` = 1),
  `empty project name` = list(`project-save_name` = "", `project-save_project` = 1),
  `opening a path outside the store` = list(`project-saved_pick` = "../../etc/passwd", `project-open_project` = 1),
  `deleting a project that does not exist` = list(`project-saved_pick` = "nope", `project-delete_project` = 1),
  `dropping a layer nobody chose` = list(`project-confirm_drop_layer` = 1),
  `omics type that does not exist` = list(`import-omics_type` = "bogus"),
  `assay type that does not exist` = list(`import-assay_type` = "bogus"),
  `upload whose file is gone` = list(`import-file` = list(datapath = "/nonexistent/x.xlsx", name = "x.xlsx", size = 0)),
  `integration with one layer` = list(`integration-rerun` = 1),
  `accepting a pairing that was never offered` = list(`integration-accept_pairing` = 1),
  `every nav link at once` = list(nav_import = 1, nav_qc = 1, nav_diff = 1, nav_enrich = 1,
                                  nav_integration = 1, nav_report = 1, nav_project = 1)
)

test_that("no hostile input ends the session or loses the last result", {
  skip_if_no_xlsx()
  store <- local_store()
  xlsx <- fixture_xlsx()
  for (nm in names(scenarios)) {
    hostile <- scenarios[[nm]]
    shiny::testServer(app_server, {
      upload_and_run(session, xlsx)
      expect_true("diff" %in% names(current_project()$bundles), label = nm)
      before <- diff_view$bundle()
      # The `bogus` omics type is answered by a warning from omicsCore
      # and then ignored; that warning is the point of the scenario.
      suppressWarnings(do.call(session$setInputs, hostile))
      expect_true(omicsCore::is_omics_project(current_project()), label = nm)
      expect_identical(diff_view$bundle(), before, label = paste(nm, "keeps the result"))
      states <- vapply(watched, function(id) render_state(output, id), character(1))
      bad <- states[!states %in% c("rendered", "waiting")]
      expect_identical(unname(bad), character(0),
                       label = sprintf("%s: outputs %s", nm, paste(names(bad), collapse = ", ")))
    })
  }
})

test_that("a hostile upload is refused in the import view, not thrown", {
  skip_if_no_xlsx()
  store <- local_store()
  garbage <- withr::local_tempfile(fileext = ".xlsx")
  writeLines("this is not a workbook", garbage)
  empty <- withr::local_tempfile(fileext = ".csv")
  file.create(empty)
  for (path in c(garbage, empty)) {
    shiny::testServer(app_server, {
      suppressWarnings(session$setInputs(
        `import-omics_type` = "proteomics",
        `import-file` = list(datapath = path, name = basename(path), size = file.size(path))))
      session$setInputs(`import-confirm` = 1)
      expect_null(current_project())
      state <- render_state(output, "import-upload_status")
      expect_true(state %in% c("rendered", "waiting"), label = state)
    })
  }
})

test_that("confirming with nothing uploaded, and restoring with nothing saved, are no-ops", {
  skip_if_no_xlsx()
  store <- local_store()
  shiny::testServer(app_server, {
    session$setInputs(`import-confirm` = 1)
    session$setInputs(`import-confirm_replace` = 1)
    session$setInputs(`project-restore_autosave` = 1)
    session$setInputs(`enrich-database` = "hallmark", `enrich-type` = "ora", `enrich-rerun` = 1)
    expect_null(current_project())
    expect_false(file.exists(file.path(store, "_autosave.omp")))
  })
})

test_that("clicking Re-run n times runs the engine n times", {
  skip_if_no_xlsx()
  store <- local_store()
  xlsx <- fixture_xlsx()
  calls <- 0L
  real <- omicsCore::run_diff
  testthat::local_mocked_bindings(
    run_diff = function(...) { calls <<- calls + 1L; real(...) },
    .package = "omicsCore")
  shiny::testServer(app_server, {
    upload_and_run(session, xlsx)
    after_first <- calls
    for (i in 2:6) session$setInputs(`diff-rerun` = i)
    # An observer created inside an observer would run once per earlier
    # click as well, and the count would be 1 + 2 + ... rather than 5.
    expect_identical(calls - after_first, 5L)
  })
})

test_that("a hostile project name is stored under a safe slug and shown as text", {
  skip_if_no_xlsx()
  store <- local_store()
  xlsx <- fixture_xlsx()
  shiny::testServer(app_server, {
    upload_and_run(session, xlsx)
    session$setInputs(`project-save_name` = "<script>alert(1)</script> cohort",
                      `project-save_project` = 1)
    saved <- list.files(store, pattern = "\\.omp$")
    saved <- saved[saved != "_autosave.omp"]
    expect_length(saved, 1L)
    expect_false(grepl("[<>]", saved))
    body <- as.character(output[["project-body"]]$html)
    expect_false(grepl("<script>alert(1)</script>", body, fixed = TRUE))
  })
})
