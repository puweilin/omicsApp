# The tutorial: an example project every step can be run on, and the
# Project view's guide through it.

tut_store <- function(env = parent.frame()) {
  dir <- withr::local_tempdir(.local_envir = env)
  withr::local_envvar(OMICSAPP_DATA_DIR = dir, .local_envir = env)
  dir
}

tut_html <- function(x) paste(as.character(x), collapse = "")

test_that("the tutorial project is a real two-layer project the integration can join", {
  p <- tutorial_project()
  expect_true(omicsCore::is_omics_project(p))
  expect_identical(names(p$experiments), c("proteomics", "rnaseq"))
  expect_identical(p$experiments$rnaseq$assay_type, "raw_count")
  prot_sym <- p$experiments$proteomics$feature_df$feature_symbol
  rna_sym <- p$experiments$rnaseq$feature_df$feature_symbol
  expect_equal(length(intersect(prot_sym, rna_sym)), length(prot_sym))
  # One control and two treatments.
  expect_setequal(unique(p$experiments$proteomics$meta_df$group),
                  c("Control", "TreatA", "TreatB"))
  # Distinct sample ids, paired by the donor column.
  expect_length(intersect(colnames(p$experiments$proteomics$expr_mat),
                          colnames(p$experiments$rnaseq$expr_mat)), 0L)
  pr <- omicsCore::sample_pairing_preview(p, "proteomics", "rnaseq")
  expect_identical(pr$source, "donor")
  expect_equal(nrow(pr$pairs), 12L)
})

test_that("each step of the tutorial finds what was put there", {
  p <- tutorial_project()
  d <- omicsCore::run_diff(p$experiments$proteomics, method = "limma",
                           group_col = "group", control_group = "Control",
                           case_group = c("TreatA", "TreatB"))
  s <- omicsCore::summarize_diff_contrasts(d, effect_cutoff = log2(1.2))
  expect_gte(s$n_up[s$comparison == "TreatA_vs_Control"], 25L)
  expect_gte(s$n_down[s$comparison == "TreatB_vs_Control"], 25L)

  a <- omicsCore::select_comparison(d, "TreatA_vs_Control")
  r <- omicsCore::run_diff(p$experiments$rnaseq, method = "edger",
                           group_col = "group", control_group = "Control",
                           case_group = "TreatA")
  i <- omicsCore::run_integration(p, "concordance", c("proteomics", "rnaseq"),
                                  diff_bundles = list(proteomics = a, rnaseq = r))
  expect_gte(i$params$method_info$n_concordant, 20L)
  expect_s3_class(omicsCore::run_integration(p, "correlation",
                                             c("proteomics", "rnaseq")),
                  "analysis_bundle")
})

test_that("building the tutorial leaves the session's random stream alone", {
  .example_cache$tutorial <- NULL
  set.seed(99)
  before <- .Random.seed
  tutorial_project()
  expect_identical(.Random.seed, before)
})

test_that("the first page welcomes, explains the workflow and offers the example", {
  tut_store()
  nav <- character(0)
  shiny::testServer(project_view_server, args = list(
    current_project = shiny::reactiveVal(NULL),
    navigate = function(v) nav <<- c(nav, v)), {
    # The autosave restore runs once on arrival; nothing to restore here.
    session$flushReact()
    html <- tut_html(output$guide)
    expect_match(html, "Welcome to omicsApp", fixed = TRUE)
    expect_match(html, "Try the example project", fixed = TRUE)
    for (step in c("Import", "Quality control", "Differential", "Enrichment",
                   "Integration", "Report")) {
      expect_match(html, step, fixed = TRUE)
    }
    session$setInputs(go_import = 1)
    expect_identical(nav, "import")

    session$setInputs(load_tutorial = 1)
    expect_true(is_tutorial_project(current_project()))
    html <- tut_html(output$guide)
    expect_match(html, "Tutorial", fixed = TRUE)
    expect_match(html, "next: Quality control", fixed = TRUE)
    session$setInputs(go_next = 1)
    expect_identical(nav, c("import", "qc"))
  })
})

test_that("the checklist follows what has been run", {
  p <- tutorial_project()
  expect_identical(WORKFLOW_STEPS$id[workflow_progress(p)$next_step], "qc")
  p$bundles <- list(qc = 1, diff = 1)
  expect_identical(WORKFLOW_STEPS$id[workflow_progress(p)$next_step], "enrich")
  p$bundles <- list(qc = 1, diff = 1, enrich = 1)
  expect_identical(WORKFLOW_STEPS$id[workflow_progress(p)$next_step], "integration")
  one <- p
  one$experiments$rnaseq <- NULL
  prog <- workflow_progress(one)
  expect_identical(WORKFLOW_STEPS$id[prog$next_step], "report")
  expect_true(prog$skip[WORKFLOW_STEPS$id == "integration"])
})

test_that("the example does not replace a loaded project without asking", {
  tut_store()
  mine <- omicsCore::omics_project("Mine", list(proteomics = example_input("proteomics")))
  shiny::testServer(project_view_server, args = list(
    current_project = shiny::reactiveVal(mine)), {
    session$flushReact()
    session$setInputs(load_tutorial = 1)
    expect_identical(current_project()$name, "Mine")
    session$setInputs(confirm_tutorial = 1)
    expect_true(is_tutorial_project(current_project()))
  })
})

test_that("deleting a saved project asks first", {
  tut_store()
  store_save_project(omicsCore::omics_project("Keep", list(
    proteomics = example_input("proteomics"))), "keep")
  shiny::testServer(project_view_server, args = list(
    current_project = shiny::reactiveVal(NULL)), {
    session$flushReact()
    session$setInputs(saved_pick = "keep", delete_project = 1)
    expect_true("keep" %in% list_saved_projects()$slug)
    session$setInputs(confirm_delete_project = 1)
    expect_false("keep" %in% list_saved_projects()$slug)
  })
})

test_that("the Project view counts every layer's samples and the analyses run", {
  tut_store()
  p <- tutorial_project()
  p$bundles <- list(qc = 1, diff = 1)
  shiny::testServer(project_view_server, args = list(
    current_project = shiny::reactiveVal(p)), {
    session$flushReact()
    html <- tut_html(output$stats)
    expect_match(html, ">24<", fixed = TRUE)
    expect_match(html, "qc · diff", fixed = TRUE)
  })
})
