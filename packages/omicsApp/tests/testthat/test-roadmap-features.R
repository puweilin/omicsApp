# The second round: contrast modes, the global test, enrichment across
# comparisons, the design recorded at import, and a project that keeps
# every comparison of a run.

rf_store <- function(env = parent.frame()) {
  dir <- withr::local_tempdir(.local_envir = env)
  withr::local_envvar(OMICSAPP_DATA_DIR = dir, .local_envir = env)
  dir
}

rf_html <- function(x) paste(as.character(x), collapse = "")

test_that("all pairs and custom comparisons run from the Differential view", {
  proj <- shiny::reactiveVal(tutorial_project())
  shiny::testServer(diff_view_server, args = list(current_project = proj), {
    session$setInputs(layer = "proteomics", group_col = "group",
                      control = "Control", method = "limma",
                      contrast_mode = "pairwise", rerun = 1)
    expect_identical(omicsCore::diff_comparisons(diff_bundle()),
                     c("TreatA_vs_Control", "TreatB_vs_Control", "TreatB_vs_TreatA"))

    session$setInputs(contrast_mode = "custom",
                      custom_contrasts = "(TreatA + TreatB)/2 - Control\nTreatB - TreatA",
                      rerun = 2)
    expect_identical(omicsCore::diff_comparisons(diff_bundle()),
                     c("(TreatA + TreatB)/2 - Control", "TreatB_vs_TreatA"))
    session$setInputs(comparison = "(TreatA + TreatB)/2 - Control")
    shown <- session$returned$bundle()
    expect_null(shown$params$case_group)
    expect_identical(shown$params$contrasts, "(TreatA + TreatB)/2 - Control")

    # A contrast that is not one says so, and nothing runs.
    session$setInputs(custom_contrasts = "TreatA + TreatB", rerun = 3)
    expect_match(diff_error(), "sum to")
  })
})

test_that("the project keeps every comparison, and which one was on screen", {
  proj <- shiny::reactiveVal(tutorial_project())
  shiny::testServer(diff_view_server, args = list(current_project = proj), {
    session$setInputs(layer = "proteomics", group_col = "group",
                      control = "Control", case = c("TreatA", "TreatB"),
                      method = "limma", rerun = 1)
    session$setInputs(comparison = "TreatB_vs_Control")
    pb <- session$returned$project_bundle()
    expect_identical(omicsCore::diff_comparisons(pb),
                     c("TreatA_vs_Control", "TreatB_vs_Control"))
    expect_identical(pb$params$shown_comparison, "TreatB_vs_Control")
    # And the overlap plot draws.
    expect_false(is.null(output$contrast_summary))
  })
})

test_that("the global test runs from the Differential view for counts and intensities", {
  proj <- shiny::reactiveVal(tutorial_project())
  shiny::testServer(diff_view_server, args = list(current_project = proj), {
    session$setInputs(layer = "proteomics", group_col = "group", method = "limma")
    session$setInputs(run_anova = 1)
    b <- anova_bundle()
    expect_identical(b$params$analysis_type, "anova")
    expect_match(rf_html(output$anova_summary), "differ between the groups")

    session$setInputs(layer = "rnaseq", method = "edger")
    expect_null(anova_bundle())
    session$setInputs(run_anova = 2)
    expect_identical(anova_bundle()$params$method, "edger")
    # The inflammatory, OXPHOS and E2F genes (102 of 252) carry signal.
    df <- anova_hits()
    expect_gte(sum(df$adj_p_value < 0.05, na.rm = TRUE), 60L)
  })
})

test_that("the Enrichment view enriches every comparison side by side", {
  skip_if_not(has_pkg("clusterProfiler"))
  p <- tutorial_project()
  full <- omicsCore::run_diff(p$experiments$proteomics, method = "limma",
                              group_col = "group", control_group = "Control",
                              case_group = c("TreatA", "TreatB"))
  shown <- omicsCore::select_comparison(full, "TreatA_vs_Control")
  shiny::testServer(enrich_view_server, args = list(
    diff_bundle = shiny::reactiveVal(shown),
    diff_all = shiny::reactiveVal(full),
    diff_thresholds = shiny::reactive(list(p_cutoff = 0.05, p_preference = "adjusted",
                                           effect_cutoff = log2(1.2)))), {
    session$setInputs(type = "ora", database = "hallmark", direction = "both")
    expect_false(is.null(output$compare_card))
    session$setInputs(run_compare = 1)
    ce <- compare_bundle()
    expect_identical(ce$analysis_name, "compare_enrichment")
    df <- ce$results$enrich_result_df
    top <- function(cmp) {
      sub <- df[df$comparison == cmp, ]
      sub$pathway_id[which.min(sub$adj_p_value)]
    }
    expect_identical(top("TreatA_vs_Control"), "HALLMARK_INFLAMMATORY_RESPONSE")
    expect_identical(top("TreatB_vs_Control"), "HALLMARK_OXIDATIVE_PHOSPHORYLATION")
  })
})

test_that("the design chosen at import is recorded and becomes the default contrast", {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("readxl")
  rf_store()
  xlsx <- tempfile(fileext = ".xlsx")
  on.exit(unlink(xlsx), add = TRUE)
  write_tiny_omics_xlsx(xlsx, n_features = 20L, n_samples = 8L)
  shiny::testServer(import_view_server, {
    session$setInputs(omics_type = "proteomics",
                      file = list(datapath = xlsx, name = "tiny.xlsx",
                                  size = file.info(xlsx)$size))
    expect_match(rf_html(output$confirm_design), "Study design", fixed = TRUE)
    session$setInputs(design_group = "group", design_reference = "G2", confirm = 1)
    inp <- session$returned()
    expect_identical(omicsCore::study_design(inp),
                     list(group_col = "group", reference = "G2"))
  })

  # The Differential view then compares against G2, not the first level.
  inp <- omicsCore::set_study_design(
    omicsCore::read_omics(xlsx, omics_type = "proteomics",
                          assay_type = "normalized_intensity")$input, "group", "G2")
  proj <- shiny::reactiveVal(omicsCore::omics_project("d", list(proteomics = inp)))
  shiny::testServer(diff_view_server, args = list(current_project = proj), {
    session$flushReact()
    expect_identical(default_contrast()$control, "G2")
    expect_identical(default_contrast()$case, "G1")
  })
})

test_that("Integration repeats a custom contrast on the partner layer", {
  p <- tutorial_project()
  full <- omicsCore::run_diff(p$experiments$proteomics, method = "limma",
                              group_col = "group",
                              contrasts = c("(TreatA + TreatB)/2 - Control",
                                            "TreatB - TreatA"))
  shown <- omicsCore::select_comparison(full, "TreatB_vs_TreatA")
  shiny::testServer(integration_view_server, args = list(
    current_project = shiny::reactiveVal(p),
    diff_bundle = shiny::reactiveVal(shown),
    diff_layer = shiny::reactiveVal("proteomics")), {
    session$flushReact()
    expect_true(isTRUE(can_run()$ok))
    sec <- sec_cache()$bundle
    expect_identical(sec$params$comparison, "TreatB_vs_TreatA")
    expect_identical(sec$params$all_contrasts,
                     c("(TreatA + TreatB)/2 - Control", "TreatB - TreatA"))
    expect_s3_class(integration_bundle(), "analysis_bundle")
  })
})

test_that("the report view lists every comparison of a run", {
  p <- tutorial_project()
  p$bundles <- list(diff = omicsCore::run_diff(
    p$experiments$proteomics, method = "limma", group_col = "group",
    control_group = "Control", case_group = c("TreatA", "TreatB")))
  shiny::testServer(report_view_server, args = list(
    current_project = shiny::reactiveVal(p)), {
    html <- rf_html(output$bundle_cards)
    expect_match(html, "2 comparisons: TreatA vs Control, TreatB vs Control", fixed = TRUE)
    expect_match(html, "Results in the report", fixed = TRUE)
  })
})

test_that("the import button is disabled until a file parses", {
  # Package sources are not there under R CMD check.
  src_file <- test_path("..", "..", "R", "mod_import_view.R")
  skip_if_not(file.exists(src_file), "package sources not beside the tests")
  src <- paste(readLines(src_file), collapse = "\n")
  expect_match(src, 'shinyjs::toggleState("confirm", condition = ok)', fixed = TRUE)
})

test_that("the page carries busy indicators", {
  html <- rf_html(app_ui())
  expect_match(html, "busy", ignore.case = TRUE)
})
