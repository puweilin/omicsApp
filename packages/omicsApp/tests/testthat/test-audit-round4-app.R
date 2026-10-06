# App defects found in the round-4 UX and UI audits.

r4_store <- function(env = parent.frame()) {
  dir <- withr::local_tempdir(.local_envir = env)
  withr::local_envvar(OMICSAPP_DATA_DIR = dir, .local_envir = env)
  dir
}

r4_layer <- function(meta, n_feat = 60L, seed = 1L, shift_col = NULL, shift_level = NULL) {
  set.seed(seed)
  ids <- rownames(meta)
  m <- matrix(stats::rnorm(n_feat * length(ids), 20, 0.3), n_feat,
              dimnames = list(paste0("GENE", seq_len(n_feat)), ids))
  if (!is.null(shift_col)) m[1:10, meta[[shift_col]] == shift_level] <-
    m[1:10, meta[[shift_col]] == shift_level] + 2
  omicsCore::omics_input(m, meta,
                         data.frame(feature_id = rownames(m), feature_symbol = rownames(m)),
                         omics_type = "proteomics", assay_type = "normalized_intensity")
}

test_that("the default control is the control in common namings", {
  expect_identical(default_control_level(c("After", "Before")), "Before")
  expect_identical(default_control_level(c("Post", "Pre")), "Pre")
  expect_identical(default_control_level(c("LPS", "PBS")), "PBS")
  expect_identical(default_control_level(c("处理组", "对照组")),
                   "对照组")
})

test_that("group columns: coded numbers are offered, nuisance columns come last", {
  meta <- data.frame(arm = rep(0:1, each = 4), batch = rep(c("b1", "b2"), 4),
                     Diet = rep(c("chow", "hfd", "chow", "hfd"), each = 2),
                     age = 31:38)
  cands <- grouping_candidates(meta)
  expect_true("arm" %in% cands)
  expect_false("age" %in% cands)
  expect_identical(utils::tail(cands, 1L), "batch")
  expect_identical(cands[[1L]], "arm")   # a stated group name wins
})

test_that("a before/after design is paired by default and finds the change", {
  ids <- paste0("s", 1:12)
  meta <- data.frame(patient = rep(paste0("P", 1:6), 2),
                     timepoint = rep(c("Before", "After"), each = 6), row.names = ids)
  inp <- r4_layer(meta, shift_col = "timepoint", shift_level = "After")
  # A large between-patient effect, which the pairing removes.
  inp$expr_mat <- inp$expr_mat + matrix(rep(stats::rnorm(6, 0, 3), 2), nrow(inp$expr_mat), 12,
                                        byrow = TRUE)
  proj <- shiny::reactiveVal(omicsCore::omics_project("paired", list(proteomics = inp)))
  shiny::testServer(diff_view_server, args = list(current_project = proj), {
    session$flushReact()
    expect_identical(pairing_candidates(), "patient")
    expect_identical(default_contrast()$control, "Before")
    session$setInputs(layer = "proteomics", group_col = "timepoint", control = "Before",
                      case = "After", method = "limma", paired_col = "patient", rerun = 1)
    b <- diff_bundle()
    expect_identical(b$params$paired_col, "patient")
    expect_gte(sum(b$results$diff_result_df$adj_p_value < 0.05, na.rm = TRUE), 8L)
  })
})

test_that("a continuous variable can be tested for a trend", {
  ids <- paste0("s", 1:15)
  meta <- data.frame(dose = rep(c(0, 1, 5, 10, 20), 3), row.names = ids)
  inp <- r4_layer(meta)
  inp$expr_mat[1:10, ] <- inp$expr_mat[1:10, ] + rep(meta$dose / 10, each = 10)
  proj <- shiny::reactiveVal(omicsCore::omics_project("dose", list(proteomics = inp)))
  shiny::testServer(diff_view_server, args = list(current_project = proj), {
    session$flushReact()
    expect_identical(continuous_cols(), "dose")
    session$setInputs(layer = "proteomics", design_mode = "continuous",
                      continuous_col = "dose", continuous_model = "linear",
                      method = "limma", rerun = 1)
    b <- diff_bundle()
    expect_identical(b$params$analysis_type, "continuous")
    expect_true(all(b$results$diff_result_df$adj_p_value[1:10] < 0.05))
  })
})

test_that("a failed run keeps the old result, labelled, and shows the engine's own sentence", {
  proj <- shiny::reactiveVal(tutorial_project())
  shiny::testServer(diff_view_server, args = list(current_project = proj), {
    session$setInputs(layer = "proteomics", group_col = "group", control = "Control",
                      case = c("TreatA", "TreatB"), method = "limma", rerun = 1)
    expect_false(is.null(diff_bundle()))
    session$setInputs(contrast_mode = "custom", custom_contrasts = "TreatA - Nope", rerun = 2)
    expect_false(is.null(diff_bundle()))
    expect_match(paste(as.character(output$header), collapse = ""), "the latest run failed",
                 fixed = TRUE)
    html <- paste(as.character(output$notices), collapse = "")
    expect_match(html, "is not a group", fixed = TRUE)
    expect_false(grepl("See the technical details below", html, fixed = TRUE))
  })
})

test_that("the full differential table downloads, every comparison included", {
  proj <- shiny::reactiveVal(tutorial_project())
  shiny::testServer(diff_view_server, args = list(current_project = proj), {
    session$setInputs(layer = "proteomics", group_col = "group", control = "Control",
                      case = c("TreatA", "TreatB"), method = "limma", rerun = 1)
    f <- output$download_csv
    df <- utils::read.csv(f)
    expect_setequal(unique(df$comparison), c("TreatA_vs_Control", "TreatB_vs_Control"))
    expect_identical(nrow(df), 2L * nrow(proj()$experiments$proteomics$expr_mat))
    expect_true("significant" %in% names(df))
  })
})

test_that("untested features are not counted as tested", {
  p <- tutorial_project()
  p$experiments$proteomics$expr_mat[1:5, 1:10] <- NA
  proj <- shiny::reactiveVal(p)
  shiny::testServer(diff_view_server, args = list(current_project = proj), {
    session$setInputs(layer = "proteomics", group_col = "group", control = "Control",
                      case = "TreatA", method = "limma", rerun = 1)
    html <- paste(as.character(output$stats), collapse = "")
    expect_match(html, "not testable", fixed = TRUE)
  })
})

test_that("an opened project shows its saved results, and never the demo pathways", {
  skip_if_not(has_pkg("clusterProfiler"))
  r4_store()
  p <- tutorial_project()
  d <- omicsCore::run_diff(p$experiments$proteomics, method = "limma", group_col = "group",
                           control_group = "Control", case_group = "TreatA")
  e <- omicsCore::run_enrichment(d, type = "ora", database = "hallmark")
  p$bundles <- list(diff = d, enrich = e)
  shiny::testServer(app_server, {
    session$setInputs(`project-load_tutorial` = 0)
    current_project(p)
    session$flushReact()
    expect_identical(omicsCore::diff_comparisons(diff_view$all_bundle()), "TreatA_vs_Control")
    expect_identical(enrich_view$bundle()$params$comparison, "TreatA_vs_Control")
    html <- paste(as.character(output$`enrich-header`), collapse = "")
    expect_false(grepl("demo data", html, fixed = TRUE))
  })
})

test_that("mouse gene names suggest the mouse gene sets", {
  expect_identical(guess_organism(c("Trp53", "Egfr", "Myc", "Il6", "Tnf", "Cd4", "Actb",
                                    "Gapdh", "Stat3", "Jun", "Fos")), "Mm")
  expect_identical(guess_organism(c("TP53", "EGFR", "MYC", "IL6", "TNF", "CD4", "ACTB",
                                    "GAPDH", "STAT3", "JUN", "FOS")), "Hs")
  expect_true(all(is_gene_symbol(c("TP53", "Trp53", "HLA-A"))))
  expect_false(any(is_gene_symbol(c("ENSG00000141510", "P04637", NA, ""))))
})

test_that("a step is done once its result exists and it has been visited", {
  p <- tutorial_project()
  p$bundles <- list(qc = 1)
  expect_identical(WORKFLOW_STEPS$id[workflow_progress(p)$next_step], "qc")
  p$visited_steps <- "qc"
  expect_identical(WORKFLOW_STEPS$id[workflow_progress(p)$next_step], "diff")
  legacy <- p
  legacy$visited_steps <- NULL
  expect_identical(WORKFLOW_STEPS$id[workflow_progress(legacy)$next_step], "diff")
})

test_that("a matrix file and a separate sample sheet are imported together", {
  r4_store()
  samp <- paste0("S", 1:6)
  set.seed(4)
  mat <- data.frame(Gene = paste0("GENE", 1:30),
                    matrix(round(stats::rnorm(180, 20), 2), 30, dimnames = list(NULL, samp)),
                    check.names = FALSE)
  mp <- withr::local_tempfile(fileext = ".csv")
  sp <- withr::local_tempfile(fileext = ".csv")
  utils::write.csv(mat, mp, row.names = FALSE)
  utils::write.csv(data.frame(sample = samp, group = rep(c("Ctrl", "KO"), each = 3)), sp,
                   row.names = FALSE)
  shiny::testServer(import_view_server, {
    session$setInputs(omics_type = "proteomics",
                      file = list(datapath = mp, name = "matrix.csv", size = file.size(mp)))
    expect_false("group" %in% names(parsed()$input$meta_df))
    html <- paste(as.character(output$confirm_design), collapse = "")
    expect_match(html, "no sample information", fixed = TRUE)
    session$setInputs(sample_file = list(datapath = sp, name = "samples.csv", size = file.size(sp)))
    expect_identical(parsed()$input$meta_df$group, rep(c("Ctrl", "KO"), each = 3))
  })
})

test_that("excluding a flagged sample replaces the layer and clears its results", {
  r4_store()
  p <- tutorial_project()
  p$experiments$proteomics$expr_mat[, 1] <- p$experiments$proteomics$expr_mat[, 1] +
    stats::rnorm(nrow(p$experiments$proteomics$expr_mat), 0, 3)
  p$bundles <- list(diff = omicsCore::run_diff(p$experiments$proteomics, method = "limma",
                                               group_col = "group", control_group = "Control",
                                               case_group = "TreatA"))
  proj <- shiny::reactiveVal(p)
  shiny::testServer(qc_view_server, args = list(current_project = proj), {
    session$setInputs(layer = "proteomics", missing_threshold = 0.5, outlier_method = "all",
                      impute_method = "none")
    flagged <- last_bundle()$results$qc_summary$outliers$flagged_samples
    skip_if(!length(flagged), "no sample flagged in this simulation")
    session$setInputs(exclude_flagged = 1)
    session$setInputs(confirm_exclude = 1)
    inp <- proj()$experiments$proteomics
    expect_false(any(flagged %in% colnames(inp$expr_mat)))
    expect_identical(inp$excluded_samples, flagged)
    expect_null(proj()$bundles$diff)
  })
})

test_that("the script downloads with the data files it reads", {
  r4_store()
  p <- tutorial_project()
  raw <- withr::local_tempfile(fileext = ".xlsx")
  writeLines("data", raw)
  p$experiments$proteomics$source_path <- raw
  shiny::testServer(report_view_server, args = list(current_project = shiny::reactiveVal(p)), {
    f <- output$download_bundle
    listing <- if (grepl("zip$", f)) utils::unzip(f, list = TRUE)$Name else utils::untar(f, list = TRUE)
    expect_true("analysis/analysis.R" %in% listing)
    expect_true(file.path("analysis", "raw", basename(raw)) %in% listing)
    expect_true("analysis/README.txt" %in% listing)
  })
})
