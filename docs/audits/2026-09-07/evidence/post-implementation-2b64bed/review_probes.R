suppressPackageStartupMessages({
  pkgload::load_all("packages/omicsCore", quiet = TRUE)
  pkgload::load_all("packages/omicsApp", quiet = TRUE)
  library(shiny)
})
options(shiny.allowoutputreads = TRUE)
emit <- function(key, value) cat(key, "=", if (length(value)) paste(value, collapse = " | ") else "<NULL>", "\n", sep = "")
case <- function(label, code) {
  cat("\nCASE ", label, "\n", sep = "")
  tryCatch(force(code), error = function(e) emit("PROBE_ERROR", conditionMessage(e)))
}
make_input <- function(seed = 9, type = "proteomics", assay = "normalized_intensity") {
  set.seed(seed)
  x <- matrix(rnorm(40 * 8, 12, 1), 40, 8,
              dimnames = list(paste0("g", 1:40), paste0("s", 1:8)))
  x[, 5:8] <- x[, 5:8] + 3
  md <- data.frame(group = rep(c("A", "B"), each = 4), row.names = colnames(x))
  fd <- data.frame(feature_id = rownames(x), feature_symbol = rownames(x), row.names = rownames(x))
  omicsCore::omics_input(x, md, fd, omics_type = type, assay_type = assay)
}
make_project <- function(name = "A", seed = 9) omicsCore::omics_project(name, list(prot = make_input(seed)))
make_diff <- function(p, tag = "prot") {
  e <- omicsCore::effective_input(p, tag)
  b <- omicsCore::run_diff(e$input, method = "limma", group_col = "group", control_group = "A", case_group = "B")
  omicsApp:::bind_result(b, p, "diff", tag, e$version)
}

case("explicit small raw intensities", {
  i <- make_input(assay = "raw_intensity")
  b <- omicsCore::run_diff(i, method = "auto", group_col = "group", control_group = "A", case_group = "B")
  d <- b$results$diff_result_df
  ix <- match(rownames(i$expr_mat), d$feature_id)
  expected_log <- rowMeans(log2(i$expr_mat[, 5:8])) - rowMeans(log2(i$expr_mat[, 1:4]))
  emit("RAW_INTENSITY_ACCEPTED", TRUE)
  emit("EFFECT_TYPE", unique(d$effect_type))
  emit("MAX_EFFECT_VS_LOG2_DIFFERENCE", max(abs(d$effect[ix] - expected_log)))
})

case("actual diff module late callback after switching projects", {
  pa <- make_project("A"); pb <- make_project("B")
  project <- reactiveVal(pa); generation <- reactiveVal(0L); pending <- NULL
  testthat::with_mocked_bindings({
    shiny::testServer(omicsApp:::diff_view_server,
      args = list(current_project = project, invalidate = generation), {
      session$setInputs(layer = "prot", group_col = "group", control = "A", case = "B", method = "limma", rerun = 1)
      session$flushReact()
      stopifnot(!is.null(pending))
      computed <- pending$func()
      project(pb); generation(1L); session$flushReact()
      emit("VIEW_EMPTY_AFTER_SWITCH", is.null(diff_bundle()))
      emit("LATE_SUCCESS_ADMITTED", pending$admissible("success")$ok)
      if (isTRUE(pending$admissible("success")$ok)) pending$on_success(computed)
      got <- diff_bundle()
      emit("A_RESULT_STAMPED_AS_B", identical(got$provenance_binding$project_id, pb$metadata$project_id))
      emit("WRONGLY_STAMPED_RESULT_VALID_IN_B", omicsCore::result_inputs_current(got$provenance_binding, pb))
      emit("LATE_ERROR_ADMITTED", pending$admissible("error")$ok)
      if (isTRUE(pending$admissible("error")$ok)) pending$on_error("late error from A")
      emit("CURRENT_ERROR", diff_error())
    })
  }, run_async = function(func, on_success, on_error, message, admissible, on_discard, ...) {
    pending <<- list(func = func, on_success = on_success, on_error = on_error, admissible = admissible)
  }, .package = "omicsApp")
})

case("hidden QC parameters survive loading and applying", {
  p <- make_project()
  saved <- omicsCore::run_qc(p$experiments$prot, impute_method = "none", outlier_sd_threshold = 7, sample_missing_threshold = 0.9)
  p <- omicsCore::set_preprocessing(p, "prot", saved$results$cleaned_input, saved$params)
  project <- reactiveVal(p)
  shiny::testServer(omicsApp:::qc_view_server, args = list(current_project = project), {
    session$setInputs(layer = "prot", missing_threshold = 0.5, outlier_method = "iqr", impute_method = "none")
    session$flushReact()
    session$setInputs(apply_to_downstream = 1)
    got <- session$getReturned()()
    emit("SAVED_OUTLIER_SD", saved$params$outlier_sd_threshold)
    emit("REAPPLIED_OUTLIER_SD", got$params$outlier_sd_threshold)
    emit("SAVED_SAMPLE_MISSING", saved$params$sample_missing_threshold)
    emit("REAPPLIED_SAMPLE_MISSING", got$params$sample_missing_threshold)
  })
})

case("QC application reaches both preprocessing and report bundle", {
  local({
    withr::local_envvar(OMICSAPP_DATA_DIR = tempfile("review-store-"))
    testthat::with_mocked_bindings({
      shiny::testServer(omicsApp:::app_server, {
        current_project(make_project()); session$flushReact()
        session$setInputs(`qc-layer` = "prot", `qc-missing_threshold` = 0.5, `qc-impute_method` = "none")
        session$setInputs(`qc-apply_to_downstream` = 1); session$flushReact()
        emit("PREPROCESSING_PRESENT", !is.null(current_project()$preprocessing$prot))
        emit("QC_BUNDLE_PRESENT", !is.null(current_project()$bundles$qc))
        emit("QC_APPLIED_MODULE_PRESENT", !is.null(qc_applied()))
      })
    }, run_async = function(...) NULL, .package = "omicsApp")
  })
})

case("saved result loads into actual analysis view", {
  local({
    store <- tempfile("review-store-"); dir.create(store)
    withr::local_envvar(OMICSAPP_DATA_DIR = store)
    p <- make_project(); b <- make_diff(p)
    p <- omicsCore::register_analysis(p, b); p$bundles <- list(diff = b)
    omicsApp:::store_save_project(p, "saved", dir = store)
    testthat::with_mocked_bindings({
      shiny::testServer(omicsApp:::app_server, {
        session$setInputs(`project-saved_pick` = "saved", `project-open_project` = 1)
        session$flushReact()
        emit("SAVED_PROJECT_DIFF_PRESENT", !is.null(current_project()$bundles$diff))
        emit("RESTORED_DIFF_VIEW_PRESENT", !is.null(diff_bundle()))
      })
    }, run_async = function(...) NULL, .package = "omicsApp")
  })
})

case("real project with no result does not display demo", {
  p <- reactiveVal(make_project())
  shiny::testServer(omicsApp:::enrich_view_server, args = list(current_project = p), {
    session$flushReact()
    emit("REAL_ENRICH_IS_DEMO", is_demo())
    emit("REAL_ENRICH_DEMO_ROWS", nrow(table_data()))
  })
  shiny::testServer(omicsApp:::integration_view_server, args = list(current_project = p), {
    session$flushReact()
    emit("REAL_INTEGRATION_IS_DEMO", is_demo())
    emit("REAL_INTEGRATION_AP_ROWS", nrow(ap_df()))
  })
})

case("actual integration keeps secondary source and exports executable script", {
  local({
    replay <- tempfile("integration-replay-"); dir.create(replay); dir.create(file.path(replay, "raw"))
    inp <- make_input()
    for (tag in c("prot", "prot2")) {
      path <- file.path(replay, "raw", paste0(tag, ".xlsx"))
      openxlsx::write.xlsx(list(expression = data.frame(feature_id = rownames(inp$expr_mat), inp$expr_mat, check.names = FALSE),
                               metadata = data.frame(sample_id = rownames(inp$meta_df), inp$meta_df),
                               features = inp$feature_df), path)
    }
    exps <- lapply(c("prot", "prot2"), function(tag) {
      path <- file.path(replay, "raw", paste0(tag, ".xlsx"))
      x <- omicsCore::read_omics(path, omics_type = "proteomics", assay_type = "normalized_intensity")$input
      x$source_path <- path; x
    }); names(exps) <- c("prot", "prot2")
    p <- omicsCore::omics_project("integration", exps)
    primary <- make_diff(p); p <- omicsCore::register_analysis(p, primary)
    project <- reactiveVal(p); upstream <- reactiveVal(primary); integrated <- NULL
    shiny::testServer(omicsApp:::integration_view_server, args = list(current_project = project, diff_bundle = upstream), {
      session$setInputs(rerun = 1); session$flushReact()
      integrated <<- integration_bundle()
      emit("INTEGRATION_ERROR", integration_error())
    })
    stopifnot(!is.null(integrated))
    emit("SECONDARY_ANALYSIS_IDS", unique(integrated$results$integration_df$analysis_id_b))
    emit("INTEGRATION_UPSTREAM_COUNT", length(integrated$provenance_binding$upstream_ids))
    emit("INTEGRATION_PARAM_SIGNATURE_MATCH", identical(integrated$provenance_binding$params_signature, omicsCore::params_signature(integrated$params)))
    p <- omicsCore::register_analysis(p, integrated)
    p$bundles <- list(diff = primary, integration = integrated)
    script <- file.path(replay, "replay.R"); omicsCore::export_script(p, script, include_plots = FALSE)
    replay_env <- new.env(parent = globalenv())
    err <- tryCatch(withr::with_dir(replay, {sys.source(script, replay_env); "NONE"}), error = function(e) conditionMessage(e))
    emit("INTEGRATION_REPLAY_ERROR", err)
    file.copy(script, "integration-replay-generated.R", overwrite = TRUE)
    emit("REPLAY_SCRIPT_PATH", "integration-replay-generated.R")
  })
})

case("unrelated project update resets QC draft controls", {
  p <- make_project()
  b <- omicsCore::run_qc(p$experiments$prot, impute_method = "min")
  p <- omicsCore::set_preprocessing(p, "prot", b$results$cleaned_input, b$params)
  project <- reactiveVal(p); updates <- list()
  testthat::with_mocked_bindings({
    shiny::testServer(omicsApp:::qc_view_server, args = list(current_project = project), {
      session$setInputs(layer = "prot", impute_method = "none")
      session$flushReact(); updates <<- list()
      same <- project(); same$name <- "unrelated title change"; project(same)
      session$flushReact()
      emit("DRAFT_IMPUTE_BEFORE_BROWSER_UPDATE", input$impute_method)
      emit("RESTORE_MESSAGES_ON_UNRELATED_CHANGE", length(updates))
      emit("RESTORE_SELECTED", vapply(updates, function(u) u$selected, character(1)))
    })
  }, updateSelectInput = function(session, inputId, label = NULL, choices = NULL, selected = NULL) {
    updates[[length(updates) + 1L]] <<- list(id = inputId, selected = selected)
  }, .package = "shiny")
})

case("same workbook with changed sheet roles can be reimported", {
  local({
    withr::local_envvar(OMICSAPP_DATA_DIR = tempfile("roles-store-"))
    i <- make_input(); x <- data.frame(feature_id = rownames(i$expr_mat), i$expr_mat, check.names = FALSE)
    y <- x; y[, -1] <- y[, -1] + 5
    path <- tempfile(fileext = ".xlsx")
    openxlsx::write.xlsx(list(ExpressionA = x, ExpressionB = y,
       Metadata = data.frame(sample_id = rownames(i$meta_df), i$meta_df), Features = i$feature_df), path)
    project <- reactiveVal(NULL)
    shiny::testServer(omicsApp:::import_view_server, args = list(current_project = project), {
      session$setInputs(file = list(datapath = path, name = "roles.xlsx", size = file.info(path)$size), omics_type = "proteomics", normalize = FALSE)
      role_overrides(c(ExpressionA = "matrix", ExpressionB = "unknown", Metadata = "metadata", Features = "feature_annot")); do_parse()
      session$setInputs(confirm = 1)
      first <- confirmed_input()
      project(omicsCore::omics_project("roles", list(proteomics = first)))
      role_overrides(c(ExpressionA = "unknown", ExpressionB = "matrix", Metadata = "metadata", Features = "feature_annot")); do_parse()
      emit("ROLE_CHANGE_CANDIDATE_MATRIX_DIFF", max(abs(parsed()$input$expr_mat - first$expr_mat)))
      session$setInputs(confirm = 2)
      emit("ROLE_CHANGE_NEW_INPUT_COMMITTED", !is.null(confirmed_input()))
      emit("ROLE_CHANGE_PROJECT_MATRIX_DIFF", max(abs(project()$experiments$proteomics$expr_mat - first$expr_mat)))
      emit("ROLE_CHANGE_FINGERPRINT_STILL_MATCHES", omicsApp:::fingerprints_match(first, stamp_fingerprint(parsed()$input)))
    })
  })
})
