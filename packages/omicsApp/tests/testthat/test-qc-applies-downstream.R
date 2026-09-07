# QC's filtering and imputation reaching the analyses -- and not reaching
# them before the user asks.
#
# Both halves matter. The preview runs from the moment the module is
# constructed, on every control change, whether or not anyone has opened the
# view; wiring that straight to the analyses would mean the default cutoff
# silently filtered every project. And until it was wired at all, a user who
# removed high-missing features and imputed the rest watched the differential
# run on the raw import while the QC panel showed the cleaned numbers.

qcd_input <- function(omics_type = "proteomics", n_features = 30, n_na = 12,
                      seed = 4) {
  set.seed(seed)
  n <- 8
  e <- matrix(rnorm(n_features * n, 20, 1), nrow = n_features,
              dimnames = list(paste0("g", seq_len(n_features)),
                              paste0("s", seq_len(n))))
  # Concentrated in the first two features, so a cutoff can single them out.
  e[1, seq_len(min(n_na, n))] <- NA
  e[2, seq_len(min(n_na %/% 2, n))] <- NA
  md <- data.frame(group = rep(c("A", "B"), each = n / 2),
                   row.names = colnames(e), stringsAsFactors = FALSE)
  fd <- data.frame(feature_id = rownames(e), feature_symbol = rownames(e),
                   row.names = rownames(e), stringsAsFactors = FALSE)
  omicsCore::omics_input(
    e, md, fd, omics_type = omics_type,
    # Log-scale on both sides: the point here is QC reaching the analyses,
    # not the scale gate, and limma refuses a linear layer.
    assay_type = if (omics_type == "rnaseq") "logcpm" else "normalized_intensity")
}

qcd_project <- function(...) {
  exps <- list(...)
  if (!length(exps)) exps <- list(prot = qcd_input())
  omicsCore::omics_project("qcd", exps)
}

apply_qc <- function(project, tag, ..., click = 1L) {
  proj <- shiny::reactiveVal(project)
  out <- NULL
  shiny::testServer(qc_view_server, args = list(current_project = proj), {
    session$setInputs(layer = tag)
    session$setInputs(...)
    session$flushReact()
    session$setInputs(apply_to_downstream = click)
    session$flushReact()
    out <<- session$getReturned()()
  })
  out
}

# --- the preview is not the answer --------------------------------------

test_that("previewing changes nothing downstream", {
  proj <- shiny::reactiveVal(qcd_project())
  shiny::testServer(qc_view_server, args = list(current_project = proj), {
    session$setInputs(layer = "prot", missing_threshold = 0.1)
    session$flushReact()
    # The preview ran -- it always does -- and produced a cleaned layer.
    expect_false(is.null(last_bundle()))
    # Nothing left the module.
    expect_null(session$getReturned()())
  })
})

test_that("the module hands back the applied layer, not the preview", {
  out <- apply_qc(qcd_project(), "prot",
                  missing_threshold = 0.1, impute_method = "min")
  expect_equal(out$tag, "prot")
  expect_s3_class(out$bundle, "analysis_bundle")
  expect_true(nzchar(out$base_version))
  expect_equal(out$params$impute_method, "min")
})

test_that("a preview that failed cannot be applied", {
  # A value the slider cannot produce but the websocket can, which is the
  # premise of every hostile-input case here. The old behaviour kept the
  # last successful bundle in place, so pressing Apply after a failure
  # applied settings that had never run.
  proj <- shiny::reactiveVal(qcd_project())
  shiny::testServer(qc_view_server, args = list(current_project = proj), {
    session$setInputs(layer = "prot", missing_threshold = 0.5)
    session$flushReact()
    expect_true(isTRUE(preview()$ok))

    session$setInputs(missing_threshold = -1)
    session$flushReact()
    expect_false(isTRUE(preview()$ok))

    session$setInputs(apply_to_downstream = 1)
    session$flushReact()
    expect_null(session$getReturned()())
  })
})

test_that("a preview from another layer cannot be applied to this one", {
  # Preview layer A, switch to B, B's preview fails, press Apply. The
  # signature is what refuses it: without one, the apply handler stamped
  # the *current* layer's tag and version onto A's result.
  proj <- shiny::reactiveVal(qcd_project(prot = qcd_input(),
                                         rna = qcd_input("rnaseq")))
  shiny::testServer(qc_view_server, args = list(current_project = proj), {
    session$setInputs(layer = "prot", missing_threshold = 0.5)
    session$flushReact()
    from_a <- preview()$bundle
    expect_false(is.null(from_a))

    session$setInputs(layer = "rna", missing_threshold = -1)
    session$flushReact()
    expect_false(isTRUE(preview()$ok))

    session$setInputs(apply_to_downstream = 1)
    session$flushReact()
    expect_null(session$getReturned()())
  })
})

# --- QC recomputes from the import, not from its own output -------------

test_that("loosening the cutoff after applying recovers the features", {
  # The reason QC reads baseline_input() rather than effective_input(). Run
  # over already-filtered values, a looser cutoff has nothing to give back.
  p <- qcd_project()
  tight <- apply_qc(p, "prot", missing_threshold = 0.1)
  p2 <- omicsCore::set_preprocessing(p, "prot",
                                     tight$bundle$results$cleaned_input,
                                     tight$params, tight$base_version)
  n_tight <- nrow(omicsCore::effective_input(p2, "prot")$input$expr_mat)
  expect_lt(n_tight, nrow(p$experiments$prot$expr_mat))

  loose <- apply_qc(p2, "prot", missing_threshold = 0.9)
  n_loose <- nrow(loose$bundle$results$cleaned_input$expr_mat)
  expect_gt(n_loose, n_tight)
  # And it is exactly what the same settings give when run against the
  # import directly -- the applied layer left no trace on the recomputation.
  from_import <- omicsCore::run_qc(p$experiments$prot, missing_threshold = 0.9)
  expect_equal(n_loose, nrow(from_import$results$cleaned_input$expr_mat))
})

test_that("applying the same settings twice does not compound", {
  p <- qcd_project()
  once <- apply_qc(p, "prot", missing_threshold = 0.1, impute_method = "min")
  p2 <- omicsCore::set_preprocessing(p, "prot",
                                     once$bundle$results$cleaned_input,
                                     once$params, once$base_version)
  twice <- apply_qc(p2, "prot", missing_threshold = 0.1, impute_method = "min")
  expect_equal(nrow(twice$bundle$results$cleaned_input$expr_mat),
               nrow(once$bundle$results$cleaned_input$expr_mat))
})

# --- the controls describe what is applied ------------------------------

test_that("arriving at a layer shows the settings it was applied with", {
  # Opening a project saved with `min` while `knn` happened to be selected
  # showed `knn`, so the page described settings the analyses were not using.
  p <- qcd_project()
  applied <- apply_qc(p, "prot", missing_threshold = 0.15, impute_method = "min")
  saved <- omicsCore::set_preprocessing(p, "prot",
                                        applied$bundle$results$cleaned_input,
                                        applied$params, applied$base_version)

  proj <- shiny::reactiveVal(NULL)
  shiny::testServer(qc_view_server, args = list(current_project = proj), {
    session$setInputs(layer = "prot", missing_threshold = 0.9,
                      impute_method = "MinProb")
    session$flushReact()
    proj(saved)                     # the project arrives with QC applied
    session$flushReact()
    expect_equal(applied_pp()$params$missing_threshold, 0.15)
    expect_equal(applied_pp()$params$impute_method, "min")
  })
})

test_that("preprocessing left behind by a replaced layer is not offered back", {
  p <- qcd_project()
  applied <- apply_qc(p, "prot", missing_threshold = 0.15)
  saved <- omicsCore::set_preprocessing(p, "prot",
                                        applied$bundle$results$cleaned_input,
                                        applied$params, applied$base_version)
  # A different layer under the same tag. add_experiment() stamps its
  # version, which is what makes the stale preprocessing detectable.
  saved <- omicsCore::remove_experiment(saved, "prot")
  saved <- omicsCore::add_experiment(saved, "prot", qcd_input(seed = 99))

  proj <- shiny::reactiveVal(saved)
  shiny::testServer(qc_view_server, args = list(current_project = proj), {
    session$setInputs(layer = "prot")
    session$flushReact()
    expect_null(applied_pp())
  })
})

# --- what the analyses actually receive ---------------------------------

test_that("the differential runs on the applied layer", {
  p <- qcd_project()
  applied <- apply_qc(p, "prot", missing_threshold = 0.1, impute_method = "min")
  p2 <- omicsCore::set_preprocessing(p, "prot",
                                     applied$bundle$results$cleaned_input,
                                     applied$params, applied$base_version)

  proj <- shiny::reactiveVal(p2)
  shiny::testServer(diff_view_server, args = list(current_project = proj), {
    session$setInputs(layer = "prot")
    session$flushReact()
    a <- active()
    expect_equal(a$source, "qc_cleaned")
    expect_equal(nrow(a$input$expr_mat),
                 nrow(applied$bundle$results$cleaned_input$expr_mat))
    expect_false(anyNA(a$input$expr_mat))
    # The raw layer had both, so this is the difference the QC page promised.
    expect_lt(nrow(a$input$expr_mat), nrow(p$experiments$prot$expr_mat))
    expect_true(anyNA(p$experiments$prot$expr_mat))
  })
})

test_that("without an applied layer the differential uses the import", {
  proj <- shiny::reactiveVal(qcd_project())
  shiny::testServer(diff_view_server, args = list(current_project = proj), {
    session$setInputs(layer = "prot")
    session$flushReact()
    expect_equal(active()$source, "imported")
  })
})

test_that("each layer keeps its own applied preprocessing", {
  p <- qcd_project(prot = qcd_input(), rna = qcd_input("rnaseq"))
  a <- apply_qc(p, "prot", missing_threshold = 0.1)
  p <- omicsCore::set_preprocessing(p, "prot", a$bundle$results$cleaned_input,
                                    a$params, a$base_version)
  b <- apply_qc(p, "rna", missing_threshold = 0.9)
  p <- omicsCore::set_preprocessing(p, "rna", b$bundle$results$cleaned_input,
                                    b$params, b$base_version)

  # A single global slot let the second apply replace the first.
  eff_p <- omicsCore::effective_input(p, "prot")
  eff_r <- omicsCore::effective_input(p, "rna")
  expect_equal(eff_p$source, "qc_cleaned")
  expect_equal(eff_r$source, "qc_cleaned")
  expect_lt(nrow(eff_p$input$expr_mat), nrow(eff_r$input$expr_mat))
})

test_that("the integration's second layer uses the applied version too", {
  # The second differential is computed inside the integration module from
  # proj$experiments, so a concordance compared a cleaned first layer
  # against an uncleaned second one.
  p <- qcd_project(prot = qcd_input(), rna = qcd_input("rnaseq"))
  b <- apply_qc(p, "rna", missing_threshold = 0.1, impute_method = "min")
  p <- omicsCore::set_preprocessing(p, "rna", b$bundle$results$cleaned_input,
                                    b$params, b$base_version)

  primary <- omicsCore::run_diff(
    omicsCore::effective_input(p, "prot")$input,
    method = "limma", analysis_type = "group", group_col = "group",
    control_group = "A", case_group = "B")

  proj <- shiny::reactiveVal(p)
  db <- shiny::reactiveVal(primary)
  shiny::testServer(integration_view_server,
                    args = list(current_project = proj, diff_bundle = db), {
    session$flushReact()
    info <- can_run()
    skip_if_not(isTRUE(info$ok), "fixture did not produce a runnable pair")
    expect_equal(info$secondary_source, "qc_cleaned")
    expect_equal(nrow(info$secondary$expr_mat),
                 nrow(b$bundle$results$cleaned_input$expr_mat))
  })
})

# --- parameters the interface does not show ------------------------------

# The decision was that QC parameters without a control keep their saved
# values. The preview built its argument list from the three visible
# controls alone, so re-applying a restored project reset every other one
# to the default -- a layer saved with outlier_sd_threshold = 7 came back
# at 3, and sample_missing_threshold vanished. The cleaned matrix restoring
# correctly says nothing about this: the loss happens on the next apply.
test_that("re-applying keeps the parameters that have no control", {
  p <- qcd_project()
  saved_params <- list(missing_threshold = 0.4, outlier_method = "iqr",
                       impute_method = "min",
                       outlier_sd_threshold = 7,
                       sample_missing_threshold = 0.9)
  b <- do.call(omicsCore::run_qc, c(list(input = p$experiments$prot),
                                    saved_params))
  p <- omicsCore::set_preprocessing(p, "prot", b$results$cleaned_input,
                                    saved_params)

  again <- apply_qc(p, "prot")
  expect_equal(again$params$outlier_sd_threshold, 7)
  expect_equal(again$params$sample_missing_threshold, 0.9)
  # And the visible ones still come from the controls.
  expect_equal(again$params$missing_threshold, 0.4)
})

test_that("a visible control still overrides what was saved", {
  p <- qcd_project()
  saved <- list(missing_threshold = 0.4, outlier_method = "iqr",
                impute_method = "min", outlier_sd_threshold = 7)
  b <- do.call(omicsCore::run_qc, c(list(input = p$experiments$prot), saved))
  p <- omicsCore::set_preprocessing(p, "prot", b$results$cleaned_input, saved)

  changed <- apply_qc(p, "prot", missing_threshold = 0.8)
  expect_equal(changed$params$missing_threshold, 0.8)
  # Overriding one does not reset the others.
  expect_equal(changed$params$outlier_sd_threshold, 7)
})

test_that("an unrelated project update does not overwrite the draft", {
  # The restore observer watched the project object, and every analysis
  # write-back updates it -- so finishing a differential, or renaming the
  # project, pushed the applied settings back over what the user was in
  # the middle of choosing.
  p <- qcd_project()
  saved <- list(missing_threshold = 0.4, outlier_method = "iqr",
                impute_method = "min")
  b <- do.call(omicsCore::run_qc, c(list(input = p$experiments$prot), saved))
  p <- omicsCore::set_preprocessing(p, "prot", b$results$cleaned_input, saved)

  proj <- shiny::reactiveVal(p)
  act <- shiny::reactiveVal(1L)
  shiny::testServer(qc_view_server,
                    args = list(current_project = proj, activation = act), {
    session$setInputs(layer = "prot")
    session$flushReact()
    # The user is part way through a different choice.
    session$setInputs(missing_threshold = 0.05)
    session$flushReact()
    expect_equal(collect_qc_params()$missing_threshold, 0.05)

    # Something else finishes and updates the project.
    renamed <- p
    renamed$name <- "renamed mid-edit"
    proj(renamed)
    session$flushReact()
    expect_equal(collect_qc_params()$missing_threshold, 0.05)
  })
})

test_that("applying does not clear the result it just produced", {
  # Applying bumps the generation, which is how the downstream views are
  # told to let go -- and the QC view watched the same signal, so it
  # cleared its own output in the same flush and the observer that files
  # results into the project found nothing.
  skip_if_not_installed("qs2")
  store <- withr::local_tempdir()
  withr::local_envvar(OMICSAPP_DATA_DIR = store)
  p <- qcd_project()
  proj <- shiny::reactiveVal(p)
  act <- shiny::reactiveVal(1L)
  gen <- shiny::reactiveVal(0L)

  shiny::testServer(qc_view_server,
                    args = list(current_project = proj, invalidate = gen,
                                activation = act), {
    session$setInputs(layer = "prot", missing_threshold = 0.1)
    session$flushReact()
    session$setInputs(apply_to_downstream = 1)
    session$flushReact()
    expect_false(is.null(session$getReturned()()))

    # What app_server does next.
    gen(1L)
    session$flushReact()
    expect_false(is.null(session$getReturned()()),
                 info = "the apply cleared its own output")

    # A genuine layer replacement still clears it.
    gen(2L)
    session$flushReact()
    expect_null(session$getReturned()())
  })
})

test_that("the QC bundle reaches the project, not just the cleaned data", {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("readxl")
  skip_if_not_installed("withr")
  store <- file.path(withr::local_tempdir(), "store")
  dir.create(store, recursive = TRUE)
  withr::local_envvar(OMICSAPP_DATA_DIR = store)
  xlsx <- tempfile(fileext = ".xlsx")
  withr::defer(unlink(xlsx))
  write_tiny_omics_xlsx(xlsx, n_features = 40L, n_samples = 8L, seed = 2)
  shiny::testServer(app_server, {
    suppressWarnings(session$setInputs(
      `import-omics_type` = "proteomics",
      `import-file` = list(datapath = xlsx, name = basename(xlsx),
                           size = file.info(xlsx)$size)))
    session$setInputs(`import-confirm` = 1)
    session$flushReact()
    session$setInputs(`qc-missing_threshold` = 0.3,
                      `qc-apply_to_downstream` = 1)
    session$flushReact()
    proj <- current_project()
    expect_false(is.null(proj$preprocessing[[names(proj$experiments)[[1]]]]),
                 info = "cleaned data did not reach the project")
    # The diagnostics and the QC parameters are a result too, and a report
    # that cannot show them is describing an analysis it did not keep.
    expect_true("qc" %in% names(proj$bundles),
                info = "the QC bundle did not reach the project")
  })
})
