# Every control, driven through every value it offers.
#
# A control that is rendered and does nothing reads as a broken app, and
# it is the failure mode this app has actually shipped: an outlier-method
# radio that did not re-run, a Re-run that greyed the page, download
# buttons that produced no file and no message. None of those errored.
#
# So these assert the thing a user checks -- something changed, or
# something rendered -- rather than that a function was reachable.

sm_input <- function(omics_type = "proteomics", n_feat = 40L) {
  ids <- paste0("S", 1:8)
  set.seed(11)
  m <- matrix(as.numeric(stats::rpois(n_feat * 8, 300)), nrow = n_feat,
              dimnames = list(paste0("F", seq_len(n_feat)), ids))
  if (omics_type == "proteomics") m[sample(length(m), 30)] <- NA
  omicsCore::omics_input(
    m,
    data.frame(sample_id = ids,
               condition = rep(c("G1", "G2"), each = 4),
               donor = paste0("D", 1:8),
               row.names = ids, stringsAsFactors = FALSE),
    data.frame(feature_id = rownames(m),
               feature_symbol = rownames(m),
               row.names = rownames(m), stringsAsFactors = FALSE),
    omics_type = omics_type,
    assay_type = if (omics_type == "rnaseq") "raw_count" else "raw_intensity")
}

sm_project <- function() {
  omicsCore::omics_project("smoke",
                           list(prot = sm_input("proteomics"),
                                rna  = sm_input("rnaseq")))
}

rendered <- function(x) {
  !is.null(x) && nzchar(paste(as.character(x), collapse = ""))
}

# testServer hands back a renderUI output as `list(html =, deps =)`, so
# as.character() on it is a two-element vector and a fixed match against
# the whole thing reads whichever element came first.
ui_html <- function(x) {
  if (is.list(x) && !is.null(x$html)) as.character(x$html) else
    paste(as.character(x), collapse = "")
}

# ---- QC ---------------------------------------------------------------

test_that("every QC control changes what QC computes", {
  proj <- shiny::reactiveVal(sm_project())
  shiny::testServer(qc_view_server, args = list(current_project = proj), {
    session$flushReact()

    for (m in c("iqr", "pca", "connectivity")) {
      session$setInputs(outlier_method = m)
      session$flushReact()
      b <- last_bundle()
      expect_false(is.null(b), info = m)
      # The bundle records the method it ran, so a radio that did not
      # re-run is visible here rather than only on the plot.
      expect_identical(b$params$outlier_method, m)
    }

    # The slider is a filter: a stricter cutoff keeps no more features.
    session$setInputs(missing_threshold = 1)
    session$flushReact()
    loose <- last_bundle()$input_info$n_features_out
    session$setInputs(missing_threshold = 0.01)
    session$flushReact()
    tight <- last_bundle()$input_info$n_features_out
    expect_lte(tight, loose)

    for (v in c("depth", "missing")) {
      session$setInputs(quality_view = v)
      session$flushReact()
      expect_identical(quality_view(), v)
    }

    for (tag in c("rna", "prot")) {
      session$setInputs(layer = tag)
      session$flushReact()
      expect_identical(active()$tag, tag)
    }
  })
})

# ---- Differential -----------------------------------------------------

test_that("the Differential controls reach the analysis", {
  proj <- shiny::reactiveVal(sm_project())
  shiny::testServer(diff_view_server, args = list(current_project = proj), {
    session$flushReact()
    out <- session$getReturned()

    session$setInputs(layer = "prot", group_col = "condition",
                      control_group = "G1", case_group = "G2",
                      p_kind = "raw", fdr_cut = 0.05, fc_cut = 0.58)
    # The cutoffs are debounced by 250ms so that dragging a slider does
    # not re-filter on every pixel. testServer's clock does not advance
    # on its own, so without this the reactive still holds its default
    # and the test reads as a control that does nothing.
    session$elapse(400)

    thr <- out$thresholds()
    expect_identical(thr$p_preference, "raw")
    expect_equal(thr$p_cutoff, 0.05)
    expect_equal(thr$effect_cutoff, 0.58)

    session$setInputs(p_kind = "adj", fdr_cut = 0.2)
    session$elapse(400)
    thr <- out$thresholds()
    expect_identical(thr$p_preference, "adjusted")
    expect_equal(thr$p_cutoff, 0.2)

    # The layer the result belongs to travels with it.
    expect_identical(out$layer(), "prot")
    session$setInputs(layer = "rna")
    session$flushReact()
    expect_identical(out$layer(), "rna")
  })
})

# ---- Enrichment -------------------------------------------------------

test_that("the Enrichment display controls change the table and the plot", {
  shiny::testServer(
    enrich_view_server,
    args = list(diff_bundle = shiny::reactiveVal(example_diff_bundle())), {
      session$setInputs(type = "ora", database = "hallmark",
                        direction = "both",
                        show_p = "adjusted", show_cutoff = 0.05)

      n_strict <- nrow(omicsCore::filter_enrich_results(
        table_data(), p_cutoff = show_cutoff(), p_preference = show_p()))

      session$setInputs(show_p = "raw", show_cutoff = 0.5)
      n_loose <- nrow(omicsCore::filter_enrich_results(
        table_data(), p_cutoff = show_cutoff(), p_preference = show_p()))

      # Raw p at 0.5 cannot keep fewer pathways than adjusted p at 0.05.
      expect_gte(n_loose, n_strict)

      # The download is the whole table regardless of either control.
      written <- utils::read.csv(output$download_table, check.names = FALSE)
      expect_equal(nrow(written), nrow(table_data()))
    })
})

# ---- Report -----------------------------------------------------------

test_that("the export button is live exactly when it can work", {
  skip_if_not_installed("openxlsx")

  # No project: the handler stops on req(proj), silently, as req() does.
  # A live-looking button that produces no file and no message is
  # indistinguishable from a broken one -- which is what it looked like.
  shiny::testServer(report_view_server,
                    args = list(current_project = shiny::reactiveVal(NULL)), {
    html <- ui_html(output$header)
    expect_match(html, "disabled")
    expect_match(ui_html(output$notices), "Import data first")
  })

  # A project with no analyses is the same case, and used to be missed:
  # the workbook would hold nothing, so there is nothing to download.
  shiny::testServer(report_view_server,
                    args = list(current_project = shiny::reactiveVal(sm_project())), {
    expect_match(ui_html(output$header), "pointer-events:none")
    expect_match(ui_html(output$sheet_card), "Nothing to export yet")
  })

  withr::with_options(list(shiny.allowoutputreads = TRUE), {
    proj <- sm_project()
    proj$bundles <- list(qc = example_qc_bundle(), diff = example_diff_bundle())
    shiny::testServer(report_view_server,
                      args = list(current_project = shiny::reactiveVal(proj)), {
      expect_no_match(ui_html(output$header), "pointer-events:none")
      # The card lists what the file will hold, from the same source the
      # download reads, so the two cannot describe different workbooks.
      card <- ui_html(output$sheet_card)
      expect_match(card, "Differential", fixed = TRUE)
      expect_match(card, "Analysis parameters", fixed = TRUE)
    })
  })
})

# The Enrich view runs one combination at a time, so everything ORA and
# GSEA have to say about a differential took fourteen visits to it. The
# grid is one press, and what it produces has to reach the workbook --
# through the project, like every other result, not through a private
# channel between two modules.
test_that("the enrichment grid runs from the Report view and reaches the sheets", {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("msigdbr")

  proj <- sm_project()
  proj$bundles <- list(diff = example_diff_bundle())
  pv <- shiny::reactiveVal(proj)

  shiny::testServer(report_view_server, args = list(current_project = pv), {
    session$flushReact()
    expect_null(grid_bundle())

    session$setInputs(run_grid = 1)
    session$flushReact()

    g <- grid_bundle()
    expect_false(is.null(g))
    expect_identical(g$analysis_name, "run_enrichment_grid")
    # Both tests attempted, whatever came back
    expect_setequal(unique(g$results$grid_log$type), c("ora", "gsea"))

    # And the workbook picks it up once it is on the project. app_server()
    # does the filing; here it is done by hand to keep the test to one
    # module.
    p2 <- proj
    p2$bundles$enrich_all <- g
    sheets <- names(omicsCore::workbook_sheets(p2))
    expect_true("Enrichment log" %in% sheets)
  })
})

test_that("the grid button is refused without a differential to enrich", {
  proj <- sm_project()            # layers, no analyses
  shiny::testServer(report_view_server,
                    args = list(current_project = shiny::reactiveVal(proj)), {
    session$flushReact()
    expect_match(ui_html(output$grid_card), "pointer-events:none")

    # And pressing it anyway says why rather than failing in a worker.
    session$setInputs(run_grid = 1)
    session$flushReact()
    expect_match(grid_error(), "differential")
  })
})

test_that("the workbook comes out with a sheet per analysis", {
  skip_if_not_installed("openxlsx")

  proj <- sm_project()
  proj$bundles <- list(qc = example_qc_bundle(), diff = example_diff_bundle())

  out <- withr::local_tempfile(fileext = ".xlsx")
  omicsCore::export_workbook(proj, out, overwrite = TRUE)
  sheets <- openxlsx::getSheetNames(out)
  expect_true(all(c("Analysis parameters", "Differential") %in% sheets))
  # and the rows survive the trip, which a workbook of headers would not
  back <- openxlsx::read.xlsx(out, sheet = "Differential")
  expect_equal(nrow(back),
               nrow(proj$bundles$diff$results$diff_result_df))
})

# ---- Project ----------------------------------------------------------

test_that("the Projects controls act on the project", {
  proj <- shiny::reactiveVal(sm_project())
  shiny::testServer(project_view_server,
                    args = list(current_project = proj), {
    session$flushReact()
    # View hands a layer to the other views.
    session$setInputs(view_layer_2 = 1)
    session$flushReact()

    session$setInputs(drop_layer_1 = 1)
    session$flushReact()
    expect_identical(pending_drop(), "prot")
    session$setInputs(confirm_drop_layer = 1)
    session$flushReact()
    expect_identical(names(current_project()$experiments), "rna")
  })
})

# The session restores the last project on arrival, which is what someone
# coming back wants and left someone starting a new study with no way to
# say so -- importing a file replaces a layer of the project they are
# already in. This is the way out.
test_that("Start a new project clears the session, once confirmed", {
  proj <- shiny::reactiveVal(sm_project())
  activated <- 0L
  shiny::testServer(
    project_view_server,
    args = list(current_project = proj,
                on_activate = function() activated <<- activated + 1L), {
    session$flushReact()

    # Asking is not doing: the modal has to be confirmed.
    session$setInputs(new_project = 1)
    session$flushReact()
    expect_false(is.null(current_project()))

    session$setInputs(confirm_new_project = 1)
    session$flushReact()
    expect_null(current_project())
    # The views are told to let go before the project disappears.
    expect_gt(activated, 0L)
  })
})

test_that("the new-project button is offered only with a project loaded", {
  # With nothing loaded this is already where the button would take you,
  # and a button that does nothing reads as a broken one.
  shiny::testServer(project_view_server,
                    args = list(current_project = shiny::reactiveVal(NULL)), {
    session$flushReact()
    expect_no_match(ui_html(output$storage), "Start a new project",
                    fixed = TRUE)
  })
  shiny::testServer(project_view_server,
                    args = list(current_project = shiny::reactiveVal(sm_project())), {
    session$flushReact()
    expect_match(ui_html(output$storage), "Start a new project",
                 fixed = TRUE)
  })
})

# ---- Integration ------------------------------------------------------

test_that("the Integration pairing control writes the pairing", {
  proj <- shiny::reactiveVal(sm_project())
  shiny::testServer(integration_view_server,
                    args = list(current_project = proj), {
    # Both layers carry a donor column, so nothing needs confirming.
    expect_identical(pairing()$source, "donor")
    expect_equal(nrow(pairing()$pairs), 8L)
    expect_null(output$pairing_action)
    expect_true(rendered(output$pairing_note))
  })
})

# ---- QC imputation ----------------------------------------------------
# Offered for proteomics, withheld for counts. A missing intensity in DIA
# usually means "below the detection limit", which is information worth
# filling in deliberately; a missing count does not mean that, and
# imputing counts feeds DESeq2 numbers its model never saw.

test_that("the imputation control reaches run_qc and changes the result", {
  proj <- shiny::reactiveVal(sm_project())
  shiny::testServer(qc_view_server, args = list(current_project = proj), {
    session$setInputs(layer = "prot")
    session$flushReact()

    session$setInputs(impute_method = "none")
    session$flushReact()
    b <- last_bundle()
    expect_identical(b$params$impute_method, "none")
    expect_true(anyNA(b$results$cleaned_input$expr_mat))

    session$setInputs(impute_method = "MinDet")
    session$flushReact()
    b <- last_bundle()
    expect_identical(b$params$impute_method, "MinDet")
    # The point of the control: the NAs are gone afterwards.
    expect_false(anyNA(b$results$cleaned_input$expr_mat))
  })
})

test_that("counts are never imputed, even with a method left selected", {
  proj <- shiny::reactiveVal(sm_project())
  shiny::testServer(qc_view_server, args = list(current_project = proj), {
    session$setInputs(layer = "prot", impute_method = "MinDet")
    session$flushReact()
    expect_true(rendered(output$ui_impute))

    # Shiny keeps an input's value when its control is removed, so the
    # hidden `half_min` would otherwise follow the user to a counts
    # layer and impute it with a control they can no longer see.
    session$setInputs(layer = "rna")
    session$flushReact()
    expect_null(output$ui_impute)
    expect_identical(last_bundle()$params$impute_method, "none")
  })
})

# The low-expression filter is the mirror of imputation: offered on raw
# counts, withheld everywhere else. Both halves matter, and the second one
# is the half that bites -- Shiny keeps a removed input's value, so a count
# threshold set here would otherwise follow the user to an intensity layer,
# where 10 is not a low number.

test_that("the low-expression control reaches run_qc and changes the result", {
  # A mixture rather than sm_input()'s uniform Poisson(300): a threshold
  # that drops every feature is refused by run_qc() -- rightly -- and the
  # preview then keeps the previous bundle, so a test built that way
  # watches nothing change and cannot tell why.
  quiet <- sm_input("rnaseq")
  quiet$expr_mat[1:10, ] <- 0
  proj <- shiny::reactiveVal(
    omicsCore::omics_project("smoke", list(prot = sm_input("proteomics"),
                                           rna = quiet)))
  shiny::testServer(qc_view_server, args = list(current_project = proj), {
    session$setInputs(layer = "rna")
    session$flushReact()
    expect_true(rendered(output$ui_low_expression))

    # Defaults: the ten silent features go, the other thirty stay.
    b <- last_bundle()
    expect_identical(b$params$min_count, 10)
    expect_equal(nrow(b$results$cleaned_input$expr_mat), 30L)

    session$setInputs(min_count = 0, min_samples = 0)
    session$flushReact()
    expect_equal(nrow(last_bundle()$results$cleaned_input$expr_mat), 40L)
  })
})

test_that("an emptied threshold box does not break the preview", {
  proj <- shiny::reactiveVal(sm_project())
  shiny::testServer(qc_view_server, args = list(current_project = proj), {
    session$setInputs(layer = "rna")
    session$flushReact()
    # numericInput reports NA while the field is empty, and NA satisfies
    # is.numeric() -- so it reached assert_number() and the preview failed
    # on the keystroke that cleared the box.
    session$setInputs(min_count = NA_real_, min_samples = NA_real_)
    session$flushReact()
    expect_null(last_error())
    expect_false(is.null(last_bundle()))
  })
})

test_that("intensity layers are never filtered by a count threshold", {
  proj <- shiny::reactiveVal(sm_project())
  shiny::testServer(qc_view_server, args = list(current_project = proj), {
    session$setInputs(layer = "rna", min_count = 250, min_samples = 8)
    session$flushReact()
    expect_true(rendered(output$ui_low_expression))

    session$setInputs(layer = "prot")
    session$flushReact()
    expect_null(output$ui_low_expression)
    b <- last_bundle()
    expect_identical(b$params$min_count, 0)
    # nothing dropped for being "low" on a scale where 10000 is meaningless
    expect_equal(nrow(b$results$cleaned_input$expr_mat), 40L)
  })
})

# An RNA-seq layer where the filter does not apply is the case that reads
# as a bug: the user is on RNA-seq, the control they were told about is
# not there, and nothing on the page says which of the two conditions
# failed. Proteomics is different -- nobody expects a count filter there --
# so that stays silent.
test_that("an RNA-seq layer that cannot be filtered says so", {
  logcpm <- sm_input("rnaseq")
  logcpm$assay_type <- "logcpm"
  proj <- shiny::reactiveVal(
    omicsCore::omics_project("smoke", list(prot = sm_input("proteomics"),
                                           rna = logcpm)))
  shiny::testServer(qc_view_server, args = list(current_project = proj), {
    session$setInputs(layer = "rna")
    session$flushReact()
    html <- ui_html(output$ui_low_expression)
    expect_match(html, "raw counts", fixed = TRUE)
    expect_match(html, "logcpm", fixed = TRUE)
    # and it is a note, not a control
    expect_false(grepl("min_count", html, fixed = TRUE))

    # proteomics says nothing: a count filter was never on offer there
    session$setInputs(layer = "prot")
    session$flushReact()
    expect_null(output$ui_low_expression)
  })
})

# The assay type decides whether two of this page's controls appear, and it
# was the one property of the layer the header did not show.
test_that("the header names the assay type", {
  proj <- shiny::reactiveVal(sm_project())
  shiny::testServer(qc_view_server, args = list(current_project = proj), {
    session$setInputs(layer = "rna")
    session$flushReact()
    expect_match(ui_html(output$header), "raw_count", fixed = TRUE)

    session$setInputs(layer = "prot")
    session$flushReact()
    expect_match(ui_html(output$header), "raw_intensity", fixed = TRUE)
  })
})

test_that("only methods whose packages are installed are offered", {
  proj <- shiny::reactiveVal(sm_project())
  shiny::testServer(qc_view_server, args = list(current_project = proj), {
    session$setInputs(layer = "prot")
    session$flushReact()
    # Grouped by MAR/MNAR for the dropdown, so flatten before asking
    # what is on offer.
    choices <- unlist(impute_choices(), use.names = FALSE)
    expect_true("none" %in% choices)
    expect_true("MinProb" %in% choices)
    # An option that errors on selection is worse than one not offered.
    needs <- c(MinProb = "imputeLCMD", QRILC = "imputeLCMD", bpca = "pcaMethods")
    for (m in names(needs)) {
      if (!has_pkg(needs[[m]])) {
        expect_false(m %in% choices, info = m)
      }
    }
  })
})

# ---- two layers of one omics type -------------------------------------
#
# The failure the review reproduced: a project with `skin` and `follicle`
# proteomics, a differential run on follicle. Both call sites resolved the
# primary layer by omics_type and took the first match -- skin -- so the
# "other" layer was follicle, which is where the primary already came
# from. Concordance compared follicle with itself and labelled it
# skin-against-follicle, and the exported script read the wrong file.

two_of_a_kind <- function() {
  mk <- function(seed) {
    set.seed(seed)
    ids <- paste0("F", 1:40); sids <- paste0("s", 1:8)
    m <- matrix(stats::rnorm(320, 20, 2), 40, dimnames = list(ids, sids))
    omicsCore::omics_input(
      m,
      data.frame(condition = rep(c("G1", "G2"), each = 4),
                 donor = paste0("D", 1:8), row.names = sids),
      data.frame(feature_id = ids, feature_symbol = ids, row.names = ids),
      omics_type = "proteomics", assay_type = "normalized_intensity")
  }
  omicsCore::omics_project("two-of-a-kind",
                           list(skin = mk(1), follicle = mk(2)))
}

test_that("integration will not pair a layer with itself", {
  proj <- two_of_a_kind()
  # A differential that says which layer it ran on.
  b <- omicsCore::run_diff(proj$experiments$follicle, method = "ttest",
                           analysis_type = "group", group_col = "condition",
                           control_group = "G1", case_group = "G2")
  b$input_info$layer_tag <- "follicle"

  shiny::testServer(
    integration_view_server,
    args = list(current_project = shiny::reactiveVal(proj),
                diff_bundle = shiny::reactiveVal(b)), {
      session$flushReact()
      info <- can_run()
      expect_true(isTRUE(info$ok))
      expect_identical(info$primary_tag, "follicle")
      expect_identical(info$secondary_tag, "skin")
      expect_false(identical(info$primary_tag, info$secondary_tag))
    })
})

test_that("an untagged result in an ambiguous project is refused, not guessed", {
  proj <- two_of_a_kind()
  b <- omicsCore::run_diff(proj$experiments$follicle, method = "ttest",
                           analysis_type = "group", group_col = "condition",
                           control_group = "G1", case_group = "G2")
  b$input_info$layer_tag <- NULL          # an older project

  shiny::testServer(
    integration_view_server,
    args = list(current_project = shiny::reactiveVal(proj),
                diff_bundle = shiny::reactiveVal(b)), {
      session$flushReact()
      info <- can_run()
      expect_false(isTRUE(info$ok))
      expect_match(info$why, "more than one proteomics layer")
    })
})

# ---- expensive work runs when the work changes ------------------------
#
# can_run() reads the whole project, so every write to it re-evaluated --
# including the write this view's own result causes. One differential
# therefore launched the integration twice: once when the prerequisites
# became true, once when the result of the first run was filed. The second
# computed the same numbers and replaced the first.

test_that("filing the integration's own result does not re-run it", {
  proj <- sm_project()
  pv <- shiny::reactiveVal(proj)
  b <- omicsCore::run_diff(proj$experiments$prot, method = "ttest",
                           analysis_type = "group", group_col = "condition",
                           control_group = "G1", case_group = "G2")
  b$input_info$layer_tag <- "prot"
  b$provenance_binding <- list(analysis_id = "ana_primary")

  runs <- 0L
  testthat::local_mocked_bindings(
    run_async = function(func, on_success, on_error, ...) {
      runs <<- runs + 1L
      invisible(NULL)          # never completes: counting launches only
    },
    .package = "omicsApp")

  shiny::testServer(
    integration_view_server,
    args = list(current_project = pv, diff_bundle = shiny::reactiveVal(b)), {
      session$flushReact()
      expect_equal(runs, 1L)

      # A write to the project that changes nothing this view depends on.
      p2 <- proj
      p2$name <- "renamed"
      p2$bundles <- list(integration = list(dummy = TRUE))
      pv(p2)
      session$flushReact()
      expect_equal(runs, 1L)
    })
})

test_that("a change that does change the work does re-run it", {
  proj <- sm_project()
  pv <- shiny::reactiveVal(proj)
  b <- omicsCore::run_diff(proj$experiments$prot, method = "ttest",
                           analysis_type = "group", group_col = "condition",
                           control_group = "G1", case_group = "G2")
  b$input_info$layer_tag <- "prot"
  b$provenance_binding <- list(analysis_id = "ana_primary")
  bv <- shiny::reactiveVal(b)

  runs <- 0L
  testthat::local_mocked_bindings(
    run_async = function(func, on_success, on_error, ...) {
      runs <<- runs + 1L
      invisible(NULL)
    },
    .package = "omicsApp")

  shiny::testServer(
    integration_view_server,
    args = list(current_project = pv, diff_bundle = bv), {
      session$flushReact()
      expect_equal(runs, 1L)

      # A different differential upstream is different work.
      b2 <- b
      b2$provenance_binding <- list(analysis_id = "ana_other")
      bv(b2)
      session$flushReact()
      expect_equal(runs, 2L)
    })
})

# This observer reads active(), which reads the project, so every write to
# the project re-ran run_qc -- a rename, a differential being filed, an
# enrichment result arriving. On the demo that is milliseconds; on
# 63k x 129 it is not, and none of those writes changes what QC would say.
test_that("an unrelated write to the project does not re-run QC", {
  proj <- sm_project()
  pv <- shiny::reactiveVal(proj)
  calls <- 0L
  real_run_qc <- omicsCore::run_qc
  testthat::local_mocked_bindings(
    run_qc = function(...) {
      calls <<- calls + 1L
      real_run_qc(...)
    },
    .package = "omicsCore")

  shiny::testServer(qc_view_server, args = list(current_project = pv), {
    session$flushReact()
    first <- calls
    expect_gt(first, 0L)

    p2 <- proj
    p2$name <- "renamed"
    p2$bundles <- list(enrich = list(dummy = TRUE))
    pv(p2)
    session$flushReact()
    expect_equal(calls, first)

    # A control the user moved is different work, and does re-run.
    session$setInputs(missing_threshold = 0.01)
    session$flushReact()
    expect_gt(calls, first)
  })
})
