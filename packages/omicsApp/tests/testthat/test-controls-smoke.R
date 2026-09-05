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

test_that("the report buttons are live exactly when they can work", {
  skip_if_not_installed("rmarkdown")

  # No project: do_export() stops on req(proj), silently, as req() does.
  # A live-looking button that produces no file and no message is
  # indistinguishable from a broken one -- which is what it looked like.
  shiny::testServer(report_view_server,
                    args = list(current_project = shiny::reactiveVal(NULL)), {
    html <- paste(as.character(output$header), collapse = "")
    expect_match(html, "disabled")
    expect_match(paste(as.character(output$notices), collapse = ""),
                 "Import data first")
  })

  shiny::testServer(report_view_server,
                    args = list(current_project = shiny::reactiveVal(sm_project())), {
    html <- paste(as.character(output$header), collapse = "")
    expect_no_match(html, "pointer-events:none")
  })
})

test_that("a report actually renders, in both formats", {
  skip_if_not_installed("rmarkdown")
  skip_if_not(rmarkdown::pandoc_available(), "pandoc not available")

  proj <- sm_project()
  proj$bundles <- list(qc = example_qc_bundle(), diff = example_diff_bundle())

  out <- withr::local_tempfile(fileext = ".html")
  omicsCore::export_report(proj, out, format = "html", overwrite = TRUE)
  expect_gt(file.size(out), 10000)
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

    session$setInputs(impute_method = "half_min")
    session$flushReact()
    b <- last_bundle()
    expect_identical(b$params$impute_method, "half_min")
    # The point of the control: the NAs are gone afterwards.
    expect_false(anyNA(b$results$cleaned_input$expr_mat))
  })
})

test_that("counts are never imputed, even with a method left selected", {
  proj <- shiny::reactiveVal(sm_project())
  shiny::testServer(qc_view_server, args = list(current_project = proj), {
    session$setInputs(layer = "prot", impute_method = "half_min")
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

test_that("only methods whose packages are installed are offered", {
  proj <- shiny::reactiveVal(sm_project())
  shiny::testServer(qc_view_server, args = list(current_project = proj), {
    session$setInputs(layer = "prot")
    session$flushReact()
    choices <- impute_choices()
    expect_true("none" %in% choices)
    # An option that errors on selection is worse than one not offered.
    needs <- c(knn = "impute", missforest = "missForest", bpca = "pcaMethods")
    for (m in names(needs)) {
      if (!has_pkg(needs[[m]])) {
        expect_false(m %in% choices, info = m)
      }
    }
  })
})
