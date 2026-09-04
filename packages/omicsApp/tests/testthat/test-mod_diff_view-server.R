# testServer harness for the Differential view (slice 3D).
#
# Three scenarios:
#   1. NULL project → demo fallback bundle from
#      `example_diff_bundle()` populates the volcano + hits.
#   2. Live project from the tiny synthetic xlsx → Re-run with
#      G2 vs G1 contrast emits a real `analysis_bundle` with
#      `diff_result_df`.
#   3. Mis-configured contrast (Control == Case) → error notice
#      surfaces in `diff_error()` instead of stop()-ing the
#      reactive graph.

test_that("with no project the view analyses the demo, for real", {
  # It used to hand back example_diff_bundle(), a bundle computed
  # elsewhere at fixed settings, so the Method dropdown it had just
  # drawn had no bearing on the result. The demo now goes through
  # run_diff() like a project does, which is also what makes the
  # rnaseq layer -- and so deseq2 -- reachable.
  current_project <- shiny::reactiveVal(NULL)
  shiny::testServer(
    diff_view_server,
    args = list(current_project = current_project),
    {
      session$setInputs(group_col = "group", control = "G1", case = "G2",
                        method = "auto", rerun = 1)
      b <- diff_bundle()
      expect_s3_class(b, "analysis_bundle")
      expect_identical(b$analysis_name, "run_diff")
      expect_true(nrow(b$results$diff_result_df) > 0L)
      expect_null(diff_error())
    }
  )
})

test_that("diff view runs run_diff() against the live project", {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("readxl")

  xlsx <- tempfile(fileext = ".xlsx")
  on.exit(unlink(xlsx), add = TRUE)
  write_tiny_omics_xlsx(xlsx)
  parsed <- omicsCore::read_omics(xlsx, omics_type = "proteomics", assay_type = "normalized_intensity")
  inp <- parsed$input
  proj <- omicsCore::omics_project(
    name        = "test",
    experiments = list(proteomics = inp)
  )
  current_project <- shiny::reactiveVal(proj)

  shiny::testServer(
    diff_view_server,
    args = list(current_project = current_project),
    {
      # Drive the contrast inputs the same way the renderUI shells
      # would. ttest is always available; pin it to bypass any
      # missing Bioconductor backends in the test env.
      session$setInputs(
        method     = "ttest",
        group_col  = "group",
        control    = "G1",
        case       = "G2",
        covariates = character(0)
      )
      # observeEvent(active()) auto-runs once on init using "auto"
      # method. Click Re-run to apply the ttest pin.
      session$setInputs(rerun = 1)
      b <- diff_bundle()
      expect_s3_class(b, "analysis_bundle")
      expect_equal(b$params$method, "ttest")
      expect_equal(b$params$control_group, "G1")
      expect_equal(b$params$case_group, "G2")
      df <- b$results$diff_result_df
      expect_true(nrow(df) > 0L)
      expect_true(all(c("feature_id", "effect", "adj_p_value") %in% names(df)))
      expect_null(diff_error())
    }
  )
})

test_that("diff view surfaces validation errors instead of crashing", {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("readxl")

  xlsx <- tempfile(fileext = ".xlsx")
  on.exit(unlink(xlsx), add = TRUE)
  write_tiny_omics_xlsx(xlsx)
  parsed <- omicsCore::read_omics(xlsx, omics_type = "proteomics", assay_type = "normalized_intensity")
  proj <- omicsCore::omics_project(
    name        = "test",
    experiments = list(proteomics = parsed$input)
  )
  current_project <- shiny::reactiveVal(proj)

  shiny::testServer(
    diff_view_server,
    args = list(current_project = current_project),
    {
      session$setInputs(
        method    = "ttest",
        group_col = "group",
        control   = "G1",
        case      = "G1",
        rerun     = 1
      )
      expect_false(is.null(diff_error()))
      expect_match(diff_error(), "Control and Case|distinct", fixed = FALSE)
    }
  )
})

# ---- method gating ----------------------------------------------------

test_that("the demo (proteomics) view hides the count engines", {
  shiny::testServer(
    diff_view_server,
    args = list(current_project = shiny::reactiveVal(NULL)),
    {
      html <- render_html(output$ui_method)
      expect_match(html, "limma", fixed = TRUE)
      # Offering DESeq2 for intensities would let a user produce a full,
      # plausible, meaningless result table with no warning.
      expect_false(grepl("deseq2", html, fixed = TRUE))
      expect_false(grepl("edger", html, fixed = TRUE))
      expect_match(render_html(output$method_note), "hidden", fixed = TRUE)
    }
  )
})

test_that("a raw-count layer offers the count engines instead", {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("readxl")
  xlsx <- tempfile(fileext = ".xlsx"); on.exit(unlink(xlsx), add = TRUE)
  write_tiny_omics_xlsx(xlsx)
  inp <- omicsCore::read_omics(xlsx, omics_type = "rnaseq",
                               assay_type = "raw_count")$input
  proj <- omicsCore::omics_project("counts",
                                   experiments = list(rnaseq = inp))
  shiny::testServer(
    diff_view_server,
    args = list(current_project = shiny::reactiveVal(proj)),
    {
      html <- render_html(output$ui_method)
      expect_match(html, "deseq2", fixed = TRUE)
      expect_match(html, "edger", fixed = TRUE)
      # limma here has no voom step, so counts are not its business.
      expect_false(grepl(">limma<", html, fixed = TRUE))
    }
  )
})

# ---- the volcano is not driven by the sliders -------------------------

test_that("the volcano does not depend on the threshold sliders", {
  shiny::testServer(
    diff_view_server,
    args = list(current_project = shiny::reactiveVal(NULL)),
    {
      session$setInputs(group_col = "group", control = "G1", case = "G2",
                        rerun = 1)
      shiny::req(diff_bundle())
      caption_of <- function() {
        p <- omicsCore::plot_volcano(
          diff_bundle(),
          top_n = if (isTRUE(input$label_top)) 20L else 0L)
        ggplot2::ggplot_build(p)$plot$labels$caption
      }
      before <- caption_of()
      session$setInputs(fdr_cut = 0.5, fc_cut = 0.1)
      # The figure is the stable reference the hit table is read against;
      # a screenshot of it must not depend on where a control was left.
      expect_identical(caption_of(), before)
      expect_match(before, "adj_p_value < 0.05", fixed = TRUE)
    }
  )
})

test_that("the sliders still drive the hit table", {
  shiny::testServer(
    diff_view_server,
    args = list(current_project = shiny::reactiveVal(NULL)),
    {
      session$setInputs(group_col = "group", control = "G1", case = "G2",
                        rerun = 1, fdr_cut = 0.05, fc_cut = 1)
      # The thresholds are debounced, so the mask only follows them once
      # the window has passed. A real drag settles the same way.
      session$elapse(300)
      strict <- sum(marked()$is_significant)
      session$setInputs(fdr_cut = 1, fc_cut = 0)
      session$elapse(300)
      # Sweeping a threshold is the useful thing to do to a table, and
      # that is where it still happens.
      expect_gt(sum(marked()$is_significant), strict)
    }
  )
})

test_that("the volcano card says the sliders do not reach it", {
  html <- render_html(diff_volcano_card(function(x) x))
  expect_match(html, "sliders filter the table", fixed = TRUE)
})
