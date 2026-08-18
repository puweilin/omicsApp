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

test_that("diff view falls back to example_diff_bundle when project is NULL", {
  current_project <- shiny::reactiveVal(NULL)
  shiny::testServer(
    diff_view_server,
    args = list(current_project = current_project),
    {
      # Demo bundle no longer auto-populates on init; Re-run triggers it.
      session$setInputs(rerun = 1)
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
  parsed <- omicsCore::read_omics(xlsx, omics_type = "proteomics")
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
  parsed <- omicsCore::read_omics(xlsx, omics_type = "proteomics")
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
