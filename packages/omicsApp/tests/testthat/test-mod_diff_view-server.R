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

# ---- the volcano follows the threshold controls -----------------------

# Where ggplotly() put the volcano's dashed reference lines (line traces).
volcano_line_x <- function(fig) {
  lines <- Filter(function(t) identical(t$mode, "lines"), fig$data)
  unlist(lapply(lines, function(t) unlist(t$x)))
}

test_that("the volcano is drawn at the thresholds the controls set", {
  shiny::testServer(
    diff_view_server,
    args = list(current_project = shiny::reactiveVal(NULL)),
    {
      session$setInputs(group_col = "group", control = "G1", case = "G2",
                        rerun = 1, fdr_cut = 0.05, fc_cut = 0.5)
      session$elapse(300)
      fig <- jsonlite::fromJSON(output$volcano, simplifyVector = FALSE)$x
      # The fold-change cut is on the figure, on both sides of zero...
      expect_true(all(c(-0.5, 0.5) %in% volcano_line_x(fig)))
      # ...and the points it colours are the hits the table lists, up
      # and down apart.
      trace_names <- vapply(fig$data, function(t) t$name %||% "", "")
      n_up   <- sum(marked()$is_significant & marked()$effect > 0)
      n_down <- sum(marked()$is_significant & marked()$effect < 0)
      expect_gt(n_up + n_down, 0L)
      if (n_up > 0L) expect_true(sprintf("up (%d)", n_up) %in% trace_names)
      if (n_down > 0L) expect_true(sprintf("down (%d)", n_down) %in% trace_names)
      # The card's own legend is the one shown; plotly's would repeat it.
      expect_false(isTRUE(fig$layout$showlegend))
      # The direction-coloured figure still converts cleanly, its traces
      # named in the legend's plain words.
      p <- omicsCore::plot_volcano(shown_bundle(), top_n = 0L, effect_threshold = 0.5)
      built <- expect_no_warning(plotly::plotly_build(plotly::ggplotly(p, tooltip = "text")))
      built_names <- vapply(built$x$data, function(t) t$name %||% "", "")
      expect_true("not significant" %in% built_names)
      expect_false(any(c("ns", "significant", "TRUE", "FALSE") %in% built_names))
      # And through WebGL as the card draws it: an invisible helper layer
      # once became a "gl" trace plotly warned about on every build.
      expect_no_warning(plotly::plotly_build(plotly::toWebGL(
        drop_hoveron(plotly::ggplotly(p, tooltip = "text")))))
      # The same when nothing passes and the hit layer is empty.
      none <- omicsCore::plot_volcano(shown_bundle(), top_n = 0L, p_threshold = 1e-300)
      expect_no_warning(plotly::plotly_build(plotly::toWebGL(
        drop_hoveron(plotly::ggplotly(none, tooltip = "text")))))

      session$setInputs(fc_cut = 1)
      session$elapse(300)
      fig <- jsonlite::fromJSON(output$volcano, simplifyVector = FALSE)$x
      expect_true(all(c(-1, 1) %in% volcano_line_x(fig)))
      expect_false(0.5 %in% volcano_line_x(fig))

      # No fold-change cut, no fold-change lines.
      session$setInputs(fc_cut = 0)
      session$elapse(300)
      fig <- jsonlite::fromJSON(output$volcano, simplifyVector = FALSE)$x
      expect_false(any(c(-1, 1, -0.5, 0.5) %in% volcano_line_x(fig)))
    }
  )
})

test_that("the volcano card states the cut, and follows the controls", {
  shiny::testServer(
    diff_view_server,
    args = list(current_project = shiny::reactiveVal(NULL)),
    {
      session$setInputs(group_col = "group", control = "G1", case = "G2",
                        rerun = 1, fdr_cut = 0.05, fc_cut = 0.263)
      session$elapse(300)
      # plotly shows no caption, so the card says what "significant" is.
      expect_identical(output$volcano_cut,
                       "significant = adjusted p < 0.05 and |log2FC| \u2265 0.263")
      session$setInputs(fdr_cut = 0.01, fc_cut = 1)
      session$elapse(300)
      expect_identical(output$volcano_cut,
                       "significant = adjusted p < 0.01 and |log2FC| \u2265 1")
      session$setInputs(p_kind = "raw")
      session$elapse(300)
      expect_match(output$volcano_cut, "significant = p < 0.01", fixed = TRUE)
      # A cutoff of 0 is no cutoff, and the card does not claim one.
      session$setInputs(fc_cut = 0)
      session$elapse(300)
      expect_identical(output$volcano_cut, "significant = p < 0.01")
      # The legend under the plot counts the same hits as the stat cards.
      legend <- output$volcano_legend$html
      n_up   <- sum(marked()$is_significant & marked()$effect > 0)
      n_down <- sum(marked()$is_significant & marked()$effect < 0)
      expect_match(legend, sprintf("up (%d)", n_up), fixed = TRUE)
      expect_match(legend, sprintf("down (%d)", n_down), fixed = TRUE)
      expect_match(legend, "not significant", fixed = TRUE)
    }
  )
})

test_that("the volcano's labels sit at the p-value it is drawn on", {
  shiny::testServer(
    diff_view_server,
    args = list(current_project = shiny::reactiveVal(NULL)),
    {
      session$setInputs(group_col = "group", control = "G1", case = "G2",
                        rerun = 1, p_kind = "raw", label_top = TRUE)
      session$elapse(300)
      fig <- jsonlite::fromJSON(output$volcano, simplifyVector = FALSE)$x
      df <- shown_bundle()$results$diff_result_df
      top <- df[order(df$p_value), ][1, ]
      ann <- fig$layout$annotations[[1]]
      expect_identical(ann$text, top$feature_symbol)
      # -log10 of the raw p, the y axis of a raw-p volcano: placed at the
      # adjusted p, every label floated below its point.
      expect_equal(ann$y, -log10(top$p_value), tolerance = 1e-6)
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

test_that("the volcano card has a line for the cut and an up/down legend", {
  html <- render_html(diff_volcano_card(function(x) x))
  expect_match(html, "volcano_cut", fixed = TRUE)
  expect_match(html, "volcano_legend", fixed = TRUE)
  expect_no_match(html, "fixed thresholds", fixed = TRUE)
  legend <- render_html(volcano_legend(412, 388))
  expect_match(legend, "up (412)", fixed = TRUE)
  expect_match(legend, "down (388)", fixed = TRUE)
  expect_match(legend, "not significant", fixed = TRUE)
  expect_match(legend, omics_colors$up, fixed = TRUE)
  expect_match(legend, omics_colors$down, fixed = TRUE)
  # Along a continuous variable the effect is a slope: positive or negative.
  expect_match(render_html(volcano_legend(1, 2, continuous = TRUE)), "negative (2)",
               fixed = TRUE)
  # Before a run there is nothing to count, but the key is still there.
  expect_match(render_html(volcano_legend()), "</span>\\s*up\\s*</span>")
})
