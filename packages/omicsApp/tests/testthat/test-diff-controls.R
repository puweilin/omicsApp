# Four reports about the Differential view, three of which were real:
# DESeq2 unreachable, no raw-p option, and a fold-change slider that
# could not express the fold change people actually use.

test_that("the layer control is what makes DESeq2 reachable", {
  # The gate itself was always right -- deseq2 and edger need rnaseq
  # raw counts. The view pinned itself to the first proteomics layer,
  # so there was no way to the side of the gate where it opens.
  shiny::testServer(diff_view_server, args = list(), {
    session$flushReact()
    expect_identical(active()$tag, "proteomics")
    expect_false("deseq2" %in%
                   omicsCore::applicable_diff_methods(active()$input))

    session$setInputs(layer = "rnaseq")
    expect_identical(active()$tag, "rnaseq")
    expect_true("deseq2" %in%
                  omicsCore::applicable_diff_methods(active()$input))
  })
})

test_that("a layer name that is not in the project is ignored", {
  shiny::testServer(diff_view_server, args = list(), {
    session$setInputs(layer = "no_such_layer")
    expect_identical(active()$tag, "proteomics")
    expect_false(is.null(active()$input))
  })
})

# ---- thresholds -------------------------------------------------------

diff_frame <- function() {
  data.frame(
    feature_id  = paste0("f", 1:4),
    effect      = c(2, 2, -2, 0.1),
    p_value     = c(0.001, 0.02, 0.001, 0.001),
    adj_p_value = c(0.01,  0.30, 0.01,  0.01)
  )
}

test_that("significance follows the chosen p column", {
  # f2 is the discriminating row: raw p 0.02 passes, adj.P 0.30 does
  # not. An exploratory screen and a confirmatory one want different
  # answers here, and only adj.P was on offer.
  shiny::testServer(diff_view_server, args = list(), {
    session$setInputs(fdr_cut = 0.05, fc_cut = 1, p_kind = "adj")
    diff_bundle(list(results = list(diff_result_df = diff_frame())))
    session$elapse(300)
    expect_identical(marked()$is_significant, c(TRUE, FALSE, TRUE, FALSE))

    session$setInputs(p_kind = "raw")
    session$elapse(300)
    expect_identical(marked()$is_significant, c(TRUE, TRUE, TRUE, FALSE))
  })
})

test_that("the label follows the column, so no figure says adj.P over raw p", {
  shiny::testServer(diff_view_server, args = list(), {
    session$setInputs(p_kind = "adj")
    expect_identical(p_label(), "adj.P")
    expect_identical(p_col(), "adj_p_value")

    session$setInputs(p_kind = "raw")
    expect_identical(p_label(), "p")
    expect_identical(p_col(), "p_value")
  })
})

test_that("the fold-change cutoff defaults to log2(1.2)", {
  # The old slider stepped 0.05 from 0, so 0.263 was one of the values
  # it could not reach -- and it is the one most often wanted.
  shiny::testServer(diff_view_server, args = list(), {
    session$elapse(300)
    expect_equal(fc_cut_d(), round(log2(1.2), 3))
  })
})

test_that("an empty or negative fold-change box falls back, not to nothing", {
  # numericInput reports NA while the box is mid-edit. NA in the
  # comparison marks every feature non-significant, with nothing on
  # screen to say why.
  shiny::testServer(diff_view_server, args = list(), {
    session$setInputs(fc_cut = NA_real_)
    session$elapse(300)
    expect_equal(fc_cut_d(), round(log2(1.2), 3))

    session$setInputs(fc_cut = -1)
    session$elapse(300)
    expect_equal(fc_cut_d(), round(log2(1.2), 3))

    session$setInputs(fc_cut = 0.5)
    session$elapse(300)
    expect_equal(fc_cut_d(), 0.5)
  })
})

test_that("a fold-change the old slider could not express now works", {
  shiny::testServer(diff_view_server, args = list(), {
    session$setInputs(fdr_cut = 0.05, p_kind = "adj",
                      fc_cut = round(log2(1.2), 3))
    df <- diff_frame()
    df$effect <- c(0.3, 0.2, -0.3, 0.1)   # straddles 0.263
    diff_bundle(list(results = list(diff_result_df = df)))
    session$elapse(300)
    expect_identical(marked()$is_significant, c(TRUE, FALSE, TRUE, FALSE))
  })
})

# ---- first paint ------------------------------------------------------

test_that("opening the view does not produce an error", {
  # do_run() fires as soon as active() settles, which is before the
  # renderUI-driven contrast controls exist. Reading input$group_col
  # there gets NULL, and refusing on NULL put "Pick a group column with
  # distinct Control and Case levels" on a view the user had just
  # opened, having touched nothing.
  shiny::testServer(diff_view_server, args = list(), {
    session$flushReact()
    session$setInputs(layer = "proteomics")   # what makes active() settle
    expect_null(diff_error())
    expect_false(is.null(diff_bundle()))
  })
})

test_that("the first run uses the same contrast the controls will show", {
  shiny::testServer(diff_view_server, args = list(), {
    session$flushReact()
    session$setInputs(layer = "proteomics")
    d <- default_contrast()
    b <- diff_bundle()
    expect_identical(b$params$control_group, d$control)
    expect_identical(b$params$case_group, d$case)
  })
})

test_that("switching layer re-runs with an engine that layer supports", {
  shiny::testServer(diff_view_server, args = list(), {
    session$flushReact()
    session$setInputs(layer = "rnaseq")
    expect_null(diff_error())
    # auto resolves to deseq2 for raw counts; the point is that it ran
    # at all, on a layer that was unreachable before.
    expect_identical(diff_bundle()$params$method, "deseq2")
  })
})

test_that("a project with no experiments says so rather than blaming the contrast", {
  empty <- omicsCore::omics_project(name = "empty", experiments = list())
  shiny::testServer(
    diff_view_server,
    args = list(current_project = shiny::reactiveVal(empty)),
    {
      session$setInputs(rerun = 1)
      expect_match(diff_error() %||% "", "no experiments", ignore.case = TRUE)
    }
  )
})
