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

# ---- changing layer ---------------------------------------------------

test_that("changing layer clears the result rather than re-running", {
  # Picking a layer is the first of several decisions -- method,
  # contrast, covariates. Spending a run on the state halfway through
  # them is work nobody asked for, on settings nobody has finished
  # choosing.
  shiny::testServer(diff_view_server, args = list(), {
    session$flushReact()
    session$setInputs(layer = "proteomics")
    expect_false(is.null(diff_bundle()))          # the one automatic run

    session$setInputs(layer = "rnaseq")
    expect_null(diff_bundle())
    expect_null(diff_error())
  })
})

test_that("keeping the old result would have been worse than clearing", {
  # A bundle from the previous layer describes different features
  # entirely, and nothing on screen would say which layer it came from.
  shiny::testServer(diff_view_server, args = list(), {
    session$flushReact()
    session$setInputs(layer = "proteomics")
    before <- diff_bundle()$input_info$omics_type

    session$setInputs(layer = "rnaseq")
    session$setInputs(rerun = 1)
    expect_identical(diff_bundle()$input_info$omics_type, "rnaseq")
    expect_false(identical(diff_bundle()$input_info$omics_type, before))
  })
})

test_that("Re-run after a layer change uses the new layer's engine", {
  shiny::testServer(diff_view_server, args = list(), {
    session$flushReact()
    session$setInputs(layer = "proteomics")
    expect_identical(diff_bundle()$params$method, "limma")

    session$setInputs(layer = "rnaseq", rerun = 1)
    expect_identical(diff_bundle()$params$method, "deseq2")
  })
})

# ---- the controls must not undo the user ------------------------------

real_shaped_project <- function(n = 40L, nf = 300L) {
  set.seed(7L)
  ids <- sprintf("RD%03d-C", seq_len(n))
  m <- matrix(rnorm(nf * n, 20, 1.5), nf, n,
              dimnames = list(sprintf("A0A%03d", seq_len(nf)), ids))
  cond <- rep(c("G1", "G2"), each = n / 2)
  m[1:30, cond == "G2"] <- m[1:30, cond == "G2"] + 3
  # `label` holds the sample names, as a real workbook's does.
  meta <- data.frame(label = ids, tissue = rep("Cheek", n),
                     condition = cond, row.names = ids,
                     stringsAsFactors = FALSE)
  omicsCore::omics_project(
    name = "real",
    experiments = list(proteomics = omicsCore::omics_input(
      m, meta, data.frame(feature_id = rownames(m)),
      omics_type = "proteomics", assay_type = "normalized_intensity")))
}

test_that("a column of sample identifiers is not offered as a contrast", {
  # `label` satisfied the old rule -- not numeric, more than one value --
  # so the view defaulted to one sample against one other. limma cannot
  # fit that, and said so as "Partial NA coefficients for 2294 probe(s)"
  # and an empty result.
  meta <- real_shaped_project()$experiments$proteomics$meta_df
  cands <- grouping_candidates(meta)
  expect_false("label" %in% cands)
  expect_identical(cands[[1L]], "condition")
})

test_that("every offered column has at least two samples per level", {
  meta <- real_shaped_project()$experiments$proteomics$meta_df
  for (nm in grouping_candidates(meta)) {
    expect_gte(min(table(meta[[nm]])), 2L, label = nm)
  }
})

test_that("a named column wins over one that merely has fewer levels", {
  meta <- data.frame(batch = rep(c("a", "b"), each = 4),
                     condition = rep(c("ctl", "trt", "x", "y"), each = 2),
                     stringsAsFactors = FALSE)
  # batch has fewer levels; condition says what it is.
  expect_identical(grouping_candidates(meta)[[1L]], "condition")
})

test_that("a finished result survives the controls re-rendering", {
  # ui_layer read active(), and active() reads input$layer, so
  # re-rendering the control re-sent its value and invalidated active()
  # -- which cleared the result. The user picked condition / G1 / G2,
  # ran it, and the view came back saying "no result yet" with `label`
  # selected.
  shiny::testServer(
    diff_view_server,
    args = list(current_project = shiny::reactiveVal(real_shaped_project())),
    {
      session$flushReact()
      session$setInputs(layer = "proteomics")
      session$setInputs(group_col = "condition", control = "G1",
                        case = "G2", rerun = 1)
      expect_false(is.null(diff_bundle()))

      invisible(output$ui_group_col)
      invisible(output$ui_layer)
      invisible(output$ui_contrast)
      session$flushReact()

      expect_false(is.null(diff_bundle()))
      expect_identical(input$group_col, "condition")
      expect_identical(input$control, "G1")
      expect_identical(input$case, "G2")
    }
  )
})
