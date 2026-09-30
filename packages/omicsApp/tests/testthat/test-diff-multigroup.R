# One control, several treatments, in the Differential view.

mg_project <- function(levels = c("DMSO", "DrugA", "DrugB", "DrugC")) {
  set.seed(5)
  groups <- rep(levels, each = 3)
  samp <- paste0("S", seq_along(groups))
  ids <- paste0("P", 1:40)
  m <- matrix(stats::rnorm(40 * length(groups), 12, 0.3), 40,
              dimnames = list(ids, samp))
  m[1:6, groups == "DrugA"] <- m[1:6, groups == "DrugA"] + 3
  m[7:12, groups == "DrugC"] <- m[7:12, groups == "DrugC"] - 3
  inp <- omicsCore::omics_input(
    m, data.frame(treatment = groups, row.names = samp, stringsAsFactors = FALSE),
    data.frame(feature_id = ids, feature_symbol = paste0("G", 1:40), row.names = ids,
               stringsAsFactors = FALSE),
    omics_type = "proteomics", assay_type = "normalized_intensity")
  omicsCore::omics_project("mg", list(proteomics = inp))
}

test_that("the control-looking level is the default reference, and every other group is compared with it", {
  expect_identical(default_control_level(c("DrugA", "DMSO", "DrugB")), "DMSO")
  expect_identical(default_control_level(c("KO", "WT")), "WT")
  expect_identical(default_control_level(c("treated", "Vehicle")), "Vehicle")
  expect_identical(default_control_level(c("G1", "G2")), "G1")
  expect_identical(default_control_level(c("A", "ctrl_24h")), "ctrl_24h")

  shiny::testServer(diff_view_server,
                    args = list(current_project = shiny::reactiveVal(mg_project())), {
    session$flushReact()
    d <- default_contrast()
    expect_identical(d$control, "DMSO")
    expect_identical(d$case, c("DrugA", "DrugB", "DrugC"))
  })
})

test_that("several groups run in one model and the view switches between them", {
  shiny::testServer(diff_view_server,
                    args = list(current_project = shiny::reactiveVal(mg_project())), {
    session$setInputs(group_col = "treatment", control = "DMSO",
                      case = c("DrugA", "DrugB", "DrugC"), method = "limma",
                      rerun = 1)
    b <- diff_bundle()
    expect_identical(omicsCore::diff_comparisons(b),
                     c("DrugA_vs_DMSO", "DrugB_vs_DMSO", "DrugC_vs_DMSO"))
    # What the rest of the app receives is one contrast.
    shown <- session$returned$bundle()
    expect_identical(shown$params$comparison, "DrugA_vs_DMSO")
    expect_true(all(shown$results$diff_result_df$comparison == "DrugA_vs_DMSO"))

    session$setInputs(comparison = "DrugC_vs_DMSO")
    shown <- session$returned$bundle()
    expect_identical(shown$params$case_group, "DrugC")
    expect_identical(shown$params$all_case_groups, c("DrugA", "DrugB", "DrugC"))

    # The summary card appears only with several comparisons.
    expect_false(is.null(output$contrast_summary))
    s <- contrast_summary_df()
    expect_identical(nrow(s), 3L)
    expect_gte(s$n_down[s$comparison == "DrugC_vs_DMSO"], 5L)
  })
})

test_that("the hit table names the p-value it was masked on", {
  shiny::testServer(diff_view_server,
                    args = list(current_project = shiny::reactiveVal(mg_project())), {
    session$setInputs(group_col = "treatment", control = "DMSO", case = "DrugA",
                      method = "limma", p_kind = "raw", rerun = 1)
    expect_identical(p_label(), "p")
    session$setInputs(p_kind = "adj")
    expect_identical(p_label(), "adj.P")
  })
})

test_that("the control cannot also be a case", {
  shiny::testServer(diff_view_server,
                    args = list(current_project = shiny::reactiveVal(mg_project())), {
    session$setInputs(group_col = "treatment", control = "DMSO", case = "DMSO",
                      method = "limma", rerun = 1)
    expect_null(diff_bundle())
    expect_match(diff_error(), "distinct from the control")
  })
})
