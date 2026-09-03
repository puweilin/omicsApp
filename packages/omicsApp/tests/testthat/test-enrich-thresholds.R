# Enrichment took only the diff bundle and applied its own defaults --
# adjusted p at 0.05, no fold-change bound -- so the genes it enriched
# were a different set from the hits on screen, and neither view said
# so. Reported as "100+ differential genes, every database comes back
# empty", which is what a raw-p hit list looks like to an adjusted-p
# selection.

split_bundle <- function(n = 3000L, n_raw = 120L) {
  set.seed(9L)
  df <- data.frame(
    feature_id     = sprintf("P%04d", seq_len(n)),
    feature_symbol = sprintf("SYM%04d", seq_len(n)),
    effect         = rnorm(n, 0, 0.5),
    # Significant by raw p, not one of them by adjusted.
    p_value     = c(runif(n_raw, 0, 0.04), runif(n - n_raw, 0.05, 1)),
    adj_p_value = c(runif(n_raw, 0.2, 0.9), runif(n - n_raw, 0.5, 1)),
    stringsAsFactors = FALSE
  )
  omicsCore::new_analysis_bundle(
    "run_diff", input_info = list(omics_type = "proteomics"),
    params = list(method = "limma", comparison = "G2_vs_G1"),
    results = list(diff_result_df = df))
}

# testServer() takes its block unevaluated, so this builds the args
# rather than wrapping the call: a helper that forwards the block turns
# it into a promise and the module's reactives go out of scope.
enrich_args <- function(bundle, ...) {
  thr <- utils::modifyList(
    list(p_cutoff = 0.05, p_preference = "adjusted", effect_cutoff = NULL),
    list(...))
  list(diff_bundle = shiny::reactiveVal(bundle),
       diff_thresholds = shiny::reactive(thr))
}

test_that("the p column follows the Differential view", {
  b <- split_bundle()
  shiny::testServer(enrich_view_server, args = enrich_args(b, p_preference = "raw"),
                    expect_equal(selected_features()$n, 120L))
  shiny::testServer(enrich_view_server, args = enrich_args(b, p_preference = "adjusted"),
                    expect_equal(selected_features()$n, 0L))
})

test_that("the fold-change bound is applied too", {
  # run_enrichment's own default is no bound at all, so a hit list the
  # user narrowed with the box was widened again on the way here.
  b <- split_bundle()
  shiny::testServer(enrich_view_server,
                    args = enrich_args(b, p_preference = "raw", effect_cutoff = 0),
                    expect_equal(selected_features()$n, 120L))
  shiny::testServer(enrich_view_server,
                    args = enrich_args(b, p_preference = "raw", effect_cutoff = 5),
                    expect_equal(selected_features()$n, 0L))
})

# ---- saying why it is empty -------------------------------------------

test_that("the panel states how many features went in", {
  html <- NULL
  shiny::testServer(enrich_view_server,
                    args = enrich_args(split_bundle(), p_preference = "raw"), {
    session$flushReact()
    html <<- paste(unlist(output$input_summary), collapse = " ")
  })
  expect_true(grepl("120 of 3000 features", html, fixed = TRUE))
  # `<` renders escaped, so assert on the parts that survive it: which
  # column was used and at what cutoff.
  expect_true(grepl("p_value", html, fixed = TRUE))
  expect_false(grepl("adj_p_value", html, fixed = TRUE))
  expect_true(grepl("0.05", html, fixed = TRUE))
})

test_that("an empty selection says so, and says where to change it", {
  # Otherwise "no pathways" reads the same as a real negative result.
  html <- NULL
  shiny::testServer(enrich_view_server,
                    args = enrich_args(split_bundle(), p_preference = "adjusted"), {
    session$flushReact()
    html <<- paste(unlist(output$input_summary), collapse = " ")
  })
  expect_true(grepl("0 of 3000 features", html, fixed = TRUE))
  expect_true(grepl("nothing to enrich", html, fixed = TRUE))
  expect_true(grepl("Differential", html, fixed = TRUE))
})

test_that("features carrying no gene symbol are called out", {
  # feature_symbol falls back to feature_id when the workbook has no
  # gene column; clusterProfiler then maps nothing and returns NULL,
  # which arrives as an empty result rather than a complaint.
  b <- split_bundle()
  b$results$diff_result_df$feature_symbol <- b$results$diff_result_df$feature_id
  html <- NULL
  shiny::testServer(enrich_view_server,
                    args = enrich_args(b, p_preference = "raw"), {
    session$flushReact()
    html <<- paste(unlist(output$input_summary), collapse = " ")
  })
  expect_true(grepl("none carry a gene symbol", html, fixed = TRUE))
})

test_that("nothing is claimed before a differential has been run", {
  shiny::testServer(enrich_view_server, args = list(), {
    expect_null(selected_features())
    expect_null(output$input_summary)
  })
})
