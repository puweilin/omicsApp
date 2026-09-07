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

# ---- GSEA is a ranked-list method -------------------------------------

test_that("the ranked list keeps features no threshold would have kept", {
  # This is the guarantee, not an implementation detail: GSEA reads the
  # whole ordered list, and pre-filtering it is not a stricter analysis
  # but a different and invalid one.
  b <- omicsApp:::example_diff_bundle()
  df <- b$results$diff_result_df
  ranked <- omicsCore::make_ranked_features(df)

  expect_equal(length(ranked), sum(!is.na(df$effect) & !is.na(df$feature_symbol)))
  # Features nowhere near significance are in it.
  dull <- df$feature_symbol[df$adj_p_value > 0.5]
  expect_true(any(names(ranked) %in% dull))
  # And it is ordered.
  expect_false(is.unsorted(rev(unname(ranked))))
})

test_that("a feature with no effect cannot be ranked and is dropped", {
  b <- omicsApp:::example_diff_bundle()
  df <- b$results$diff_result_df
  n_before <- length(omicsCore::make_ranked_features(df))
  df$effect[1:7] <- NA
  expect_equal(length(omicsCore::make_ranked_features(df)), n_before - 7L)
})

test_that("GSEA mode reports the ranking, not a selection", {
  # "120 of 3000 selected" would describe a step GSEA does not take.
  html <- NULL
  shiny::testServer(enrich_view_server,
                    args = enrich_args(split_bundle(), p_preference = "adjusted"), {
    session$setInputs(type = "gsea")
    session$flushReact()
    s <- selected_features()
    expect_true(isTRUE(s$gsea))
    expect_equal(s$n, 3000L)          # every feature, not the 0 that pass adj.P
    html <<- paste(unlist(output$input_summary), collapse = " ")
  })
  expect_true(grepl("3000 of 3000 features ranked", html, fixed = TRUE))
  expect_true(grepl("do not apply", html, fixed = TRUE))
  expect_false(grepl("nothing to enrich", html, fixed = TRUE))
})

test_that("ORA mode still reports the selection", {
  html <- NULL
  shiny::testServer(enrich_view_server,
                    args = enrich_args(split_bundle(), p_preference = "raw"), {
    session$setInputs(type = "ora")
    session$flushReact()
    expect_false(isTRUE(selected_features()$gsea))
    html <<- paste(unlist(output$input_summary), collapse = " ")
  })
  expect_true(grepl("120 of 3000 features", html, fixed = TRUE))
  expect_false(grepl("ranked", html, fixed = TRUE))
})

# --- what the run actually asks for ---------------------------------------

# The threshold tests above hand the module a synthetic result table, which
# is the right way to test the display rules and the wrong way to find out
# what the run requested. The stored bundle was being cut at the display
# threshold during computation, so switching the table to raw p had nothing
# left to switch to -- and no test that injected its own table could see it.
#
# So this one inspects the call. Running MSigDB enrichment for real belongs
# in omicsCore, where test-output-p-cutoff.R does it; what the view owes is
# the right arguments.
test_that("both enrichment types are asked for the whole pathway table", {
  asked <- list()
  testthat::local_mocked_bindings(
    run_enrichment = function(diff_bundle, type, ...) {
      dots <- list(...)
      asked[[type]] <<- dots
      omicsCore::new_analysis_bundle(
        "run_enrichment",
        params = c(list(type = type), dots),
        results = list(enrich_result_df = data.frame(
          database = "hallmark", result_type = type, comparison = "c",
          pathway_id = "P1", pathway_name = "P1", direction = "up",
          n_genes = 3L, gene_ratio = 0.1, p_value = 0.01,
          adj_p_value = 0.2, enrichment_score = 1, genes = "A/B/C",
          source_label = "x", stringsAsFactors = FALSE)))
    },
    .package = "omicsCore"
  )

  db <- shiny::reactiveVal(split_bundle())
  shiny::testServer(enrich_view_server, args = list(diff_bundle = db), {
    for (type in c("ora", "gsea")) {
      session$setInputs(type = type, database = "hallmark", rerun = 1)
      session$flushReact()
    }
  })

  # ORA needs both numbers, and they are different questions: which
  # features go in, and how much of the table comes back.
  expect_equal(asked$ora$output_p_cutoff, 1)
  expect_false(is.null(asked$ora$p_cutoff))

  # GSEA reads the whole ranked list, so it takes no feature threshold --
  # but it needs the output bound just as much, and that was the half
  # that got missed.
  expect_equal(asked$gsea$output_p_cutoff, 1)
  expect_null(asked$gsea$p_cutoff)
})
