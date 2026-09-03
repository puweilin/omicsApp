# Which pathways count as significant is a reading of a finished result,
# not an input to it. It used to be neither choosable nor separable:
# run_enrichment()'s p_cutoff bounded the stored table *and*, for ORA,
# selected the features, so a bundle could only hold what one number
# admitted and there was nothing to switch between.

gsea_frame <- function(n = 200L, n_adj = 9L) {
  set.seed(5L)
  data.frame(
    database    = "go_bp",
    result_type = "gsea",
    comparison  = "G2_vs_G1",
    pathway_id  = sprintf("GO%04d", seq_len(n)),
    pathway_name = sprintf("PATHWAY %d", seq_len(n)),
    effect      = rnorm(n, 0, 1.5),
    effect_type = "nes",
    direction   = "both",
    # 9 survive adjustment, 60 pass at raw p -- the shape of the real
    # result this came from.
    p_value     = c(runif(60L, 0, 0.049), runif(n - 60L, 0.05, 1)),
    adj_p_value = c(runif(n_adj, 0, 0.049), runif(n - n_adj, 0.05, 1)),
    q_value     = NA_real_,
    gene_set_size = 100L, overlap_size = 20L,
    overlap_features = NA_character_, leading_features = NA_character_,
    source_label = "gsea_both_go_bp",
    stringsAsFactors = FALSE
  )
}

gsea_bundle <- function(...) {
  omicsCore::new_analysis_bundle(
    "run_enrichment", input_info = list(omics_type = "proteomics"),
    params = list(type = "gsea", database = "go_bp", organism = "Hs",
                  direction = "both", comparison = "G2_vs_G1"),
    results = list(enrich_result_df = gsea_frame(...)))
}

view_args <- function(bundle) {
  list(diff_bundle = shiny::reactiveVal(NULL),
       diff_thresholds = shiny::reactive(list(
         p_cutoff = 0.05, p_preference = "adjusted", effect_cutoff = NULL)))
}

test_that("the display threshold selects on the column it names", {
  df <- gsea_frame()
  expect_equal(nrow(omicsCore::filter_enrich_results(
    df, p_cutoff = 0.05, p_preference = "adjusted")), 9L)
  expect_equal(nrow(omicsCore::filter_enrich_results(
    df, p_cutoff = 0.05, p_preference = "raw")), 60L)
})

test_that("the table is sorted and labelled by the column in use", {
  b <- gsea_bundle()
  for (pref in c("adjusted", "raw")) {
    html <- NULL
    shiny::testServer(enrich_view_server, args = view_args(b), {
      enrich_bundle(b); is_demo(FALSE)
      session$setInputs(show_p = pref, show_cutoff = 0.05)
      session$flushReact()
      html <<- paste(unlist(output$hits), collapse = " ")
    })
    if (identical(pref, "raw")) {
      expect_true(grepl(">p<", html, fixed = TRUE), info = pref)
    } else {
      expect_true(grepl("adj.P", html, fixed = TRUE), info = pref)
    }
  }
})

test_that("an empty or nonsensical cutoff falls back rather than emptying", {
  # numericInput reports NA mid-edit, and NA in the comparison hides
  # every pathway with nothing on screen to explain it.
  b <- gsea_bundle()
  shiny::testServer(enrich_view_server, args = view_args(b), {
    session$setInputs(show_cutoff = NA_real_)
    expect_equal(show_cutoff(), 0.05)
    session$setInputs(show_cutoff = 0)
    expect_equal(show_cutoff(), 0.05)
    session$setInputs(show_cutoff = 2)
    expect_equal(show_cutoff(), 0.05)
    session$setInputs(show_cutoff = 0.1)
    expect_equal(show_cutoff(), 0.1)
  })
})

test_that("adjusted is the default, so nothing is loosened by accident", {
  shiny::testServer(enrich_view_server, args = view_args(gsea_bundle()), {
    expect_identical(show_p(), "adjusted")
    expect_equal(show_cutoff(), 0.05)
  })
})
