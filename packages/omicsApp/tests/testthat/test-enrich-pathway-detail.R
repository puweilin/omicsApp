# A row of the Enriched sets table opens the pathway under it: GSEA's
# running-score curve, or for ORA -- which has no curve -- the pathway's
# genes that were in the list, with their log2FC.

# A differential result on real Hallmark gene names: inflammatory genes
# up in G2, oxidative phosphorylation down, everything else noise -- so
# both tests find pathways, offline.
pathway_diff <- function() {
  t2g <- omicsCore:::build_term_tables("hallmark", "Hs")$term2gene
  up <- unique(t2g$gene[t2g$term == "HALLMARK_INFLAMMATORY_RESPONSE"])
  down <- setdiff(unique(t2g$gene[t2g$term == "HALLMARK_OXIDATIVE_PHOSPHORYLATION"]), up)
  other <- setdiff(unique(t2g$gene), c(up, down))[1:1500]
  genes <- c(up, down, other)
  set.seed(11)
  n <- 4L
  m <- matrix(stats::rnorm(length(genes) * 2L * n, mean = 20, sd = 0.5),
              nrow = length(genes),
              dimnames = list(genes, paste0("S", seq_len(2L * n))))
  g2 <- seq(n + 1L, 2L * n)
  m[genes %in% up, g2] <- m[genes %in% up, g2] + 1.5
  m[genes %in% down, g2] <- m[genes %in% down, g2] - 1.5
  meta <- data.frame(sample_id = colnames(m), group = rep(c("G1", "G2"), each = n),
                     row.names = colnames(m))
  feat <- data.frame(feature_id = genes, feature_symbol = genes, row.names = genes)
  x <- omicsCore::omics_input(m, meta, feat, omics_type = "proteomics",
                              assay_type = "normalized_intensity")
  omicsCore::run_diff(x, method = "limma", group_col = "group",
                      control_group = "G1", case_group = "G2")
}

pathway_fixture <- local({
  cache <- list()
  function(type) {
    if (is.null(cache[[type]])) {
      d <- pathway_diff()
      e <- suppressWarnings(omicsCore::run_enrichment(
        d, type = type, database = "hallmark",
        direction = if (type == "ora") "separate" else "both"))
      cache[[type]] <<- list(diff = d, enrich = e)
    }
    cache[[type]]
  }
})

detail_html <- function(output) paste(unlist(output$pathway_detail), collapse = " ")

test_that("selecting a GSEA pathway draws its running-score curve", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("msigdbr")
  fx <- pathway_fixture("gsea")
  expect_gt(nrow(fx$enrich$results$enrich_result_df), 0L)
  shiny::testServer(enrich_view_server,
                    args = list(diff_bundle = shiny::reactiveVal(fx$diff)), {
    enrich_bundle(fx$enrich); is_demo(FALSE)
    session$setInputs(show_p = "adjusted", show_cutoff = 1)
    # Nothing selected: no card (the table's header says a click opens
    # one).
    expect_false(grepl("Selected pathway", detail_html(output), fixed = TRUE))
    session$setInputs(hits_rows_selected = 1L)
    expect_match(detail_html(output), "Selected pathway", fixed = TRUE)
    row <- results$selected_pathway()
    first <- enrich_hits_rows(fx$enrich$results$enrich_result_df, "adjusted", 1)
    expect_identical(row$pathway_id, first$pathway_id[[1L]])
    html <- detail_html(output)
    expect_match(html, "gsea_curve", fixed = TRUE)
    expect_match(html, "How to read it", fixed = TRUE)
    p <- results$gsea_curve()
    expect_s3_class(p, "ggplot")
    expect_match(p$labels$title, substr(row$pathway_name, 1L, 20L), fixed = TRUE)
    expect_match(p$labels$subtitle, "NES", fixed = TRUE)
    expect_false(is.null(output$gsea_curve$src))
  })
})

test_that("the curve is drawn from the differential result when the run kept no objects", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("msigdbr")
  fx <- pathway_fixture("gsea")
  bare <- fx$enrich
  bare$results$enrich_object <- NULL
  shiny::testServer(enrich_view_server,
                    args = list(diff_bundle = shiny::reactiveVal(fx$diff)), {
    enrich_bundle(bare); is_demo(FALSE)
    session$setInputs(show_p = "adjusted", show_cutoff = 1, hits_rows_selected = 1L)
    expect_s3_class(results$gsea_curve(), "ggplot")
  })
})

test_that("selecting an ORA pathway lists its genes in the list, with their log2FC", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("msigdbr")
  fx <- pathway_fixture("ora")
  df <- fx$enrich$results$enrich_result_df
  expect_gt(nrow(df), 0L)
  shiny::testServer(enrich_view_server,
                    args = list(diff_bundle = shiny::reactiveVal(fx$diff)), {
    enrich_bundle(fx$enrich); is_demo(FALSE)
    session$setInputs(show_p = "adjusted", show_cutoff = 1, hits_rows_selected = 1L)
    row <- results$selected_pathway()
    html <- detail_html(output)
    expect_match(html, "no", fixed = TRUE)
    expect_match(html, "running score to draw", fixed = TRUE)
    expect_match(html, "log2FC", fixed = TRUE)
    genes <- strsplit(row$overlap_features, "/", fixed = TRUE)[[1L]]
    expect_match(html, genes[[1L]], fixed = TRUE)
    # Each gene's value from the differential result, signed.
    d <- fx$diff$results$diff_result_df
    v <- d$effect[match(genes[[1L]], d$feature_symbol)]
    expect_match(html, sprintf("%+.2f", v), fixed = TRUE)
  })
})

test_that("the ORA genes come without values once the Differential view has moved on", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("msigdbr")
  fx <- pathway_fixture("ora")
  other <- fx$diff
  other$params$comparison <- "something_vs_else"
  shiny::testServer(enrich_view_server,
                    args = list(diff_bundle = shiny::reactiveVal(other)), {
    enrich_bundle(fx$enrich); is_demo(FALSE)
    session$setInputs(show_p = "adjusted", show_cutoff = 1, hits_rows_selected = 1L)
    expect_match(detail_html(output), "another comparison", fixed = TRUE)
  })
})

test_that("a demo pathway says why it has nothing to show", {
  shiny::testServer(enrich_view_server,
                    args = list(diff_bundle = shiny::reactiveVal(NULL)), {
    session$setInputs(show_p = "adjusted", show_cutoff = 1, hits_rows_selected = 1L)
    expect_true(isTRUE(is_demo()))
    expect_match(detail_html(output), "Demo pathways only", fixed = TRUE)
  })
})

test_that("the comparison plot grows with its rows, and more on a phone", {
  expect_equal(compare_plot_px(5L), 480)
  expect_gt(compare_plot_px(17L), 480)
  expect_gt(compare_plot_px(17L, narrow = TRUE), compare_plot_px(17L))
})
