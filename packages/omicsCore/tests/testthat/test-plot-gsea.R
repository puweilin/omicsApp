# plot_gsea(): the running-score curve of one pathway, drawn by omicsCore
# itself (it was enrichplot::gseaplot2(), in bright green and three
# panels whose axis titles ran together at phone width).

gsea_fixture <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      d <- realistic_diff_bundle()
      g <- suppressWarnings(run_enrichment(d, type = "gsea", database = "hallmark"))
      cache <<- list(diff = d, gsea = g)
    }
    cache
  }
})

test_that("the running score is GSEA's, as clusterProfiler computes it", {
  stats <- c(a = 3, b = 2, c = 1, d = -1, e = -2, f = -3)
  rs <- gsea_running_score(stats, c("a", "c"))
  # Hits add their share of the summed |stat| (3/4, 1/4); misses take
  # 1/4 each.
  expect_equal(rs$score, c(0.75, 0.5, 0.75, 0.5, 0.25, 0))
  expect_identical(rs$hit, c(TRUE, FALSE, TRUE, FALSE, FALSE, FALSE))

  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("DOSE")
  fx <- gsea_fixture()
  obj <- fx$gsea$results$enrich_object[[1L]]
  id <- fx$gsea$results$enrich_result_df$pathway_id[[1L]]
  ours <- gsea_running_score(obj@geneList, obj@geneSets[[id]])
  ns <- asNamespace("DOSE")
  skip_if_not(exists("gseaScores", envir = ns, inherits = FALSE))
  theirs <- get("gseaScores", envir = ns)(obj@geneList, obj@geneSets[[id]], fortify = TRUE)
  expect_equal(ours$score, theirs$runningScore, tolerance = 1e-10)
  # Its extreme is the enrichment score fgsea reported.
  es <- obj@result$enrichmentScore[obj@result$ID == id]
  expect_equal(ours$score[which.max(abs(ours$score))], es, tolerance = 1e-6)
})

test_that("plot_gsea titles the pathway, says NES and adjusted p, and colours by direction", {
  skip_if_not_installed("clusterProfiler")
  fx <- gsea_fixture()
  df <- fx$gsea$results$enrich_result_df
  row <- df[which.max(df$effect), ]
  p <- plot_gsea(fx$gsea, row$pathway_id)
  expect_s3_class(p, "ggplot")
  expect_identical(p$labels$title, wrap_label(row$pathway_name, width = 40L, max_lines = 2L))
  expect_match(p$labels$subtitle, sprintf("NES %.2f", row$effect), fixed = TRUE)
  expect_match(p$labels$subtitle, "adjusted p", fixed = TRUE)
  expect_identical(p$labels$y, "running enrichment score")
  expect_identical(p$labels$x, "genes ranked by test statistic")
  b <- ggplot2::ggplot_build(p)
  line <- b$data[[which(vapply(p$layers, function(l) inherits(l$geom, "GeomLine"), logical(1)))]]
  expect_identical(unique(line$colour), omics_colors$up)
  # The strip under the ticks names the two ends of the ranking.
  txt <- unlist(lapply(b$data, function(d) d$label))
  expect_true(all(c("up", "down") %in% txt))
  # By name as well as by ID.
  expect_s3_class(plot_gsea(fx$gsea, row$pathway_name), "ggplot")
})

test_that("plot_gsea redraws from the differential result when the objects are gone", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("msigdbr")
  fx <- gsea_fixture()
  id <- fx$gsea$results$enrich_result_df$pathway_id[[1L]]
  stored <- plot_gsea(fx$gsea, id)
  bare <- fx$gsea
  bare$results$enrich_object <- NULL
  expect_error(plot_gsea(bare, id), "diff_bundle")
  rebuilt <- plot_gsea(bare, id, diff_bundle = fx$diff)
  expect_equal(rebuilt$data$score, stored$data$score, tolerance = 1e-10)
  expect_error(plot_gsea(fx$gsea, "NOT_A_PATHWAY"), "No gene set")
})

test_that("plot_gsea refuses an ORA bundle", {
  b <- new_analysis_bundle("run_enrichment", params = list(type = "ora"),
                           results = list(enrich_result_df = data.frame()))
  expect_error(plot_gsea(b, "x"), "type = 'gsea'")
})
