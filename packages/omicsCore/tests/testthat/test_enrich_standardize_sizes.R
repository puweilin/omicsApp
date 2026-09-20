test_that("ORA reports background pathway size separately from query overlap", {
  skip_if_not_installed("clusterProfiler")
  universe <- paste0("g", 1:100)
  # The target has 40 annotated genes, 20 in the universe, and 5 hits.
  terms <- data.frame(
    term = rep(c("target", "background"), c(40, 100)),
    gene = c(paste0("g", c(1:20, 101:120)), universe)
  )
  obj <- clusterProfiler::enricher(
    gene = paste0("g", c(1:5, 51:55)), universe = universe,
    TERM2GENE = terms, pvalueCutoff = 1, qvalueCutoff = 1
  )
  raw <- as.data.frame(obj)
  out <- standardize_enrich_result(obj, "example", "ora", "test")
  target <- out[out$pathway_id == "target", , drop = FALSE]

  expect_equal(nrow(target), 1L)
  expect_equal(target$gene_set_size, 20)
  expect_equal(target$overlap_size, 5)
  expect_identical(out$p_value, raw$pvalue)
  expect_identical(out$adj_p_value, raw$p.adjust)
  expect_identical(out$q_value, raw$qvalue)
})

test_that("ORA can recover both sizes from ratios without Count", {
  raw <- data.frame(ID = "pathway", BgRatio = "20/100", GeneRatio = "5/10")
  out <- standardize_enrich_result(raw, "example", "ora", "test")
  expect_equal(out$gene_set_size, 20)
  expect_equal(out$overlap_size, 5)
})

test_that("missing or unreadable pathway sizes never fall back to hit counts", {
  raw <- data.frame(ID = c("a", "b", "c"), Count = c(5, 6, 7))
  out <- standardize_enrich_result(raw, "example", "ora", "test")
  expect_true(all(is.na(out$gene_set_size)))
  expect_equal(out$overlap_size, raw$Count)

  raw$BgRatio <- c(NA_character_, "garbage", "")
  out <- standardize_enrich_result(raw, "example", "ora", "test")
  expect_true(all(is.na(out$gene_set_size)))
  expect_equal(out$overlap_size, raw$Count)
})

test_that("GSEA retains its setSize and leading genes", {
  raw <- data.frame(
    ID = "pathway", setSize = 40, NES = 1.5,
    core_enrichment = "g1/g2/g3"
  )
  out <- standardize_enrich_result(raw, "example", "gsea", "test")
  expect_equal(out$gene_set_size, 40)
  expect_true(is.na(out$overlap_size))
  expect_identical(out$leading_features, "g1/g2/g3")
})
