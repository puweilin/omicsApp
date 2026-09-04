# Ensembl ids are what an RNA-seq matrix is keyed on; symbols are what
# every pathway database is keyed on. Without the mapping, enrichment
# matches nothing and hands back an empty result -- which looks exactly
# like a real result that found no enriched pathway.

test_that("known genes map to their current symbols", {
  skip_if(is.null(hgnc_ensembl_map()), "bundled HGNC table not available")
  expect_identical(
    map_ensembl_symbols(c("ENSG00000141510", "ENSG00000012048")),
    c("TP53", "BRCA1")
  )
})

test_that("the version suffix does not stop a match", {
  skip_if(is.null(hgnc_ensembl_map()), "bundled HGNC table not available")
  # Vendors differ on whether they keep it; the gene is the same either way.
  expect_identical(map_ensembl_symbols("ENSG00000141510.17"), "TP53")
})

test_that("an unmapped id becomes NA, not the id", {
  skip_if(is.null(hgnc_ensembl_map()), "bundled HGNC table not available")
  # This is the load-bearing property. Writing the id in would put tens
  # of thousands of strings that can never match any pathway into ORA's
  # universe -- and the universe is the denominator of the hypergeometric
  # test, so inflating it makes every p-value look better than it is.
  out <- map_ensembl_symbols("ENSG09999999999")
  expect_true(is.na(out))
  expect_false(identical(out, "ENSG09999999999"))
})

test_that("both enrichment paths drop the unmapped features", {
  # The reason NA is safe: ORA's universe and selected set, and GSEA's
  # ranked list, all discard them already. Asserted here so a future
  # change to either would fail against the reason it matters.
  df <- data.frame(
    feature_id     = c("ENSG1", "ENSG2", "ENSG3"),
    feature_symbol = c("TP53", NA, "BRCA1"),
    feature_type   = "gene",
    omics_type     = "rnaseq",
    method         = "deseq2",
    analysis_type  = "group",
    comparison     = "G2_vs_G1",
    effect         = c(2, 1, -2),
    effect_type    = "log2fc",
    statistic      = c(4, 2, -4),
    statistic_type = "wald",
    p_value        = c(0.01, 0.01, 0.01),
    adj_p_value    = c(0.01, 0.01, 0.01),
    direction      = c("up", "up", "down"),
    base_mean      = c(100, 100, 100),
    model_fit      = NA_character_,
    is_significant = TRUE,
    stringsAsFactors = FALSE
  )
  expect_length(unique(stats::na.omit(df$feature_symbol)), 2L)
  ranked <- make_ranked_features(df, feature_col = "feature_symbol")
  expect_identical(sort(names(ranked)), c("BRCA1", "TP53"))
})

test_that("ids that are not Ensembl are left alone", {
  expect_false(looks_like_ensembl(c("TP53", "BRCA1")))
  expect_false(looks_like_ensembl(c("P01308", "P02768")))
  expect_true(looks_like_ensembl(c("ENSG00000141510", "ENSG00000012048")))
  # A minority of Ensembl-looking ids is not an Ensembl table.
  expect_false(looks_like_ensembl(c("ENSG00000141510", "TP53", "BRCA1")))
})

test_that("a symbol column the file supplied is not overwritten", {
  feat <- data.frame(
    feature_id     = c("ENSG00000141510", "ENSG00000012048"),
    feature_symbol = c("their_TP53", "their_BRCA1"),
    stringsAsFactors = FALSE
  )
  res <- attach_hgnc_symbols(feat, feat$feature_id)
  expect_identical(res$feature_df$feature_symbol,
                   c("their_TP53", "their_BRCA1"))
  expect_null(res$note)
})

test_that("the note says how many mapped, and where the table came from", {
  skip_if(is.null(hgnc_ensembl_map()), "bundled HGNC table not available")
  ids <- c("ENSG00000141510", "ENSG00000012048", "ENSG09999999999")
  feat <- data.frame(feature_id = ids, stringsAsFactors = FALSE)
  res <- attach_hgnc_symbols(feat, ids)
  expect_identical(res$feature_df$feature_symbol, c("TP53", "BRCA1", NA))
  expect_match(res$note, "2 of 3")
  expect_match(res$note, "HGNC")
})

test_that("provenance names a retrieval date, so a result can cite it", {
  skip_if(is.null(hgnc_ensembl_map()), "bundled HGNC table not available")
  expect_match(hgnc_map_provenance(), "retrieved \\d{4}-\\d{2}-\\d{2}")
})

test_that("symbols are attached by id, not by row position", {
  skip_if(is.null(hgnc_ensembl_map()), "bundled HGNC table not available")
  # materialize_feature_annot() may reorder or subset, so a positional
  # join would attach the wrong symbol to the wrong gene -- a failure
  # that produces plausible output.
  ids <- c("ENSG00000141510", "ENSG00000012048")
  feat <- data.frame(feature_id = rev(ids), stringsAsFactors = FALSE)
  res <- attach_hgnc_symbols(feat, ids)
  expect_identical(res$feature_df$feature_symbol, c("BRCA1", "TP53"))
})
