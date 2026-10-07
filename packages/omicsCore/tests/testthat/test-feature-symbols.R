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
  res <- attach_gene_symbols(feat, feat$feature_id)
  expect_identical(res$feature_df$feature_symbol,
                   c("their_TP53", "their_BRCA1"))
  expect_null(res$note)
})

test_that("the note says how many mapped, and where the table came from", {
  skip_if(is.null(hgnc_ensembl_map()), "bundled HGNC table not available")
  ids <- c("ENSG00000141510", "ENSG00000012048", "ENSG09999999999")
  feat <- data.frame(feature_id = ids, stringsAsFactors = FALSE)
  res <- attach_gene_symbols(feat, ids)
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
  res <- attach_gene_symbols(feat, ids)
  expect_identical(res$feature_df$feature_symbol, c("BRCA1", "TP53"))
})

# ---- mouse -----------------------------------------------------------------

test_that("mouse Ensembl ids map to MGI symbols, versioned or not", {
  skip_if_not_installed("babelgene")
  expect_identical(map_ensembl_symbols(c("ENSMUSG00000059552", "ENSMUSG00000059552.14")),
                   c("Trp53", "Trp53"))
  # A gene without a human counterpart, or no gene at all, gets NA.
  expect_true(is.na(map_ensembl_symbols("ENSMUSG09999999999")))
})

test_that("human and mouse ids in one vector each map from their own table", {
  skip_if_not_installed("babelgene")
  skip_if(is.null(hgnc_ensembl_map()), "bundled HGNC table not available")
  expect_identical(map_ensembl_symbols(c("ENSG00000141510", "ENSMUSG00000059552", "P04637")),
                   c("TP53", "Trp53", NA))
})

test_that("a mouse table is recognised, and the note names MGI and the species", {
  skip_if_not_installed("babelgene")
  expect_identical(ensembl_species(c("ENSMUSG00000059552", "ENSMUSG00000000001.4")), "mouse")
  expect_identical(ensembl_species(c("ENSG00000141510", "ENSG00000012048")), "human")
  expect_true(is.na(ensembl_species(c("TP53", "Trp53"))))
  expect_true(looks_like_ensembl(c("ENSMUSG00000059552", "ENSMUSG00000000001")))
  ids <- c("ENSMUSG00000059552", "ENSMUSG00000000001", "ENSMUSG09999999999")
  res <- attach_gene_symbols(data.frame(feature_id = ids), ids)
  expect_identical(res$feature_df$feature_symbol, c("Trp53", "Gnai3", NA))
  expect_match(res$note, "2 of 3 features matched an MGI symbol", fixed = TRUE)
  expect_match(res$note, "Choose mouse as the species in Enrichment", fixed = TRUE)
  expect_match(res$note, "babelgene")
})

test_that("without babelgene, mouse ids keep no symbol and the note says why", {
  local_mocked_bindings(is_installed = function(pkg) !identical(pkg, "babelgene"))
  ids <- c("ENSMUSG00000059552", "ENSMUSG00000000001")
  res <- attach_gene_symbols(data.frame(feature_id = ids), ids)
  expect_null(res$feature_df$feature_symbol)
  expect_match(res$note, "babelgene) is not installed", fixed = TRUE)
  expect_identical(map_ensembl_symbols(ids), c(NA_character_, NA_character_))
})

test_that("a mouse counts table is imported with MGI symbols", {
  skip_if_not_installed("babelgene")
  dir <- withr::local_tempdir()
  ids <- c("ENSMUSG00000059552.14", "ENSMUSG00000000001.4", sprintf("ENSMUSG%011d.1", 900000 + 1:6))
  set.seed(5)
  m <- matrix(stats::rpois(8 * 4, 30), 8, dimnames = list(ids, paste0("S", 1:4)))
  path <- file.path(dir, "mouse.csv")
  utils::write.csv(data.frame(gene_id = ids, m, check.names = FALSE), path, row.names = FALSE)
  res <- read_omics(path, omics_type = "rnaseq", assay_type = "raw_count")
  expect_identical(res$input$feature_df$feature_symbol[1:2], c("Trp53", "Gnai3"))
  # The ids themselves are left as the file wrote them.
  expect_identical(rownames(res$input$expr_mat), ids)
  expect_true(any(grepl("MGI symbol", res$report$warnings)))
})

# ---- GENCODE's _PAR_Y ------------------------------------------------------

test_that("the _PAR_Y copy of a gene maps to the gene it copies", {
  skip_if(is.null(hgnc_ensembl_map()), "bundled HGNC table not available")
  x <- map_ensembl_symbols("ENSG00000182378.15")
  expect_false(is.na(x))
  expect_identical(map_ensembl_symbols(c("ENSG00000182378.15_PAR_Y", "ENSG00000182378_PAR_Y")),
                   c(x, x))
  expect_true(looks_like_ensembl(c("ENSG00000182378.15_PAR_Y", "ENSG00000141510.17")))
})

test_that("version stripping keeps the _PAR_Y copy apart from its X twin", {
  # tx2gene matching compares unversioned ids; the copies are different
  # rows there and must stay so.
  expect_identical(strip_ensembl_version(c("ENST00000381192.10_PAR_Y", "ENST00000381192.10",
                                           "ENSG00000141510")),
                   c("ENST00000381192_PAR_Y", "ENST00000381192", "ENSG00000141510"))
  expect_identical(ensembl_base_id("ENSG00000182378.15_PAR_Y"), "ENSG00000182378")
})

test_that("empty _PAR_Y rows are dropped, and ones with values kept apart", {
  m <- rbind(a = c(1, 2), b_PAR_Y = c(0, 0), c_PAR_Y = c(NA, 0), d_PAR_Y = c(0, 3))
  res <- par_y_rows(m)
  expect_identical(res$keep, c(TRUE, FALSE, FALSE, TRUE))
  expect_length(res$notes, 2L)
  expect_match(res$notes[[1L]], "Dropped 2 row(s) whose ids end in '_PAR_Y'", fixed = TRUE)
  expect_match(res$notes[[2L]], "Kept 1 row(s)", fixed = TRUE)
  # Nothing to do, nothing said.
  expect_identical(par_y_rows(m[1, , drop = FALSE])$notes, character(0))
  # Never every row.
  expect_true(all(par_y_rows(m[2:3, , drop = FALSE])$keep))
})

test_that("a GENCODE counts table imports with its _PAR_Y rows handled", {
  skip_if(is.null(hgnc_ensembl_map()), "bundled HGNC table not available")
  dir <- withr::local_tempdir()
  ids <- c("ENSG00000141510.17", "ENSG00000012048.23",
           "ENSG00000182378.15", "ENSG00000182378.15_PAR_Y",   # X copy, empty Y copy
           "ENSG00000178605.13", "ENSG00000178605.13_PAR_Y",   # Y copy with reads
           sprintf("ENSG%011d.1", 900000 + 1:4))
  set.seed(6)
  m <- matrix(stats::rpois(length(ids) * 4, 30) + 1, length(ids),
              dimnames = list(ids, paste0("S", 1:4)))
  m["ENSG00000182378.15_PAR_Y", ] <- 0
  path <- file.path(dir, "gencode.tsv")
  utils::write.table(data.frame(gene_id = ids, m, check.names = FALSE), path, sep = "\t",
                     quote = FALSE, row.names = FALSE)
  res <- read_omics(path, omics_type = "rnaseq", assay_type = "raw_count")
  inp <- res$input
  expect_false(is.null(inp))
  expect_false("ENSG00000182378.15_PAR_Y" %in% rownames(inp$expr_mat))
  expect_true("ENSG00000178605.13_PAR_Y" %in% rownames(inp$expr_mat))
  expect_identical(nrow(inp$expr_mat), length(ids) - 1L)
  expect_false(anyDuplicated(inp$feature_df$feature_id) > 0L)
  # One gene, one symbol: the X copy carries it.
  sym <- stats::setNames(inp$feature_df$feature_symbol, inp$feature_df$feature_id)
  expect_false(is.na(sym[["ENSG00000178605.13"]]))
  expect_true(is.na(sym[["ENSG00000178605.13_PAR_Y"]]))
  expect_true(any(grepl("Dropped 1 row(s) whose ids end in '_PAR_Y'", res$report$warnings,
                        fixed = TRUE)))
  expect_true(any(grepl("Kept 1 row(s) whose ids end in '_PAR_Y'", res$report$warnings,
                        fixed = TRUE)))
})

test_that("ids that collide once a pipeline has stripped the suffix still import", {
  # Some pipelines strip versions and the _PAR_Y suffix themselves, which
  # leaves the X and Y copies under one id.
  dir <- withr::local_tempdir()
  ids <- c("ENSG00000182378", "ENSG00000182378", "ENSG00000141510", "ENSG00000012048",
           sprintf("ENSG%011d", 900000 + 1:4))
  set.seed(7)
  m <- matrix(stats::rpois(length(ids) * 4, 30), length(ids))
  path <- file.path(dir, "stripped.csv")
  utils::write.csv(data.frame(gene_id = ids, S1 = m[, 1], S2 = m[, 2], S3 = m[, 3],
                              S4 = m[, 4]), path, row.names = FALSE)
  res <- read_omics(path, omics_type = "rnaseq", assay_type = "raw_count")
  expect_false(is.null(res$input))
  expect_identical(rownames(res$input$expr_mat)[1:2], c("ENSG00000182378", "ENSG00000182378_1"))
})

test_that("quantification files drop their empty _PAR_Y rows, lengths included", {
  dir <- withr::local_tempdir()
  ids <- c(sprintf("ENST%011d.1", 1:6), "ENST00000381192.10", "ENST00000381192.10_PAR_Y")
  paths <- vapply(c("a_1", "b_1"), function(s) {
    p <- file.path(dir, paste0(s, ".quant.sf"))
    reads <- c(1:7 * 10 + 0.5, 0)
    utils::write.table(data.frame(Name = ids, Length = 1000, EffectiveLength = 850,
                                  TPM = reads / sum(reads) * 1e6, NumReads = reads),
                       p, sep = "\t", quote = FALSE, row.names = FALSE)
    p
  }, character(1))
  res <- read_quant_files(unname(paths))
  expect_false("ENST00000381192.10_PAR_Y" %in% rownames(res$input$expr_mat))
  expect_identical(rownames(res$input$misc$tximport$length), rownames(res$input$expr_mat))
  expect_true(any(grepl("_PAR_Y", res$report$warnings)))
})
