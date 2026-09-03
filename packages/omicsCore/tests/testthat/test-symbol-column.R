# Proteomics workbooks are keyed by accession, and the gene symbol sits
# in a column next to it. Nothing read that column, so feature_symbol
# fell back to the accession -- and clusterProfiler answered every
# enrichment with "No gene can be mapped", which is not an error and so
# reached the user as an empty result.

annot <- function(symbol_col = "Gene names",
                  symbols = c("NUDT4B", "IGLV4-69", "IGLV8-61")) {
  df <- data.frame(
    Accession   = c("A0A024RBG1", "A0A075B6I0", "A0A075B6I9"),
    .sym        = symbols,
    Description = paste("Protein", 1:3),
    Coverage    = c(12.1, 8.3, 44.0),
    stringsAsFactors = FALSE
  )
  names(df)[2] <- symbol_col
  df
}

materialized <- function(...) {
  a <- annot(...)
  materialize_feature_annot(a, a$Accession)
}

test_that("the accession still becomes feature_id", {
  out <- materialized()
  expect_identical(attr(out, "id_column"), "Accession")
  expect_identical(out$feature_id[[1L]], "A0A024RBG1")
})

test_that("the gene column is found whatever the vendor calls it", {
  # Matched after lowercasing and dropping non-letters, so all of these
  # collapse onto the same entry.
  for (nm in c("Gene names", "Genes", "Gene.Name", "gene_symbol",
               "Symbol", "PG.Genes", "HGNC symbol")) {
    out <- materialized(symbol_col = nm)
    expect_identical(attr(out, "symbol_column"), nm, info = nm)
    expect_identical(out$feature_symbol[[1L]], "NUDT4B", info = nm)
  }
})

test_that("only the leading gene is kept from a packed cell", {
  # MaxQuant writes "NUDT4B;NUDT4"; the whole string matches no gene set.
  out <- materialized(symbols = c("NUDT4B;NUDT4", "IGLV4-69", "A,B"))
  expect_identical(out$feature_symbol, c("NUDT4B", "IGLV4-69", "A"))
})

test_that("an empty symbol becomes NA rather than an empty string", {
  # "" would be carried into a plot label and an enrichment query as a
  # feature named nothing.
  out <- materialized(symbols = c("NUDT4B", "", "  "))
  expect_identical(out$feature_symbol, c("NUDT4B", NA_character_, NA_character_))
})

test_that("an annotation with no gene column is left as it was", {
  # The fallback in diff-standardize.R then sets feature_symbol to
  # feature_id, which is the old behaviour and still the right one when
  # there is nothing better.
  a <- data.frame(Accession = c("A0A024RBG1", "A0A075B6I0"),
                  Coverage = c(1, 2), stringsAsFactors = FALSE)
  out <- materialize_feature_annot(a, a$Accession)
  expect_true(is.na(attr(out, "symbol_column")))
  expect_false("feature_symbol" %in% colnames(out))
})

test_that("an existing feature_symbol column is not overwritten", {
  a <- data.frame(feature_id = c("P1", "P2"),
                  feature_symbol = c("TP53", "EGFR"),
                  Genes = c("WRONG1", "WRONG2"),
                  stringsAsFactors = FALSE)
  out <- materialize_feature_annot(a, a$feature_id)
  expect_identical(out$feature_symbol, c("TP53", "EGFR"))
})

test_that("the id column is never also taken as the symbol", {
  # A workbook whose id column is literally called "Genes" would
  # otherwise have feature_symbol and feature_id set from the same
  # values, which is the fallback dressed up as a real mapping.
  a <- data.frame(Genes = c("TP53", "EGFR"), Coverage = c(1, 2),
                  stringsAsFactors = FALSE)
  out <- materialize_feature_annot(a, a$Genes)
  expect_identical(attr(out, "id_column"), "Genes")
  expect_true(is.na(attr(out, "symbol_column")))
})

test_that("features absent from the annotation still line up", {
  # materialize_feature_annot() left-joins onto the matrix rows; the
  # padded rows must carry NA symbols, not be dropped or misaligned.
  a <- annot()
  out <- materialize_feature_annot(a, c(a$Accession, "P99999"))
  expect_equal(nrow(out), 4L)
  expect_identical(out$feature_id[[4L]], "P99999")
  expect_true(is.na(out$feature_symbol[[4L]]))
})

test_that("symbols reach the differential result, which is the point", {
  set.seed(3L)
  ids <- c("A0A024RBG1", "A0A075B6I0", "A0A075B6I9")
  m <- matrix(rnorm(3 * 6, 20, 1), 3, 6,
              dimnames = list(ids, sprintf("S%d", 1:6)))
  meta <- data.frame(condition = rep(c("G1", "G2"), each = 3),
                     row.names = colnames(m))
  feat <- materialize_feature_annot(annot(), ids)
  inp <- omics_input(m, meta, feat, omics_type = "proteomics",
                     assay_type = "normalized_intensity")
  b <- suppressWarnings(run_diff(inp, method = "ttest", analysis_type = "group",
                                 group_col = "condition",
                                 control_group = "G1", case_group = "G2"))
  expect_setequal(b$results$diff_result_df$feature_symbol,
                  c("NUDT4B", "IGLV4-69", "IGLV8-61"))
})
