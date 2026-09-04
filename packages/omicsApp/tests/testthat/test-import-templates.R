# A template that does not import is worse than no template: it teaches a
# shape the reader then has to debug. These round-trip it through the real
# reader rather than checking the sheets look right.

skip_if_no_xlsx <- function() {
  testthat::skip_if_not_installed("openxlsx")
  testthat::skip_if_not_installed("readxl")
}

test_that("each template has the sheets that kind of vendor actually ships", {
  skip_if_no_xlsx()
  data_sheets <- function(t) setdiff(names(import_template_sheets(t)), "README")

  # Two for RNA-seq: symbols are derived from the Ensembl id, so no
  # feature sheet is needed. Three for proteomics: a UniProt accession
  # has no such mapping here, so the symbol has to come from the file.
  expect_identical(data_sheets("rnaseq"), c("expression", "sample_info"))
  expect_identical(data_sheets("proteomics"),
                   c("expression", "sample_info", "feature_info"))
})

test_that("both templates import cleanly and end up with symbols", {
  skip_if_no_xlsx()
  for (t in c("proteomics", "rnaseq")) {
    path <- withr::local_tempfile(fileext = ".xlsx")
    write_import_template(path, t)

    res <- omicsCore::read_omics(
      path, omics_type = t,
      assay_type = if (t == "rnaseq") "raw_count" else "raw_intensity")

    expect_false(is.null(res$input), info = t)
    expect_equal(ncol(res$input$expr_mat), 6L, info = t)
    # Every example feature resolves to a symbol -- for RNA-seq that is
    # the mapping working, which is half of what the template shows.
    expect_true(all(!is.na(res$input$feature_df$feature_symbol)), info = t)
    expect_true(all(c("donor", "condition") %in%
                      colnames(res$input$meta_df)), info = t)
  }
})

test_that("the README sheet is ignored rather than mistaken for data", {
  skip_if_no_xlsx()
  path <- withr::local_tempfile(fileext = ".xlsx")
  write_import_template(path, "rnaseq")
  res <- omicsCore::read_omics(path, omics_type = "rnaseq",
                               assay_type = "raw_count")
  sheets <- res$report$sheets
  expect_identical(sheets$role[sheets$name == "README"], "unknown")
})

test_that("the two templates share donors but not sample ids", {
  skip_if_no_xlsx()
  # This is the whole integration story: the same six people, measured
  # twice, under two sets of sample names. A template where the ids
  # matched would teach that they have to.
  p <- import_template_sheets("proteomics")$sample_info
  r <- import_template_sheets("rnaseq")$sample_info

  expect_identical(p$donor, r$donor)
  expect_length(intersect(p$sample_id, r$sample_id), 0L)
})

test_that("the RNA-seq template's ids are real enough to map", {
  skip_if_no_xlsx()
  # A made-up ENSG would map to nothing and teach the opposite of what
  # the template is for.
  ids <- import_template_sheets("rnaseq")$expression$feature_id
  expect_true(all(!is.na(omicsCore::map_ensembl_symbols(ids))))
})

test_that("the download handlers produce importable workbooks", {
  skip_if_no_xlsx()
  shiny::testServer(import_view_server, {
    session$setInputs(omics_type = "rnaseq")
    path <- output$template_rnaseq
    expect_true(file.exists(path))
    expect_gt(file.size(path), 0)
    expect_true("sample_info" %in% readxl::excel_sheets(path))
  })
})
