# The import view takes compressed text (counts.csv.gz), macro-enabled
# workbooks (.xlsm) and SummarizedExperiment / DESeqDataSet objects saved
# as .rds.

test_that("the upload fields accept .gz and .xlsm", {
  html <- as.character(import_upload_card(shiny::NS("x")))
  expect_match(html, "accept=\"[^\"]*\\.gz[^\"]*\\.rds")
  expect_match(html, "accept=\"[^\"]*\\.xlsm")
  # The sample sheet field too.
  expect_match(html, "\\.xlsm,\\.xls,\\.csv,\\.tsv,\\.txt,\\.gz\"")
})

test_that("a .csv.gz upload imports, whatever Shiny names its copy", {
  set.seed(5)
  m <- matrix(stats::rpois(10 * 4, 50), 10)
  plain <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(gene = sprintf("G%02d", 1:10), ctrl_1 = m[, 1], ctrl_2 = m[, 2],
                              ko_1 = m[, 3], ko_2 = m[, 4]), plain, row.names = FALSE)
  # Shiny keeps only the last extension: "counts.csv.gz" arrives as "0.gz".
  path <- file.path(withr::local_tempdir(), "0.gz")
  con <- gzfile(path, "wb")
  writeBin(readBin(plain, "raw", file.size(plain)), con)
  close(con)
  shiny::testServer(import_view_server, {
    session$setInputs(omics_type = "rnaseq",
                      file = list(datapath = path, name = "counts.csv.gz",
                                  size = file.size(path)))
    expect_true(parse_ok())
    expect_identical(colnames(parsed()$input$expr_mat), c("ctrl_1", "ctrl_2", "ko_1", "ko_2"))
    expect_identical(parsed()$input$assay_type, "raw_count")
  })
})

test_that("a SummarizedExperiment of read counts imports as RNA-seq counts", {
  skip_if_not_installed("SummarizedExperiment")
  set.seed(6)
  m <- matrix(stats::rnbinom(80 * 4, mu = 200, size = 5), 80,
              dimnames = list(sprintf("ENSG%011d", 1:80), c("ctrl_1", "ctrl_2", "ko_1", "ko_2")))
  se <- SummarizedExperiment::SummarizedExperiment(
    list(counts = m), colData = S4Vectors::DataFrame(group = c("ctrl", "ctrl", "ko", "ko"),
                                                     row.names = colnames(m)))
  path <- file.path(withr::local_tempdir(), "0.rds")
  saveRDS(se, path)
  shiny::testServer(import_view_server, {
    # The radio still says proteomics; the counts decide.
    session$setInputs(omics_type = "proteomics",
                      file = list(datapath = path, name = "experiment.rds",
                                  size = file.size(path)))
    expect_true(parse_ok())
    inp <- parsed()$input
    expect_identical(inp$omics_type, "rnaseq")
    expect_identical(inp$assay_type, "raw_count")
    expect_identical(inp$meta_df$group, c("ctrl", "ctrl", "ko", "ko"))
    expect_match(as.character(output$schema_warnings$html),
                 "read as RNA-seq, not proteomics", fixed = TRUE)
  })
})
