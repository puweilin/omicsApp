# The import view takes several Salmon / RSEM / kallisto files at once and
# merges them into one layer; reads Windows-1252 and GBK text; and shows
# the value scale it inferred, with the reason, before anything is
# imported.

write_salmon <- function(path, ids, reads, efflen = rep(1000, length(ids))) {
  tpm <- reads / efflen
  utils::write.table(data.frame(Name = ids, Length = efflen + 150, EffectiveLength = efflen,
                                TPM = tpm / sum(tpm) * 1e6, NumReads = reads),
                     path, sep = "\t", quote = FALSE, row.names = FALSE)
  path
}

# As Shiny hands several uploads over: one row per file, each in its own
# temporary file named by position.
uploads <- function(paths, names) {
  data.frame(name = names, size = file.size(paths), type = "",
             datapath = paths, stringsAsFactors = FALSE)
}

test_that("several quantification files import as one RNA-seq layer", {
  dir <- withr::local_tempdir()
  set.seed(3)
  ids <- sprintf("ENST%011d", 1:50)
  names <- c("ctrl_1.quant.sf", "ctrl_2.quant.sf", "ko_1.quant.sf", "ko_2.quant.sf")
  paths <- vapply(seq_along(names), function(i) {
    write_salmon(file.path(dir, sprintf("%d.sf", i - 1L)), ids,
                 stats::rpois(50, 80) + 0.5, stats::runif(50, 500, 2000))
  }, character(1))

  shiny::testServer(import_view_server, {
    session$setInputs(omics_type = "rnaseq", file = uploads(paths, names))
    expect_true(parse_ok())
    inp <- parsed()$input
    expect_identical(colnames(inp$expr_mat), c("ctrl_1", "ctrl_2", "ko_1", "ko_2"))
    expect_identical(inp$assay_type, "raw_count")
    expect_false(is.null(inp$misc$tximport$length))
    expect_true(any(grepl("Read 4 Salmon quantification file", parsed()$report$warnings)))
    expect_match(as.character(output$upload_status$html), "4 files")

    session$setInputs(confirm = 1)
    done <- session$returned()
    expect_s3_class(done, "omics_input")
    # The lengths travel with the layer into the project.
    expect_identical(dim(done$misc$tximport$length), c(50L, 4L))
    expect_match(done$source_fingerprint, "\\+")
  })
})

test_that("one quantification file is one sample, named after the upload", {
  dir <- withr::local_tempdir()
  path <- write_salmon(file.path(dir, "0.sf"), paste0("T", 1:20), 1:20 + 0.25)
  shiny::testServer(import_view_server, {
    session$setInputs(omics_type = "rnaseq",
                      file = list(datapath = path, name = "wt_3.quant.sf",
                                  size = file.size(path)))
    expect_true(parse_ok())
    expect_identical(colnames(parsed()$input$expr_mat), "wt_3")
    expect_true(any(grepl("one sample cannot be compared", parsed()$report$warnings)))
  })
})

test_that("several files that are not quantification files are refused in words", {
  dir <- withr::local_tempdir()
  a <- file.path(dir, "0.csv")
  b <- file.path(dir, "1.csv")
  writeLines(c("gene,S1,S2", "TP53,1,2", "EGFR,3,4"), a)
  writeLines(c("gene,S3,S4", "TP53,1,2", "EGFR,3,4"), b)
  shiny::testServer(import_view_server, {
    session$setInputs(omics_type = "rnaseq",
                      file = uploads(c(a, b), c("counts_a.csv", "counts_b.csv")))
    expect_false(parse_ok())
    w <- parsed()$report$warnings
    expect_true(any(grepl("Several files can be imported together only", w)))
    expect_true(any(grepl("'counts_a.csv' is not a Salmon", w, fixed = TRUE)))
  })
})

test_that("a TPM table is labelled TPM, with the reason, and not as counts", {
  set.seed(8)
  genes <- sprintf("ENSG%011d", 1:200)
  counts <- matrix(stats::rpois(200 * 4, 100), 200)
  rate <- counts / stats::runif(200, 500, 3000)
  tpm <- sweep(rate, 2L, colSums(rate), "/") * 1e6
  path <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(gene_id = genes, A_1 = tpm[, 1], A_2 = tpm[, 2],
                              B_1 = tpm[, 3], B_2 = tpm[, 4]), path, row.names = FALSE)
  shiny::testServer(import_view_server, {
    session$setInputs(omics_type = "rnaseq",
                      file = list(datapath = path, name = "tpm.csv", size = file.size(path)))
    expect_true(parse_ok())
    expect_identical(parsed()$input$assay_type, "tpm")
    picker <- as.character(output$assay_type_picker$html)
    expect_match(picker, "Guessed from the values")
    expect_match(picker, "adds up to about a million")
    expect_match(picker, "DESeq2 and edgeR need read counts")
    expect_false("deseq2" %in% omicsCore::applicable_diff_methods(parsed()$input))
    # The user can still overrule it.
    session$setInputs(assay_type = "fpkm")
    expect_identical(parsed()$input$assay_type, "fpkm")
    expect_match(as.character(output$assay_type_picker$html), "you changed it")
  })
})

test_that("a Windows-1252 CSV imports, and the page says how it was read", {
  genes <- c("TP53", "EGFR", "MYC", "KRAS", "BRCA1", "PTEN", "AKT1", "GAPDH",
             "ACTB", "VEGFA", "IL6", "TNF")
  set.seed(4)
  vals <- matrix(round(2^stats::rnorm(12 * 4, 18, 1)), 12)
  lines <- c("gene,Café_1,Café_2,Thé_1,Thé_2",
             apply(cbind(genes, vals), 1L, paste, collapse = ","))
  path <- tempfile(fileext = ".csv")
  writeBin(iconv(list(charToRaw(enc2utf8(paste0(paste(lines, collapse = "\r\n"), "\r\n")))),
                 "UTF-8", "CP1252", toRaw = TRUE)[[1L]], path)
  shiny::testServer(import_view_server, {
    session$setInputs(omics_type = "proteomics",
                      file = list(datapath = path, name = "excel.csv", size = file.size(path)))
    expect_true(parse_ok())
    expect_identical(colnames(parsed()$input$expr_mat),
                     enc2utf8(c("Café_1", "Café_2", "Thé_1", "Thé_2")))
    expect_match(as.character(output$schema_warnings$html),
                 "The file was read as Windows-1252 (not UTF-8).", fixed = TRUE)
  })
})
