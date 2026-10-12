# Text files saved by Excel on Windows are not UTF-8: Windows-1252 in
# Western Europe and the Americas, GBK (GB18030) on Chinese Windows. They
# used to stop the import with "input string 2 is invalid UTF-8". The
# fixtures below are real bytes in those encodings, written from UTF-8
# text through iconv(), and every one is read on both the fread path and
# the read.table fallback.

# Writes `lines` (UTF-8 text) to a file in `encoding`, byte for byte.
encoded_fixture <- function(lines, encoding, ext = ".csv") {
  text <- paste0(paste(enc2utf8(lines), collapse = "\n"), "\n")
  bytes <- if (identical(encoding, "UTF-8")) charToRaw(text)
           else iconv(list(charToRaw(text)), "UTF-8", encoding, toRaw = TRUE)[[1L]]
  path <- tempfile(fileext = ext)
  writeBin(bytes, path)
  path
}

GENES <- c("TP53", "EGFR", "MYC", "KRAS", "BRCA1", "PTEN", "AKT1", "GAPDH",
           "ACTB", "VEGFA", "IL6", "TNF", "CD4", "CD8A", "STAT3", "JUN",
           "FOS", "SOX2", "NOTCH1", "ESR1")

matrix_lines <- function(samples, sep = ",", extra = NULL) {
  set.seed(11)
  vals <- matrix(round(2^stats::rnorm(length(GENES) * length(samples), 18, 1), 1),
                 length(GENES))
  body <- apply(cbind(GENES, vals, extra), 1L, paste, collapse = sep)
  c(paste(c("gene", samples, if (!is.null(extra)) "unit"), collapse = sep), body)
}

read_both_ways <- function(path, ...) {
  lapply(c(fread = TRUE, read.table = FALSE), function(fast) {
    withr::with_options(list(omicsCore.use_fread = fast),
                        suppressWarnings(read_omics(path, orientation = "features_in_rows", ...)))
  })
}

test_that("the encoding is told from the bytes", {
  latin <- c("gene,Café_1", "TP53,1")
  expect_identical(detect_text_encoding(encoded_fixture(c("gene,S1", "TP53,1"), "UTF-8")), "UTF-8")
  expect_identical(detect_text_encoding(encoded_fixture(latin, "UTF-8")), "UTF-8")
  bom <- tempfile(fileext = ".csv")
  writeBin(c(as.raw(c(0xef, 0xbb, 0xbf)), charToRaw(enc2utf8("gene,Café\nTP53,1\n"))), bom)
  expect_identical(detect_text_encoding(bom), "UTF-8")
  expect_identical(detect_text_encoding(encoded_fixture(latin, "CP1252")), "CP1252")
  expect_identical(detect_text_encoding(encoded_fixture(
    c("gene,对照_1,处理_1", "TP53,1,2"), "GB18030")), "GB18030")
  # A byte Windows-1252 leaves undefined: Latin-1, which defines them all.
  odd <- tempfile(fileext = ".csv")
  writeBin(charToRaw("gene,S\x81\nTP53,1\n"), odd)
  expect_identical(detect_text_encoding(odd), "latin1")
})

test_that("Western text is not mistaken for GBK", {
  # "µg" (B5 67) and "°C" (B0 43) are valid GBK pairs, decoding to rare
  # ideographs; validity alone would have read this file as Chinese.
  path <- encoded_fixture(c("gene,unit,note", "TP53,µg,37°C", "EGFR,µl,4°C"),
                          "CP1252")
  expect_true(all(as.raw(c(0xb5, 0xb0)) %in% readBin(path, "raw", 100L)))
  expect_identical(detect_text_encoding(path), "CP1252")
})

test_that("a Windows-1252 matrix file reads, with a note saying so", {
  samples <- c("Café_1", "Café_2", "Ctrl_1", "Ctrl_2")
  path <- encoded_fixture(matrix_lines(samples, extra = rep(c("µg", "°C"), 10)),
                          "CP1252")
  bytes <- readBin(path, "raw", file.size(path))
  expect_true(all(as.raw(c(0xe9, 0xb5, 0xb0)) %in% bytes))
  expect_false(validUTF8(rawToChar(bytes)))
  res <- read_both_ways(path, omics_type = "proteomics", assay_type = "raw_intensity")
  for (r in res) {
    expect_s3_class(r$input, "omics_input")
    expect_identical(colnames(r$input$expr_mat), enc2utf8(samples))
    expect_true(all(validUTF8(colnames(r$input$expr_mat))))
    expect_true("The file was read as Windows-1252 (not UTF-8)." %in% r$report$warnings)
    expect_identical(r$report$suggested_input$encoding, "CP1252")
  }
  expect_identical(res$fread$input$expr_mat, res$read.table$input$expr_mat)
})

test_that("a GBK matrix file reads, with Chinese sample names intact", {
  samples <- c("对照_1", "对照_2", "处理_1", "处理_2")
  path <- encoded_fixture(matrix_lines(samples, sep = "\t"), "GB18030", ext = ".tsv")
  expect_false(validUTF8(rawToChar(readBin(path, "raw", file.size(path)))))
  res <- read_both_ways(path, omics_type = "rnaseq", assay_type = "raw_count")
  for (r in res) {
    expect_identical(colnames(r$input$expr_mat), enc2utf8(samples))
    expect_true("The file was read as Chinese GB18030/GBK (not UTF-8)." %in% r$report$warnings)
  }
  expect_identical(res$fread$input$expr_mat, res$read.table$input$expr_mat)
})

test_that("a UTF-8 file gets no note, and a UTF-16 one is decoded", {
  samples <- c("Café_1", "Café_2", "Ctrl_1", "Ctrl_2")
  utf8 <- encoded_fixture(matrix_lines(samples), "UTF-8")
  r <- suppressWarnings(read_omics(utf8, omics_type = "proteomics",
                                   assay_type = "raw_intensity", orientation = "features_in_rows"))
  expect_false(any(grepl("was read as", r$report$warnings)))
  expect_identical(colnames(r$input$expr_mat), enc2utf8(samples))

  # Excel's "Unicode text" export: UTF-16LE with a byte-order mark.
  text <- paste0(paste(matrix_lines(samples, sep = "\t"), collapse = "\r\n"), "\r\n")
  utf16 <- tempfile(fileext = ".txt")
  writeBin(c(as.raw(c(0xff, 0xfe)),
             iconv(list(charToRaw(enc2utf8(text))), "UTF-8", "UTF-16LE", toRaw = TRUE)[[1L]]),
           utf16)
  r <- suppressWarnings(read_omics(utf16, omics_type = "proteomics",
                                   assay_type = "raw_intensity", orientation = "features_in_rows"))
  expect_identical(colnames(r$input$expr_mat), enc2utf8(samples))
  expect_true("The file was read as UTF-16 (not UTF-8)." %in% r$report$warnings)
})

test_that("a sample sheet in Windows-1252 or GBK is decoded too", {
  samples <- sprintf("S%d", 1:4)
  matrix_path <- encoded_fixture(matrix_lines(samples), "UTF-8")
  cases <- list(
    CP1252 = list(groups = c("Contrôle", "Contrôle", "Traité", "Traité"),
                  note = "The sample sheet was read as Windows-1252 (not UTF-8)."),
    GB18030 = list(groups = c("对照", "对照", "处理", "处理"),
                   note = "The sample sheet was read as Chinese GB18030/GBK (not UTF-8).")
  )
  for (enc in names(cases)) {
    groups <- cases[[enc]]$groups
    sheet <- encoded_fixture(c("sample_id,group", paste(samples, groups, sep = ",")), enc)
    for (fast in c(TRUE, FALSE)) {
      r <- withr::with_options(list(omicsCore.use_fread = fast), suppressWarnings(
        read_omics(matrix_path, omics_type = "proteomics", assay_type = "raw_intensity",
                   orientation = "features_in_rows", sample_sheet = sheet)))
      expect_identical(r$input$meta_df$group, enc2utf8(groups), info = enc)
      expect_true(cases[[enc]]$note %in% r$report$warnings, info = enc)
    }
  }
})

test_that("the decoded copy is removed after reading", {
  path <- encoded_fixture(matrix_lines(c("Café_1", "B", "C", "D")), "CP1252")
  before <- list.files(tempdir(), pattern = "^omics-utf8-")
  suppressWarnings(read_omics(path, omics_type = "proteomics", assay_type = "raw_intensity"))
  expect_identical(list.files(tempdir(), pattern = "^omics-utf8-"), before)
})
