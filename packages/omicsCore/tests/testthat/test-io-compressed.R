# Gzip-compressed text (counts.csv.gz, quant.sf.gz) and macro-enabled
# workbooks (.xlsm): read like the plain files they wrap.

gzip_file <- function(src, dest) {
  con <- gzfile(dest, "wb")
  writeBin(readBin(src, "raw", file.size(src)), con)
  close(con)
  dest
}

write_counts_csv <- function(path, sep = ",") {
  set.seed(11)
  m <- matrix(stats::rpois(8 * 4, 40), 8,
              dimnames = list(sprintf("G%02d", 1:8), c("ctrl_1", "ctrl_2", "ko_1", "ko_2")))
  utils::write.table(data.frame(gene = rownames(m), m, check.names = FALSE), path,
                     sep = sep, quote = FALSE, row.names = FALSE)
  path
}

test_that("a .csv.gz, .tsv.gz and .txt.gz read exactly as the plain file", {
  dir <- withr::local_tempdir()
  for (ext in c("csv", "tsv", "txt")) {
    plain <- write_counts_csv(file.path(dir, paste0("counts.", ext)),
                              sep = if (ext == "csv") "," else "\t")
    gz <- gzip_file(plain, paste0(plain, ".gz"))
    a <- read_omics(plain, omics_type = "rnaseq", assay_type = "raw_count")
    b <- read_omics(gz, omics_type = "rnaseq", assay_type = "raw_count")
    expect_identical(b$input$expr_mat, a$input$expr_mat)
    expect_identical(b$input$meta_df, a$input$meta_df)
    # The sheet is named after the file, without either extension.
    expect_identical(b$report$sheets$name, "counts")
  }
})

test_that("read.table's path reads compressed text too", {
  withr::local_options(omicsCore.use_fread = FALSE)
  dir <- withr::local_tempdir()
  plain <- write_counts_csv(file.path(dir, "counts.csv"))
  gz <- gzip_file(plain, file.path(dir, "counts.csv.gz"))
  b <- read_omics(gz, omics_type = "rnaseq", assay_type = "raw_count")
  expect_identical(dim(b$input$expr_mat), c(8L, 4L))
  # Arguments for read.table also take the slow path.
  c2 <- read_omics(gz, omics_type = "rnaseq", assay_type = "raw_count", strip.white = TRUE)
  expect_identical(c2$input$expr_mat, b$input$expr_mat)
})

test_that("a browser upload named only '0.gz' is read as text", {
  dir <- withr::local_tempdir()
  plain <- write_counts_csv(file.path(dir, "counts.csv"))
  gz <- gzip_file(plain, file.path(dir, "0.gz"))
  expect_identical(detect_file_type(gz), "csv")
  b <- read_omics(gz, omics_type = "rnaseq", assay_type = "raw_count")
  expect_identical(colnames(b$input$expr_mat), c("ctrl_1", "ctrl_2", "ko_1", "ko_2"))
})

test_that("the encoding of compressed text is detected as for plain text", {
  dir <- withr::local_tempdir()
  genes <- c("TP53", "EGFR", "MYC", "KRAS", "BRCA1", "PTEN", "AKT1", "GAPDH",
             "ACTB", "VEGFA", "IL6", "TNF")
  set.seed(4)
  vals <- matrix(round(2^stats::rnorm(12 * 4, 18, 1)), 12)
  lines <- c("gene,Café_1,Café_2,Thé_1,Thé_2",
             apply(cbind(genes, vals), 1L, paste, collapse = ","))
  plain <- file.path(dir, "excel.csv")
  writeBin(iconv(list(charToRaw(enc2utf8(paste0(paste(lines, collapse = "\r\n"), "\r\n")))),
                 "UTF-8", "CP1252", toRaw = TRUE)[[1L]], plain)
  gz <- gzip_file(plain, file.path(dir, "excel.csv.gz"))
  b <- read_omics(gz, omics_type = "proteomics", assay_type = "raw_intensity")
  expect_identical(colnames(b$input$expr_mat),
                   enc2utf8(c("Café_1", "Café_2", "Thé_1", "Thé_2")))
  expect_true(any(grepl("read as Windows-1252", b$report$warnings, fixed = TRUE)))
})

test_that("a compressed sample sheet is read", {
  dir <- withr::local_tempdir()
  counts <- write_counts_csv(file.path(dir, "counts.csv"))
  sheet <- file.path(dir, "samples.csv")
  utils::write.csv(data.frame(sample_id = c("ctrl_1", "ctrl_2", "ko_1", "ko_2"),
                              group = c("ctrl", "ctrl", "ko", "ko")), sheet, row.names = FALSE)
  sheet_gz <- gzip_file(sheet, file.path(dir, "samples.csv.gz"))
  b <- read_omics(counts, omics_type = "rnaseq", assay_type = "raw_count",
                  sample_sheet = sheet_gz)
  expect_identical(b$input$meta_df$group, c("ctrl", "ctrl", "ko", "ko"))
})

test_that("compressed quantification files merge, named without the .gz", {
  dir <- withr::local_tempdir()
  ids <- sprintf("ENST%011d", 1:10)
  paths <- vapply(c("ctrl_1", "ko_1"), function(s) {
    plain <- file.path(dir, paste0(s, ".quant.sf"))
    utils::write.table(data.frame(Name = ids, Length = 1000, EffectiveLength = 850,
                                  TPM = 1e5, NumReads = 1:10 + 0.5),
                       plain, sep = "\t", quote = FALSE, row.names = FALSE)
    gzip_file(plain, paste0(plain, ".gz"))
  }, character(1))
  expect_identical(unname(quant_sample_names(paths)), c("ctrl_1", "ko_1"))
  expect_identical(quant_sample_names("salmon/ko_2/quant.sf.gz"), "ko_2")
  res <- read_quant_files(unname(paths))
  expect_identical(colnames(res$input$expr_mat), c("ctrl_1", "ko_1"))
  # And one of them through read_omics(), which recognises the format.
  one <- read_omics(paths[[1L]], omics_type = "rnaseq")
  expect_identical(one$report$suggested_input$quant_format, "salmon")
})

test_that("unpacking stops at the size limit, even when the file understates it", {
  dir <- withr::local_tempdir()
  plain <- file.path(dir, "big.csv")
  writeLines(c("gene,S1,S2", sprintf("G%07d,%d,%d", 1:120000, 1L, 2L)), plain)
  gz <- gzip_file(plain, file.path(dir, "big.csv.gz"))
  withr::local_options(omicsCore.max_unpacked_mb = 1)
  # The size recorded in the file is over the limit: refused up front.
  expect_error(read_omics(gz, omics_type = "rnaseq"), "unpacks to")
  # Two gzip members, the second tiny: the recorded size is the tiny
  # one's, and only counting the bytes as they come out catches it.
  tiny <- file.path(dir, "tiny.csv")
  writeLines("G9999999,1,2", tiny)
  two <- file.path(dir, "two.csv.gz")
  writeBin(c(readBin(gz, "raw", file.size(gz)),
             readBin(gzip_file(tiny, file.path(dir, "tiny.gz")), "raw", 100L)), two)
  expect_silent(guard_archive(two))
  err <- tryCatch(read_omics(two, omics_type = "rnaseq"), error = conditionMessage)
  expect_match(err, "This compressed file unpacks to more than the 1 MB", fixed = TRUE)
  expect_match(err, "raise options(omicsCore.max_unpacked_mb)", fixed = TRUE)
  expect_false(grepl("zlib|gzfile|connection", err))
  # No temporary copy is left behind.
  expect_length(list.files(tempdir(), pattern = "^omics-gunzip-"), 0L)
})

test_that("a damaged compressed file is refused in words", {
  dir <- withr::local_tempdir()
  plain <- file.path(dir, "counts.csv")
  writeLines(c("gene,S1,S2", sprintf("G%05d,%d,%d", 1:5000, 1L, 2L)), plain)
  gz <- gzip_file(plain, file.path(dir, "counts.csv.gz"))
  bytes <- readBin(gz, "raw", file.size(gz))
  cut <- file.path(dir, "cut.csv.gz")
  writeBin(bytes[seq_len(length(bytes) %/% 2L)], cut)
  # Cut short: the size it records is whatever bytes it stops on.
  err <- tryCatch(read_omics(cut, omics_type = "rnaseq"), error = conditionMessage)
  expect_match(err, "or it is damaged", fixed = TRUE)
  # Intact header, garbage body.
  bad <- file.path(dir, "bad.csv.gz")
  set.seed(2)
  writeBin(c(bytes[1:20], as.raw(sample(0:255, 500, TRUE)),
             as.raw(c(0, 0, 0, 0, 10, 0, 0, 0))), bad)
  err <- tryCatch(read_omics(bad, omics_type = "rnaseq"), error = conditionMessage)
  expect_match(err, "could not be unpacked; it may be damaged", fixed = TRUE)
})

test_that("a .gz that holds something other than text is refused by name", {
  dir <- withr::local_tempdir()
  f <- file.path(dir, "book.xlsx.gz")
  writeBin(as.raw(c(0x1f, 0x8b, 0, 0)), f)
  expect_error(detect_file_type(f), "holds '.xlsx'", fixed = TRUE)
})

test_that("a macro-enabled workbook (.xlsm) is read as a workbook", {
  skip_if_not_installed("openxlsx")
  dir <- withr::local_tempdir()
  xlsx <- file.path(dir, "book.xlsx")
  m <- matrix(round(2^stats::rnorm(6 * 4, 20, 1)), 6,
              dimnames = list(paste0("P", 1:6), paste0("S", 1:4)))
  openxlsx::write.xlsx(list(
    expression = data.frame(protein = rownames(m), m, check.names = FALSE),
    samples = data.frame(sample_id = colnames(m), group = c("A", "A", "B", "B"))), xlsx)
  xlsm <- file.path(dir, "book.xlsm")
  file.copy(xlsx, xlsm)
  expect_identical(detect_file_type(xlsm), "excel")
  a <- read_omics(xlsx, omics_type = "proteomics", assay_type = "raw_intensity")
  b <- read_omics(xlsm, omics_type = "proteomics", assay_type = "raw_intensity")
  expect_identical(b$input$expr_mat, a$input$expr_mat)
  expect_identical(b$input$meta_df$group, c("A", "A", "B", "B"))
  # As a sample sheet too.
  csv <- file.path(dir, "m.csv")
  utils::write.csv(data.frame(protein = rownames(m), m, check.names = FALSE), csv,
                   row.names = FALSE)
  s <- read_omics(csv, omics_type = "proteomics", assay_type = "raw_intensity",
                  sample_sheet = xlsm)
  expect_identical(s$input$meta_df$group, c("A", "A", "B", "B"))
  # Text calling itself .xlsm is read as the text it is.
  fake <- file.path(dir, "fake.xlsm")
  file.copy(csv, fake)
  expect_identical(detect_file_type(fake), "csv")
})
