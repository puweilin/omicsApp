# Delimited text is read with data.table::fread() when it is installed,
# and with read.table otherwise. The fast path is only allowed to be
# faster: on every file here the two must hand back the same data frame,
# down to column types, missing values and names. Where fread would
# disagree with read.table, the reader is expected to notice and defer.

# Written as UTF-8 bytes whatever the locale of the test process, as a
# user's file would be.
fixture <- function(lines, ext = ".csv", bom = FALSE) {
  path <- tempfile(fileext = ext)
  con <- file(path, "wb")
  if (bom) writeBin(as.raw(c(0xef, 0xbb, 0xbf)), con)
  writeBin(charToRaw(paste0(paste(enc2utf8(lines), collapse = "\n"), "\n")), con)
  close(con)
  path
}

read_with <- function(path, fast) {
  withr::with_options(list(omicsCore.use_fread = fast),
                      tryCatch(read_delimited_table(path),
                               error = function(e) paste("error:", conditionMessage(e))))
}

# TRUE when read_delimited_fast() produced the table, FALSE when it
# deferred to read.table.
took_fast_path <- function(path) {
  real <- read_delimited_fast
  took <- NA
  testthat::local_mocked_bindings(read_delimited_fast = function(...) {
    out <- real(...)
    took <<- !is.null(out)
    out
  })
  withr::with_options(list(omicsCore.use_fread = TRUE), read_delimited_table(path))
  took
}

expect_same_both_ways <- function(path) {
  fast <- read_with(path, TRUE)
  slow <- read_with(path, FALSE)
  expect_identical(fast, slow)
  invisible(slow)
}

counts_lines <- function(sep = ",", n = 25L) {
  set.seed(3)
  ids <- sprintf("ENSG%011d", seq_len(n))
  m <- matrix(stats::rpois(n * 4L, 150), n)
  c(paste(c("gene_id", sprintf("S%02d", 1:4)), collapse = sep),
    apply(cbind(ids, m), 1L, paste, collapse = sep))
}

# ---- the fixtures the fast path must take, and agree on ---------------

test_that("featureCounts: the command line is skipped the same way", {
  skip_if_not_installed("data.table")
  path <- fixture(c(
    "# Program:featureCounts v2.0.1; Command:\"featureCounts\" -a genes.gtf -o counts.txt",
    "Geneid\tChr\tStart\tEnd\tStrand\tLength\tS1.bam\tS2.bam",
    "G1\tchr1;chr1\t11;50\t20;90\t+;+\t51\t5\t6",
    "G2\tchr2\t100\t200\t-\t101\t0\t1",
    "G3\tchr3\t5\t9\t+\t5\t12\t30"), ext = ".txt")
  expect_true(took_fast_path(path))
  out <- expect_same_both_ways(path)
  expect_identical(colnames(out$df)[1:2], c("Geneid", "Chr"))
  expect_type(out$df$Length, "integer")
  expect_identical(out$df$Start, c("11;50", "100", "5"))
})

test_that("decimal commas under ';' are read as decimals, the same way", {
  skip_if_not_installed("data.table")
  path <- fixture(c("protein;S1;S2;S3",
                    "P1;1,5;12,345;0,001",
                    "P2;NA;3,25;4",
                    "P3;7;;1e3"))
  expect_true(took_fast_path(path))
  out <- expect_same_both_ways(path)
  expect_equal(out$df$S2, c(12.345, 3.25, NA))
  expect_equal(out$df$S1, c(1.5, NA, 7))
})

test_that("an apostrophe in a tab-separated file is a letter, not a quote", {
  skip_if_not_installed("data.table")
  path <- fixture(c("gene_id\tdescription\tS1\tS2",
                    "G1\t5'-nucleotidase, ecto\t1\t2",
                    "G2\t3' repair exonuclease\t3\t4",
                    "G3\t\"quoted\" in a TSV\t5\t6"), ext = ".tsv")
  expect_true(took_fast_path(path))
  out <- expect_same_both_ways(path)
  expect_identical(out$df$description,
                   c("5'-nucleotidase, ecto", "3' repair exonuclease",
                     "\"quoted\" in a TSV"))
})

test_that("a quoted header name holding the delimiter is one column", {
  skip_if_not_installed("data.table")
  path <- fixture(c("gene,\"Sample, A\",\"Sample, B\",S3",
                    sprintf("G%d,%d.5,%d.25,%d", 1:30, 1:30, 31:60, 61:90)))
  expect_true(took_fast_path(path))
  out <- expect_same_both_ways(path)
  expect_identical(colnames(out$df), c("gene", "Sample, A", "Sample, B", "S3"))
  expect_length(out$duplicated_headers, 0L)
  # It used to stop with "line 1 did not have 6 elements"; now it is a
  # matrix whose sample names are the ones in the file.
  for (fast in c(TRUE, FALSE)) {
    res <- withr::with_options(list(omicsCore.use_fread = fast),
      read_omics(path, omics_type = "proteomics", assay_type = "normalized_intensity"))
    expect_identical(colnames(res$input$expr_mat), c("Sample, A", "Sample, B", "S3"))
  }
})

test_that("missing values and spreadsheet errors come out as read.table has them", {
  skip_if_not_installed("data.table")
  # fread reads "#N/A" as NA and "#DIV/0!" as NaN; read.table keeps them
  # as text, which is what the matrix cleanup counts and reports.
  path <- fixture(c("id,a,b,c,d,e,f",
                    "x,#DIV/0!,1,\"NA\",T,,0x1A",
                    "y,2.50,#N/A,4,F,NA,2",
                    "NA,NaN,3,5,TRUE,1,3",
                    "\"NA\",Inf,4,6,FALSE,2,4"))
  out <- expect_same_both_ways(path)
  expect_identical(out$df$a, c("#DIV/0!", "2.50", "NaN", "Inf"))
  expect_identical(out$df$b, c("1", "#N/A", "3", "4"))
  expect_identical(out$df$c, c(NA, 4L, 5L, 6L))
  expect_identical(out$df$d, c(TRUE, FALSE, TRUE, FALSE))
  expect_identical(out$df$e, c(NA, NA, 1L, 2L))
  expect_identical(out$df$f, c(26, 2, 3, 4))
  expect_identical(out$df$id, c("x", "y", NA, NA))
})

test_that("gaps in every column do not cost the fast path its answer", {
  skip_if_not_installed("data.table")
  set.seed(5)
  m <- matrix(round(stats::rlnorm(40 * 5, 18, 2), 3), 40)
  m[sample(length(m), 30)] <- NA
  lines <- c("protein,S1,S2,S3,S4,S5",
             apply(cbind(sprintf("P%03d", 1:40), ifelse(is.na(m), "", m)), 1L,
                   paste, collapse = ","))
  path <- fixture(lines)
  expect_true(took_fast_path(path))
  expect_same_both_ways(path)
  # The same gaps written as "NA", and a "#" in an annotation column that
  # makes the reader check the raw text of every gap.
  lines2 <- c("protein,note,S1,S2,S3,S4,S5",
              apply(cbind(sprintf("P%03d", 1:40), "isoform #2",
                          ifelse(is.na(m), "NA", m)), 1L, paste, collapse = ","))
  path2 <- fixture(lines2)
  expect_true(took_fast_path(path2))
  expect_same_both_ways(path2)
})

test_that("a byte-order mark and non-ASCII names read the same way", {
  skip_if_not_installed("data.table")
  path <- fixture(c("gene,样本1,échantillon",
                    "广,1,2", "G2,3,4"), bom = TRUE)
  expect_true(took_fast_path(path))
  out <- expect_same_both_ways(path)
  expect_identical(colnames(out$df), c("gene", "样本1", "échantillon"))
  expect_identical(Encoding(colnames(out$df))[2:3], c("UTF-8", "UTF-8"))
})

test_that("an empty or repeated header name is kept as read.table keeps it", {
  skip_if_not_installed("data.table")
  path <- fixture(c("gene,S01,S01,,S03", "G1,1,2,3,4", "G2,5,6,7,8"))
  expect_true(took_fast_path(path))
  out <- expect_same_both_ways(path)
  expect_identical(colnames(out$df), c("gene", "S01", "S01", "", "S03"))
  expect_identical(out$duplicated_headers, "S01")
})

# ---- the files it must hand to read.table -----------------------------

test_that("a ragged row fails exactly as it did, instead of being dropped", {
  skip_if_not_installed("data.table")
  # fread stops early with a warning and keeps the rows before; that
  # would be a matrix silently missing its tail.
  lines <- counts_lines()
  lines[12] <- paste0(lines[12], ",99")
  path <- fixture(lines)
  fast <- read_with(path, TRUE)
  expect_match(fast, "did not have", fixed = TRUE)
  expect_identical(fast, read_with(path, FALSE))

  short <- counts_lines()
  short[3] <- "ENSG00000000002,1"
  path <- fixture(short)
  expect_identical(read_with(path, TRUE), read_with(path, FALSE))
})

test_that("a header one name short (row names first) defers to read.table", {
  skip_if_not_installed("data.table")
  path <- fixture(c("S1,S2", "G1,1,2", "G2,3,4"))
  expect_false(took_fast_path(path))
  expect_same_both_ways(path)
})

test_that("a quote inside an unquoted field defers to read.table", {
  skip_if_not_installed("data.table")
  path <- fixture(c("id,a,note", "x,1,he said \"hi\"", "y,2,plain"))
  expect_false(took_fast_path(path))
  expect_same_both_ways(path)
})

test_that("dates are text, as read.table leaves them", {
  skip_if_not_installed("data.table")
  path <- fixture(c("id,when,stamp,S1", "x,2020-01-01,2020-01-01T10:00:00,1",
                    "y,2021-02-03,2021-02-03T11:30:00,2"))
  out <- expect_same_both_ways(path)
  expect_type(out$df$when, "character")
  expect_type(out$df$stamp, "character")
})

test_that("read.table arguments passed through `...` are honoured by read.table", {
  skip_if_not_installed("data.table")
  path <- fixture(c("id,S1,S2", "x,1,-", "y,2,3"))
  out <- read_delimited_table(path, na.strings = c("NA", "-"))
  expect_identical(out$df$S2, c(NA, 3L))
})

test_that("without data.table the reader is read.table, unchanged", {
  path <- fixture(counts_lines())
  testthat::local_mocked_bindings(
    is_installed = function(pkg) pkg != "data.table" && requireNamespace(pkg, quietly = TRUE)
  )
  expect_false(took_fast_path(path))
  expect_same_both_ways(path)
})

# ---- end to end -------------------------------------------------------

test_that("read_omics builds the same input either way, warnings included", {
  skip_if_not_installed("data.table")
  lines <- counts_lines(sep = "\t")
  lines[1] <- "gene_id\tS01\tS02\tS02\tS04"
  path <- fixture(lines, ext = ".tsv")
  fast <- withr::with_options(list(omicsCore.use_fread = TRUE),
    read_omics(path, omics_type = "rnaseq", assay_type = "raw_count"))
  slow <- withr::with_options(list(omicsCore.use_fread = FALSE),
    read_omics(path, omics_type = "rnaseq", assay_type = "raw_count"))
  expect_identical(fast$input, slow$input)
  expect_identical(fast$report$warnings, slow$report$warnings)
  expect_true(any(grepl("'S02'", fast$report$warnings, fixed = TRUE)))
  expect_false(anyDuplicated(colnames(fast$input$expr_mat)) > 0L)
})
