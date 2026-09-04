# Files as the instruments and pipelines actually write them.
#
# Every import bug this project has hit had the same shape: a real file
# that the heuristics turned into a plausible, wrong matrix, or refused
# outright. A .xls that was tab-separated text. A gene called
# 5'-nucleotidase that opened a quote. Both units of a vendor report
# side by side in one matrix. And, found while writing this file: a
# MaxQuant proteinGroups.txt with six samples is 70% annotation and QC
# columns, so it was "no sheet looked like an expression matrix"; a
# Spectronaut report that says "Filtered" in one cell lost that whole
# sample; a featureCounts table started with its command line and
# stopped at "more columns than column names".
#
# The generators below write each shape into a temp file at test time.
# Nothing binary is committed, and every case asserts the same
# contract: the samples and only the samples become columns, the ids
# become rows, the symbol column is found, the numbers survive, and
# whatever the reader decided on the user's behalf is in the report.

set.seed(1907)
N_FEAT <- 40L
N_SAMP <- 6L
SAMPLES <- sprintf("S%02d", seq_len(N_SAMP))
PROT_IDS <- sprintf("P%05d", seq_len(N_FEAT))
GENE_IDS <- sprintf("ENSG%011d", seq_len(N_FEAT))
GENES <- c("TP53;TP53BP1", "MYC", "EGFR", "CDK1", "CCNB1", "PLK1", "AURKA",
           "BUB1", "MKI67", "TOP2A", "E2F1", "MCM2", sprintf("GENE%02d", 13:N_FEAT))
LFQ <- matrix(round(2^stats::rnorm(N_FEAT * N_SAMP, 22, 1.5)), N_FEAT, N_SAMP,
              dimnames = list(NULL, SAMPLES))
COUNTS <- matrix(stats::rpois(N_FEAT * N_SAMP, 300), N_FEAT, N_SAMP,
                 dimnames = list(NULL, SAMPLES))

# Written by hand rather than with write.table, so the bytes are exactly
# what the vendor writes: no quoting, chosen line endings, optional BOM
# and leading lines.
write_delim <- function(df, path, sep = "\t", eol = "\n", bom = FALSE,
                        leading = character(0)) {
  con <- file(path, "wb")
  on.exit(close(con))
  if (bom) writeBin(as.raw(c(0xef, 0xbb, 0xbf)), con)
  for (line in leading) writeLines(line, con, sep = eol)
  writeLines(paste(colnames(df), collapse = sep), con, sep = eol)
  for (i in seq_len(nrow(df))) {
    writeLines(paste(vapply(df[i, ], as.character, character(1)), collapse = sep),
               con, sep = eol)
  }
  invisible(path)
}

vendor_file <- function(df, ext, ...) {
  path <- tempfile(fileext = ext)
  write_delim(df, path, ...)
  path
}

maxquant_table <- function() {
  df <- data.frame(
    `Protein IDs` = PROT_IDS, `Majority protein IDs` = PROT_IDS,
    `Protein names` = paste("protein", seq_len(N_FEAT)), `Gene names` = GENES,
    `Fasta headers` = paste0(">sp|", PROT_IDS),
    Peptides = sample(2:30, N_FEAT, TRUE),
    `Razor + unique peptides` = sample(1:20, N_FEAT, TRUE),
    `Unique peptides` = sample(1:20, N_FEAT, TRUE),
    `Sequence coverage [%]` = round(stats::runif(N_FEAT, 5, 80), 1),
    `Mol. weight [kDa]` = round(stats::runif(N_FEAT, 10, 200), 2),
    Score = round(stats::runif(N_FEAT, 5, 300), 2),
    Intensity = rowSums(LFQ),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  for (s in SAMPLES) df[[paste("Intensity", s)]] <- LFQ[, s] * 2
  for (s in SAMPLES) df[[paste("LFQ intensity", s)]] <- LFQ[, s]
  df[["MS/MS count"]] <- sample(1:100, N_FEAT, TRUE)
  df[["Only identified by site"]] <- ""
  df$Reverse <- ""
  df[["Potential contaminant"]] <- c("+", rep("", N_FEAT - 1L))
  df
}

spectronaut_table <- function() {
  df <- data.frame(
    PG.ProteinGroups = PROT_IDS, PG.Genes = sub(";.*", "", GENES),
    PG.ProteinDescriptions = paste("desc", seq_len(N_FEAT)),
    PG.Organisms = "Homo sapiens",
    check.names = FALSE, stringsAsFactors = FALSE
  )
  for (j in seq_len(N_SAMP)) {
    df[[sprintf("[%d] %s.raw.PG.Quantity", j, SAMPLES[j])]] <- as.character(LFQ[, j])
  }
  df[[5L]][2L] <- "Filtered"
  df
}

diann_table <- function() {
  df <- data.frame(
    Protein.Group = PROT_IDS, Protein.Ids = PROT_IDS,
    Protein.Names = paste0(sub(";.*", "", GENES), "_HUMAN"),
    Genes = sub(";.*", "", GENES),
    First.Protein.Description = paste("desc", seq_len(N_FEAT)),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  for (j in seq_len(N_SAMP)) df[[sprintf("D:\\data\\run1\\%s.raw", SAMPLES[j])]] <- LFQ[, j]
  df
}

featurecounts_table <- function() {
  df <- data.frame(
    Geneid = GENE_IDS, Chr = "chr1;chr1", Start = "100;900", End = "500;1300",
    Strand = "+;+", Length = sample(500:5000, N_FEAT),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  counts <- COUNTS
  colnames(counts) <- paste0(SAMPLES, ".bam")
  cbind(df, counts)
}

dual_unit_table <- function() {
  df <- data.frame(gene_id = GENE_IDS, gene_name = sub(";.*", "", GENES),
                   check.names = FALSE, stringsAsFactors = FALSE)
  for (s in SAMPLES) {
    df[[paste0(s, "_count")]] <- COUNTS[, s]
    df[[paste0(s, "_FPKM")]] <- round(COUNTS[, s] / 7.3, 2)
  }
  df
}

plain_counts_table <- function() {
  data.frame(gene_id = GENE_IDS, COUNTS, check.names = FALSE,
             stringsAsFactors = FALSE)
}

read_prot <- function(path) {
  read_omics(path, omics_type = "proteomics", assay_type = "raw_intensity")
}
read_counts <- function(path) {
  read_omics(path, omics_type = "rnaseq", assay_type = "raw_count")
}

expect_sample_block <- function(out, rows, values = NULL) {
  m <- out$input$expr_mat
  expect_false(is.null(m), info = paste(out$report$warnings, collapse = " | "))
  expect_identical(dim(m), c(N_FEAT, N_SAMP))
  expect_identical(colnames(m), SAMPLES)
  expect_identical(rownames(m), rows)
  if (!is.null(values)) {
    expect_equal(unname(m[!is.na(m)]), unname(values[!is.na(m)]))
  }
}

expect_symbols_from <- function(out, expected_first) {
  expect_identical(utils::head(out$input$feature_df$feature_symbol, 4L),
                   expected_first)
}

# ---- proteomics search engines ---------------------------------------------

test_that("MaxQuant proteinGroups.txt: the LFQ block is the matrix, nothing else is", {
  out <- read_prot(vendor_file(maxquant_table(), ".txt"))
  expect_sample_block(out, PROT_IDS, LFQ)
  expect_symbols_from(out, c("TP53", "MYC", "EGFR", "CDK1"))
  expect_true(any(grepl("LFQ intensity", out$report$warnings, fixed = TRUE)))
  expect_true(any(grepl("dropped", out$report$warnings)))
})

test_that("Spectronaut report: 'Filtered' is a missing value, and the sample keeps its name", {
  out <- read_prot(vendor_file(spectronaut_table(), ".tsv"))
  expect_sample_block(out, PROT_IDS)
  m <- out$input$expr_mat
  expect_true(is.na(m[2L, "S01"]))
  expect_identical(sum(is.na(m)), 1L)
  expect_equal(unname(m[1L, "S01"]), unname(LFQ[1L, "S01"]))
  expect_symbols_from(out, c("TP53", "MYC", "EGFR", "CDK1"))
  expect_true(any(grepl("PG.Quantity", out$report$warnings, fixed = TRUE)))
})

test_that("DIA-NN pg_matrix: raw-file paths become sample names", {
  out <- read_prot(vendor_file(diann_table(), ".tsv"))
  expect_sample_block(out, PROT_IDS, LFQ)
  expect_symbols_from(out, c("TP53", "MYC", "EGFR", "CDK1"))
  expect_true(any(grepl("S01.raw", out$report$warnings, fixed = TRUE)))
})

# ---- sequencing pipelines ----------------------------------------------------

test_that("featureCounts: the command line is skipped and the annotation columns are not samples", {
  path <- vendor_file(featurecounts_table(), ".txt",
                      leading = "# Program:featureCounts v2.0.1; Command:featureCounts -a genes.gtf -o counts.txt")
  out <- read_counts(path)
  expect_sample_block(out, GENE_IDS, COUNTS)
  expect_true(any(grepl("featureCounts", out$report$warnings, fixed = TRUE)))
  expect_true(all(out$input$expr_mat == round(out$input$expr_mat)))
})

test_that("a dual-unit report named .xls is read as text and keeps the counts", {
  out <- read_counts(vendor_file(dual_unit_table(), ".xls"))
  expect_sample_block(out, GENE_IDS, COUNTS)
  expect_symbols_from(out, c("TP53", "MYC", "EGFR", "CDK1"))
  expect_true(any(grepl("count", out$report$warnings)))
})

test_that("a merged quantification table with versioned Ensembl ids maps to symbols", {
  skip_if(is.null(hgnc_ensembl_map()), "bundled HGNC table not available")
  ids <- GENE_IDS
  ids[1:3] <- c("ENSG00000141510.17", "ENSG00000136997.21", "ENSG00000146648.19")
  tpm <- round(COUNTS / 3.7, 3)
  df <- data.frame(Name = ids, tpm, check.names = FALSE, stringsAsFactors = FALSE)
  out <- read_omics(vendor_file(df, ".tsv"), omics_type = "rnaseq", assay_type = "tpm")
  m <- out$input$expr_mat
  expect_identical(dim(m), c(N_FEAT, N_SAMP))
  expect_equal(unname(m[1L, ]), unname(tpm[1L, ]))
  expect_identical(out$input$feature_df$feature_symbol[1:3], c("TP53", "MYC", "EGFR"))
})

# ---- the same table, every way it can be delivered ---------------------------

test_that("csv, tsv, CRLF, BOM, .xlsx and TSV-named-.xls all give the identical matrix", {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("readxl")
  base <- plain_counts_table()
  paths <- list(
    csv  = vendor_file(base, ".csv", sep = ","),
    tsv  = vendor_file(base, ".tsv"),
    crlf = vendor_file(base, ".tsv", eol = "\r\n"),
    bom  = vendor_file(base, ".tsv", bom = TRUE),
    xls  = vendor_file(base, ".xls"),
    semi = vendor_file(base, ".csv", sep = ";")
  )
  xlsx <- tempfile(fileext = ".xlsx")
  openxlsx::write.xlsx(base, xlsx)
  paths$xlsx <- xlsx
  mats <- lapply(paths, function(p) read_counts(p)$input$expr_mat)
  for (nm in names(mats)) {
    expect_equal(mats[[nm]], mats$csv, info = nm)
    expect_identical(colnames(mats[[nm]]), SAMPLES, info = nm)
    expect_identical(rownames(mats[[nm]]), GENE_IDS, info = nm)
  }
})

# ---- numbers written for people ------------------------------------------------

test_that("every spelling of 'missing' is NA, and the column stays a sample", {
  df <- plain_counts_table()
  df$S01 <- as.character(df$S01)
  tokens <- c("NA", "", "NaN", "#N/A", "Filtered", "n.d.")
  df$S01[seq_along(tokens)] <- tokens
  out <- read_counts(vendor_file(df, ".csv", sep = ","))
  m <- out$input$expr_mat
  expect_identical(colnames(m), SAMPLES)
  expect_true(all(is.na(m[seq_along(tokens), "S01"])))
  expect_equal(unname(m[-seq_along(tokens), "S01"]),
               unname(COUNTS[-seq_along(tokens), "S01"]))
  expect_identical(sum(is.na(m)), length(tokens))
})

test_that("thousands separators are read as the numbers they group", {
  df <- plain_counts_table()
  df$S02 <- format(df$S02 * 1000L, big.mark = ",", trim = TRUE)
  out <- read_counts(vendor_file(df, ".tsv"))
  m <- out$input$expr_mat
  expect_identical(colnames(m), SAMPLES)
  expect_equal(unname(m[, "S02"]), unname(COUNTS[, "S02"]) * 1000)
})

test_that("a numeric column stored as text in a workbook is still numeric", {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("readxl")
  df <- plain_counts_table()
  df$S03 <- as.character(df$S03)
  xlsx <- tempfile(fileext = ".xlsx")
  openxlsx::write.xlsx(df, xlsx)
  out <- read_counts(xlsx)
  expect_equal(unname(out$input$expr_mat[, "S03"]), unname(COUNTS[, "S03"]))
})

# ---- layout and identifiers ------------------------------------------------------

test_that("samples in rows are recognised and transposed", {
  df <- data.frame(sample = SAMPLES, t(LFQ), check.names = FALSE,
                   stringsAsFactors = FALSE)
  colnames(df)[-1L] <- PROT_IDS
  out <- read_prot(vendor_file(df, ".csv", sep = ","))
  expect_sample_block(out, PROT_IDS, LFQ)
  expect_identical(out$report$suggested_input$orientation, "samples_in_rows")
})

test_that("duplicate and empty ids are made unique rather than dropped", {
  df <- plain_counts_table()
  df$gene_id[2L] <- df$gene_id[1L]
  df$gene_id[3L] <- ""
  out <- read_counts(vendor_file(df, ".csv", sep = ","))
  m <- out$input$expr_mat
  expect_identical(nrow(m), N_FEAT)
  expect_false(anyDuplicated(rownames(m)) > 0L)
  expect_identical(rownames(m)[1:2], c(GENE_IDS[1L], paste0(GENE_IDS[1L], "_1")))
  expect_true(rownames(m)[3L] != "")
})

test_that("a symbol column mangled by Excel is still picked up, as written", {
  # SEPT2 and MARCH1 come back from a spreadsheet as dates. Nothing here
  # can undo that; what it must not do is fail, or lose the column.
  df <- data.frame(gene_id = GENE_IDS,
                   gene_symbol = c("2-Sep", "1-Mar", sub(";.*", "", GENES)[3:N_FEAT]),
                   COUNTS, check.names = FALSE, stringsAsFactors = FALSE)
  out <- read_counts(vendor_file(df, ".csv", sep = ","))
  expect_sample_block(out, GENE_IDS, COUNTS)
  expect_identical(out$input$feature_df$feature_symbol[1:3], c("2-Sep", "1-Mar", "EGFR"))
})

test_that("a plain matrix with a few annotation columns is still a matrix", {
  # Six samples beside four text columns is 60% numeric -- under the
  # old 80% rule this was not a matrix at all.
  df <- data.frame(gene_id = GENE_IDS, gene_name = sub(";.*", "", GENES),
                   description = paste("desc", seq_len(N_FEAT)),
                   chromosome = "chr1", COUNTS,
                   check.names = FALSE, stringsAsFactors = FALSE)
  out <- read_counts(vendor_file(df, ".tsv"))
  expect_sample_block(out, GENE_IDS, COUNTS)
  expect_symbols_from(out, c("TP53", "MYC", "EGFR", "CDK1"))
})

test_that("a sample sheet is not mistaken for a matrix", {
  # The relaxed matrix rule must not catch metadata: sample ids in the
  # first column, a couple of numeric covariates, and a group column.
  meta <- data.frame(sample_id = sprintf("S%03d", 1:30),
                     group = rep(c("G1", "G2"), 15),
                     age = round(stats::runif(30, 20, 70)),
                     bmi = round(stats::runif(30, 18, 32), 1),
                     batch = rep(c("A", "B", "C"), 10),
                     stringsAsFactors = FALSE)
  cls <- classify_sheet_role(meta, name = "samples")
  expect_identical(cls$role, "metadata")
})
