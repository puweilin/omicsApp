# Salmon, RSEM and kallisto write one file per sample. Read as a matrix,
# a quant.sf became four "samples" (Length, EffectiveLength, TPM,
# NumReads). They are recognised by their columns, and read_quant_files()
# merges several into one counts layer with the effective lengths kept
# where DESeq2 and edgeR look for tximport's.

write_quant <- function(path, format, ids, reads, efflen, tpm = NULL, gene_id = NULL) {
  if (is.null(tpm)) {
    tpm <- reads / efflen
    tpm <- tpm / sum(tpm) * 1e6
  }
  df <- switch(format,
    salmon = data.frame(Name = ids, Length = round(efflen) + 150, EffectiveLength = efflen,
                        TPM = tpm, NumReads = reads),
    rsem_genes = data.frame(gene_id = ids, `transcript_id(s)` = paste0(ids, "-T"),
                            length = round(efflen) + 150, effective_length = efflen,
                            expected_count = reads, TPM = tpm, FPKM = tpm / 2,
                            check.names = FALSE),
    rsem_isoforms = data.frame(transcript_id = ids, gene_id = gene_id,
                               length = round(efflen) + 150, effective_length = efflen,
                               expected_count = reads, TPM = tpm, FPKM = tpm / 2,
                               IsoPct = 100),
    kallisto = data.frame(target_id = ids, length = round(efflen) + 150, eff_length = efflen,
                          est_counts = reads, tpm = tpm)
  )
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  utils::write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)
  path
}

test_that("each format is recognised by its columns", {
  expect_identical(detect_quant_format(c("Name", "Length", "EffectiveLength", "TPM", "NumReads")),
                   "salmon")
  expect_identical(detect_quant_format(c("gene_id", "transcript_id(s)", "length",
                                         "effective_length", "expected_count", "TPM", "FPKM")),
                   "rsem")
  expect_identical(detect_quant_format(c("transcript_id", "gene_id", "length", "effective_length",
                                         "expected_count", "TPM", "FPKM", "IsoPct")), "rsem")
  expect_identical(detect_quant_format(c("target_id", "length", "eff_length", "est_counts", "tpm")),
                   "kallisto")
  expect_identical(detect_quant_format(c("gene_id", "S1", "S2", "TPM")), NA_character_)
})

test_that("sample names come from the file, or its directory", {
  expect_identical(
    quant_sample_names(c("/x/ctrl_1/quant.sf", "ko_2.quant.sf", "a/wt_3.genes.results",
                         "wt_4.isoforms.results", "out/ko_5/abundance.tsv", "s6_abundance.tsv",
                         "C:\\runs\\s7\\quant.sf", "quant.sf")),
    c("ctrl_1", "ko_2", "wt_3", "wt_4", "ko_5", "s6", "s7", ""))
})

test_that("Salmon files merge into one counts layer, lengths and names kept", {
  dir <- withr::local_tempdir()
  ids <- sprintf("ENST%011d.1", 1:6)
  reads <- list(c(10, 0, 5.5, 100, 3, 8), c(12, 1, 7.25, 90, 0, 9), c(30, 2, 1, 50, 4, 10))
  efflen <- list(c(100, 200, 300, 400, 500, 600), c(110, 210, 310, 410, 510, 610),
                 c(120, 220, 320, 420, 520, 620))
  samples <- c("ctrl_1", "ctrl_2", "ko_1")
  paths <- vapply(1:3, function(i) write_quant(file.path(dir, samples[i], "quant.sf"),
                                               "salmon", ids, reads[[i]], efflen[[i]]),
                  character(1))
  res <- read_quant_files(paths)
  inp <- res$input
  expect_s3_class(inp, "omics_input")
  expect_identical(inp$omics_type, "rnaseq")
  expect_identical(inp$assay_type, "raw_count")
  expect_identical(colnames(inp$expr_mat), samples)
  expect_identical(rownames(inp$expr_mat), ids)
  expect_equal(unname(inp$expr_mat), do.call(cbind, reads))
  len <- inp$misc$tximport$length
  expect_equal(unname(len), do.call(cbind, efflen))
  expect_identical(dimnames(len), dimnames(inp$expr_mat))
  expect_identical(inp$misc$tximport$counts_from_abundance, "no")
  # Exactly where the engines look.
  expect_equal(get_tximport_info(inp)$length, len)
  expect_true(any(grepl("Read 3 Salmon quantification file", res$report$warnings)))
  # And it survives a QC exclusion.
  expect_equal(subset_omics_samples(inp, samples[1:2])$misc, inp$misc)
})

test_that("RSEM and kallisto files merge the same way", {
  dir <- withr::local_tempdir()
  ids <- sprintf("ENSG%011d", 1:5)
  r1 <- c(1, 2.5, 3, 4, 50)
  r2 <- c(2, 3.5, 0, 8, 40)
  e1 <- c(500, 600, 700, 800, 900)
  e2 <- c(510, 610, 0, 810, 910)
  rsem <- c(write_quant(file.path(dir, "wt_1.genes.results"), "rsem_genes", ids, r1, e1),
            write_quant(file.path(dir, "wt_2.genes.results"), "rsem_genes", ids, r2, e2))
  res <- read_quant_files(rsem)
  expect_identical(colnames(res$input$expr_mat), c("wt_1", "wt_2"))
  expect_equal(unname(res$input$expr_mat), cbind(r1, r2), ignore_attr = TRUE)
  # RSEM writes a zero effective length for a feature shorter than the
  # fragments; it becomes 1, as tximport's importer does, so the offsets
  # stay finite.
  expect_equal(unname(res$input$misc$tximport$length), cbind(e1, c(510, 610, 1, 810, 910)),
               ignore_attr = TRUE)

  tx <- sprintf("ENST%011d", 1:5)
  iso <- write_quant(file.path(dir, "wt_1.isoforms.results"), "rsem_isoforms", tx, r1, e1,
                     gene_id = ids)
  expect_identical(read_quant_files(iso)$input$feature_df$gene_id, ids)

  kal <- c(write_quant(file.path(dir, "k1", "abundance.tsv"), "kallisto", tx, r1, e1),
           write_quant(file.path(dir, "k2", "abundance.tsv"), "kallisto", tx, r2, e1))
  res <- read_quant_files(kal, sample_names = c("A_1", "B_1"))
  expect_identical(colnames(res$input$expr_mat), c("A_1", "B_1"))
  expect_equal(unname(res$input$expr_mat[, "B_1"]), r2)
})

test_that("files that do not name their samples are numbered, and told", {
  dir <- withr::local_tempdir()
  ids <- paste0("T", 1:4)
  # Two files whose names and directories are the same.
  a <- write_quant(file.path(dir, "a", "salmon", "quant.sf"), "salmon", ids, 1:4, rep(100, 4))
  b <- write_quant(file.path(dir, "b", "salmon", "quant.sf"), "salmon", ids, 4:1, rep(100, 4))
  res <- read_quant_files(c(a, b))
  expect_identical(colnames(res$input$expr_mat), c("sample_1", "sample_2"))
  expect_true(any(grepl("do not tell the samples apart", res$report$warnings)))
  # An upload's temporary copies are named by the names they came with.
  res <- read_quant_files(c(a, b), file_names = c("ctrl_1.quant.sf", "ko_1.quant.sf"))
  expect_identical(colnames(res$input$expr_mat), c("ctrl_1", "ko_1"))
  expect_identical(res$report$source, "ctrl_1.quant.sf, ko_1.quant.sf")
  expect_error(read_quant_files(c(a, b), file_names = "x.sf"), "`file_names` must name each")
  expect_error(read_quant_files(c(a, b), sample_names = c("s", "s")), "unique")
  expect_error(read_quant_files(a, sample_names = c("s", "t")), "each of the 1 file")
})

test_that("tx2gene sums transcripts to genes as tximport does", {
  dir <- withr::local_tempdir()
  tx <- c("ENST01.1", "ENST02.2", "ENST03.1", "ENST04.1", "ENST05.3", "ENST06.1")
  # Transcript 6 is not in tx2gene.
  t2g <- data.frame(tx = c("ENST01", "ENST02", "ENST03", "ENST04", "ENST05"),
                    gene = c("G1", "G1", "G2", "G3", "G4"))
  s1 <- write_quant(file.path(dir, "s1.quant.sf"), "salmon", tx,
                    reads = c(10, 30, 5, 0, 0, 9), efflen = c(100, 300, 500, 50, 80, 10),
                    tpm = c(100, 300, 600, 0, 0, 1))
  s2 <- write_quant(file.path(dir, "s2.quant.sf"), "salmon", tx,
                    reads = c(20, 0, 7, 0, 4, 9), efflen = c(200, 400, 500, 70, 90, 10),
                    tpm = c(300, 0, 700, 0, 10, 1))
  res <- read_quant_files(c(s1, s2), tx2gene = t2g)
  inp <- res$input
  expect_identical(rownames(inp$expr_mat), c("G1", "G2", "G3", "G4"))
  # Counts summed.
  expect_equal(unname(inp$expr_mat), rbind(c(40, 20), c(5, 7), c(0, 0), c(0, 4)))
  len <- inp$misc$tximport$length
  # G1: TPM-weighted mean of 100 and 300 (weights 100, 300) = 250 in s1;
  # in s2 transcript 2 has no TPM, so it is transcript 1's 200.
  expect_equal(unname(len["G1", ]), c((100 * 100 + 300 * 300) / 400, 200))
  expect_equal(unname(len["G2", ]), c(500, 500))
  # G3 has no TPM anywhere: the mean of its transcript's lengths over
  # samples, (50 + 70) / 2.
  expect_equal(unname(len["G3", ]), c(60, 60))
  # G4 has none in s1: the geometric mean of its length elsewhere.
  expect_equal(unname(len["G4", ]), c(90, 90))
  expect_true(any(grepl("ignoring their version suffixes", res$report$warnings)))
  expect_true(any(grepl("1 of 6 transcripts are not in tx2gene", res$report$warnings)))
  expect_true(any(grepl("Summed 5 transcripts to 4 genes", res$report$warnings)))
})

test_that("files that cannot be merged are refused in words", {
  dir <- withr::local_tempdir()
  a <- write_quant(file.path(dir, "a.quant.sf"), "salmon", paste0("T", 1:4), 1:4, rep(100, 4))
  b <- write_quant(file.path(dir, "b.quant.sf"), "salmon", paste0("T", 1:5), 1:5, rep(100, 5))
  k <- write_quant(file.path(dir, "k", "abundance.tsv"), "kallisto", paste0("T", 1:4), 1:4,
                   rep(100, 4))
  m <- file.path(dir, "matrix.tsv")
  writeLines(c("gene\tS1\tS2", "G1\t1\t2"), m)
  expect_error(read_quant_files(c(a, b)), "different references")
  expect_error(read_quant_files(c(a, k)), "different quantifiers \\(Salmon, kallisto\\)")
  expect_error(read_quant_files(c(a, m)), "not a Salmon, RSEM or kallisto quantification file")
  expect_error(read_quant_files(a, tx2gene = data.frame(tx = "X", gene = "Y")),
               "None of the transcript ids")
  expect_error(read_quant_files(character(0)), "`paths`")
})

test_that("read_omics() reads one quantification file as one sample, and says why it is not enough", {
  dir <- withr::local_tempdir()
  ids <- sprintf("ENST%011d", 1:8)
  path <- write_quant(file.path(dir, "ctrl_1", "quant.sf"), "salmon", ids, 1:8 + 0.5, rep(1000, 8))
  res <- read_omics(path, omics_type = "rnaseq")
  expect_identical(dim(res$input$expr_mat), c(8L, 1L))
  expect_identical(colnames(res$input$expr_mat), "ctrl_1")
  expect_false(any(c("Length", "EffectiveLength", "TPM", "NumReads") %in%
                     colnames(res$input$expr_mat)))
  expect_true(any(grepl("one sample cannot be compared", res$report$warnings)))
  expect_identical(res$report$suggested_input$quant_format, "salmon")
  # kallisto's abundance.tsv the same way, whatever omics type was asked.
  kal <- write_quant(file.path(dir, "ko_1", "abundance.tsv"), "kallisto", ids, 1:8, rep(900, 8))
  res <- read_omics(kal, omics_type = "proteomics")
  expect_identical(res$input$omics_type, "rnaseq")
  expect_true(any(grepl("read as RNA-seq, not proteomics", res$report$warnings)))
})

test_that("a merged layer runs through edgeR and DESeq2 with its length offsets", {
  skip_if_not_installed("edgeR")
  dir <- withr::local_tempdir()
  set.seed(5)
  n <- 400
  ids <- sprintf("ENSG%011d", seq_len(n))
  samples <- c(paste0("ctrl_", 1:3), paste0("ko_", 1:3))
  base_len <- stats::runif(n, 800, 3000)
  paths <- vapply(samples, function(s) {
    ko <- startsWith(s, "ko")
    len <- base_len
    # The first 40 genes switch to an isoform twice as long in ko: the
    # same number of molecules gives twice the reads. Nothing changed in
    # expression, and only the length offsets can tell.
    if (ko) len[1:40] <- len[1:40] * 2
    mu <- 200 * len / 1000
    reads <- stats::rnbinom(n, mu = mu, size = 50) + 0.25
    write_quant(file.path(dir, s, "quant.sf"), "salmon", ids, reads, len)
  }, character(1))
  res <- read_quant_files(paths, sample_sheet = {
    sheet <- file.path(dir, "samples.csv")
    utils::write.csv(data.frame(sample_id = samples, group = rep(c("ctrl", "ko"), each = 3)),
                     sheet, row.names = FALSE)
    sheet
  })
  inp <- res$input
  expect_identical(inp$meta_df$group, rep(c("ctrl", "ko"), each = 3))
  run <- function(x, method) {
    run_diff(x, method = method, group_col = "group", control_group = "ctrl",
             case_group = "ko")$results$diff_result_df
  }
  with_len <- run(inp, "edger")
  no_len <- inp
  no_len$misc <- NULL
  without <- run(no_len, "edger")
  switched <- with_len$feature_id %in% ids[1:40]
  # Without offsets the longer isoform reads as a 2-fold increase; with
  # them, as no change.
  expect_gt(stats::median(without$effect[without$feature_id %in% ids[1:40]]), 0.8)
  expect_lt(abs(stats::median(with_len$effect[switched])), 0.25)
  skip_if_not_installed("DESeq2")
  # DESeq2 takes the fractional counts through its tximport route, which
  # it refuses without the lengths.
  ds <- suppressMessages(run(inp, "deseq2"))
  expect_lt(abs(stats::median(ds$effect[ds$feature_id %in% ids[1:40]])), 0.25)
})
