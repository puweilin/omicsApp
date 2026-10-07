# SummarizedExperiment and DESeqDataSet objects, handed over in R or
# uploaded as .rds -- and, for an upload, read without running anything
# the file carries.

skip_if_not_installed("SummarizedExperiment")

se_counts <- function(n = 60, samples = 4, seed = 1) {
  set.seed(seed)
  matrix(stats::rnbinom(n * samples, mu = rep(c(10, 100, 1000), length.out = n), size = 5),
         n, dimnames = list(sprintf("ENSG%011d", seq_len(n)), sprintf("S%d", seq_len(samples))))
}

make_se <- function(assays, col = NULL, row = NULL, ...) {
  m <- assays[[1L]]
  if (is.null(col)) {
    col <- S4Vectors::DataFrame(group = rep(c("ctrl", "ko"), length.out = ncol(m)),
                                row.names = colnames(m))
  }
  SummarizedExperiment::SummarizedExperiment(assays, colData = col, rowData = row, ...)
}

save_rds <- function(x) {
  path <- tempfile(fileext = ".rds")
  saveRDS(x, path)
  path
}

test_that("a counts SummarizedExperiment becomes an RNA-seq layer", {
  m <- se_counts()
  m[1, ] <- c(5L, 7L, 9L, 11L)
  rownames(m)[1] <- "ENSG00000141510"
  se <- make_se(list(counts = m, tpm = m / 3),
                row = S4Vectors::DataFrame(gene_name = c("their_TP53", rep(NA, nrow(m) - 1L)),
                                           biotype = "protein_coding"))
  res <- read_summarized_experiment(se)
  inp <- res$input
  expect_s3_class(inp, "omics_input")
  expect_identical(inp$omics_type, "rnaseq")
  expect_identical(inp$assay_type, "raw_count")
  expect_identical(unname(inp$expr_mat), unname(m * 1))
  expect_identical(dimnames(inp$expr_mat), dimnames(m))
  expect_identical(inp$meta_df$group, c("ctrl", "ko", "ctrl", "ko"))
  expect_identical(rownames(inp$meta_df), colnames(m))
  expect_identical(inp$meta_df$sample_id, colnames(m))
  # The symbol column of rowData is the symbol, as in a workbook.
  expect_identical(inp$feature_df$feature_symbol[1], "their_TP53")
  expect_identical(inp$feature_df$biotype[1], "protein_coding")
  expect_true(any(grepl("holds 2 assays (counts, tpm); 'counts' was read as the matrix (the read counts)",
                        res$report$warnings, fixed = TRUE)))
  expect_identical(res$report$suggested_input$se_assay, "counts")
  expect_identical(res$report$sheets$role, c("matrix", "unknown", "metadata", "feature_annot"))
})

test_that("without a counts assay the first is read, and the choice is recorded", {
  set.seed(2)
  m <- matrix(stats::rnorm(40 * 4, 22, 2), 40,
              dimnames = list(paste0("P", 1:40), paste0("S", 1:4)))
  se <- make_se(list(log2_intensity = m, raw = 2^m))
  res <- read_summarized_experiment(se, omics_type = "proteomics")
  expect_identical(res$input$omics_type, "proteomics")
  expect_identical(res$input$assay_type, infer_assay_type(m, "proteomics"))
  expect_true(any(grepl("'log2_intensity' was read as the matrix (the first)",
                        res$report$warnings, fixed = TRUE)))
  # Another assay, by name or position.
  by_name <- read_summarized_experiment(se, omics_type = "proteomics", assay = "raw")
  by_pos <- read_summarized_experiment(se, omics_type = "proteomics", assay = 2)
  expect_equal(by_name$input$expr_mat, 2^m)
  expect_identical(by_pos$input$expr_mat, by_name$input$expr_mat)
  expect_true(any(grepl("(as asked)", by_name$report$warnings, fixed = TRUE)))
  expect_error(read_summarized_experiment(se, assay = "counts"),
               "`assay` must name one of the object's assays ('log2_intensity', 'raw'), not 'counts'.",
               fixed = TRUE)
  # With no modality, nothing is built and the report says so.
  none <- read_summarized_experiment(se)
  expect_null(none$input)
  expect_true(any(grepl("`omics_type` not specified", none$report$warnings, fixed = TRUE)))
})

test_that("read counts are RNA-seq whatever modality was asked for", {
  se <- make_se(list(counts = se_counts()))
  res <- read_summarized_experiment(se, omics_type = "proteomics", assay_type = "raw_intensity")
  expect_identical(res$input$omics_type, "rnaseq")
  expect_identical(res$input$assay_type, "raw_count")
  expect_true(any(grepl("read as RNA-seq, not proteomics", res$report$warnings, fixed = TRUE)))
})

test_that("an uploaded .rds of a SummarizedExperiment reads as the object does", {
  se <- make_se(list(counts = se_counts(), tpm = se_counts(seed = 9) / 7))
  from_r <- read_summarized_experiment(se)
  path <- save_rds(se)
  up <- read_omics(path, omics_type = "rnaseq")
  expect_identical(up$input$expr_mat, from_r$input$expr_mat)
  expect_identical(up$input$meta_df, from_r$input$meta_df)
  expect_identical(up$input$feature_df, from_r$input$feature_df)
  expect_identical(up$report$source, path)
})

test_that("a ranged object's feature names come from its rowRanges", {
  m <- se_counts(20)
  gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(seq_len(20) * 10, width = 5),
                               gene_name = paste0("g", 1:20))
  names(gr) <- rownames(m)
  rse <- SummarizedExperiment::SummarizedExperiment(list(counts = unname(m)), rowRanges = gr,
                                                    colData = S4Vectors::DataFrame(row.names = colnames(m)))
  res <- read_omics(save_rds(rse), omics_type = "rnaseq")
  expect_identical(rownames(res$input$expr_mat), rownames(m))
  expect_identical(res$input$feature_df$feature_symbol, paste0("g", 1:20))
  expect_identical(colnames(res$input$expr_mat), colnames(m))
})

test_that("tximeta's length assay travels as the effective lengths", {
  m <- se_counts(30) + 0.5
  len <- matrix(stats::runif(length(m), 500, 2000), nrow(m), dimnames = dimnames(m))
  se <- make_se(list(counts = m, abundance = m / 10, length = len))
  res <- read_omics(save_rds(se), omics_type = "rnaseq")
  expect_equal(res$input$misc$tximport$length, len)
  expect_identical(res$input$misc$tximport$counts_from_abundance, "no")
  expect_true(any(grepl("'length' assay was kept as the features' effective lengths",
                        res$report$warnings, fixed = TRUE)))
})

test_that("a sparse assay is spelled out", {
  skip_if_not_installed("Matrix")
  m <- se_counts(25)
  m[m < 50] <- 0L
  sp <- Matrix::Matrix(m * 1, sparse = TRUE)
  expect_s4_class(sp, "dgCMatrix")
  se <- make_se(list(counts = sp))
  res <- read_omics(save_rds(se), omics_type = "rnaseq")
  expect_identical(res$input$expr_mat, m * 1)
})

test_that("_PAR_Y rows and duplicated names are handled as in a table", {
  m <- se_counts(10)
  rownames(m)[1:3] <- c("ENSG00000182378.15", "ENSG00000182378.15_PAR_Y", "ENSG00000182378.15")
  m[2, ] <- 0L
  res <- read_summarized_experiment(make_se(list(counts = m)))
  ids <- rownames(res$input$expr_mat)
  expect_false("ENSG00000182378.15_PAR_Y" %in% ids)
  expect_false(anyDuplicated(ids) > 0L)
  expect_identical(nrow(res$input$feature_df), 9L)
  expect_true(any(grepl("repeated or empty and were made unique", res$report$warnings)))
})

test_that("a DESeqDataSet is read as its counts, without DESeq2's fit", {
  skip_if_not_installed("DESeq2")
  m <- se_counts(200)
  dds <- suppressMessages(DESeq2::DESeqDataSetFromMatrix(
    m, colData = data.frame(group = factor(c("a", "a", "b", "b")), row.names = colnames(m)),
    design = ~group))
  dds <- suppressMessages(DESeq2::DESeq(dds, fitType = "mean", quiet = TRUE))
  res <- read_omics(save_rds(dds), omics_type = "proteomics", assay_type = "raw_intensity")
  inp <- res$input
  expect_identical(inp$omics_type, "rnaseq")
  expect_identical(inp$assay_type, "raw_count")
  expect_identical(unname(inp$expr_mat), unname(m * 1))
  # Size factors and the fitted columns are DESeq2's, re-made when it runs.
  expect_false("sizeFactor" %in% names(inp$meta_df))
  expect_identical(as.character(inp$meta_df$group), c("a", "a", "b", "b"))
  expect_false(any(c("baseMean", "dispersion", "WaldPvalue_group_b_vs_a") %in%
                     names(inp$feature_df)))
  expect_true(any(grepl("DESeq2 had computed", res$report$warnings, fixed = TRUE)))
  expect_identical(res$report$suggested_input$se_class, "DESeqDataSet")
  # The same object handed over in R.
  expect_identical(read_summarized_experiment(dds)$input$expr_mat, inp$expr_mat)
})

# ---- an uploaded object is data, never code --------------------------------

test_that("an uploaded DESeqDataSet's functions and formula are never run", {
  skip_if_not_installed("DESeq2")
  m <- se_counts(50)
  dds <- suppressMessages(DESeq2::DESeqDataSetFromMatrix(
    m, colData = data.frame(group = factor(c("a", "a", "b", "b")), row.names = colnames(m)),
    design = ~group))
  flag <- tempfile()
  # A dispersion function that, if anything called it, leaves a trace,
  # and a design formula whose environment holds one too.
  trap <- function(...) {
    writeLines("ran", flag)
    stop("this must never run")
  }
  attr(dds, "dispersionFunction") <- trap
  env <- new.env()
  assign("trap", trap, envir = env)
  delayedAssign("bomb", trap(), assign.env = env)
  f <- ~group
  environment(f) <- env
  attr(dds, "design") <- f
  attr(dds, "metadata") <- list(hook = trap, env = env)
  path <- save_rds(dds)
  res <- read_omics(path, omics_type = "rnaseq")
  expect_false(file.exists(flag))
  expect_s3_class(res$input, "omics_input")
  # Nothing of it reached the layer.
  expect_silent(check_project_structure(res$input))
})

test_that("a function hidden in a sample column is left behind, and said", {
  se <- make_se(list(counts = se_counts()))
  cd <- attr(se, "colData")
  ld <- attr(cd, "listData")
  ld$hook <- function() stop("never")
  attr(cd, "listData") <- ld
  attr(se, "colData") <- cd
  res <- read_omics(save_rds(se), omics_type = "rnaseq")
  expect_s3_class(res$input, "omics_input")
  expect_false("hook" %in% names(res$input$meta_df))
  expect_true(any(grepl("not plain values: hook", res$report$warnings, fixed = TRUE)))
})

test_that("code hidden inside a column that looks plain is refused at the gate", {
  se <- make_se(list(counts = se_counts()))
  cd <- attr(se, "colData")
  ld <- attr(cd, "listData")
  g <- factor(ld$group)
  attr(g, "hook") <- function() stop("never")
  ld$group <- g
  attr(cd, "listData") <- ld
  attr(se, "colData") <- cd
  res <- read_omics(save_rds(se), omics_type = "rnaseq")
  expect_null(res$input)
  expect_match(res$report$warnings, "a function", fixed = TRUE)
  expect_match(res$report$warnings, "which a data table never holds", fixed = TRUE)
  # Levels that are not names: refused before anything builds a table
  # from them, in words.
  g2 <- structure(1:4, levels = list(function() stop("never")), class = "factor")
  ld$group <- g2
  attr(cd, "listData") <- ld
  attr(se, "colData") <- cd
  res <- read_omics(save_rds(se), omics_type = "rnaseq")
  expect_null(res$input)
  expect_match(res$report$warnings, "which a data table never holds", fixed = TRUE)
})

test_that("assays that are not matrices of numbers are refused in words", {
  se <- make_se(list(counts = se_counts()))
  a <- attr(se, "assays")
  d <- attr(a, "data")
  ld <- attr(d, "listData")
  ld$counts <- function() stop("never")
  attr(d, "listData") <- ld
  attr(a, "data") <- d
  attr(se, "assays") <- a
  res <- read_omics(save_rds(se), omics_type = "rnaseq")
  expect_null(res$input)
  expect_match(res$report$warnings, "could not be read", fixed = TRUE)
  expect_match(res$report$warnings, "Turn it into an ordinary matrix", fixed = TRUE)

  # The assays container of an old release is an environment.
  se2 <- make_se(list(counts = se_counts()))
  attr(se2, "assays") <- new.env()
  res2 <- read_omics(save_rds(se2), omics_type = "rnaseq")
  expect_null(res2$input)
  expect_match(res2$report$warnings, "updateObject()", fixed = TRUE)
})

test_that("an assay kept out of memory is read in R, and refused from a file", {
  skip_if_not_installed("DelayedArray")
  m <- se_counts(20)
  se <- make_se(list(counts = DelayedArray::DelayedArray(m)))
  # In R the object's own accessors are trusted.
  expect_identical(read_summarized_experiment(se)$input$expr_mat, m * 1)
  res <- read_omics(save_rds(se), omics_type = "rnaseq")
  expect_null(res$input)
  expect_match(res$report$warnings, "is stored as 'DelayedMatrix'", fixed = TRUE)
})

test_that("other S4 classes are still refused, a lookalike class name included", {
  df_path <- save_rds(S4Vectors::DataFrame(a = 1:3))
  res <- read_omics(df_path, omics_type = "rnaseq")
  expect_null(res$input)
  expect_match(res$report$warnings, "class 'DFrame'", fixed = TRUE)
  se <- make_se(list(counts = se_counts()))
  attr(se, "class") <- structure("MySummarizedExperiment", package = "elsewhere")
  res <- read_omics(save_rds(se), omics_type = "rnaseq")
  expect_null(res$input)
  expect_match(res$report$warnings, "class 'MySummarizedExperiment'", fixed = TRUE)
})

test_that("reading an uploaded object loads no package and runs no method", {
  skip_if_not_installed("callr")
  se <- make_se(list(counts = se_counts()))
  path <- save_rds(se)
  pkg <- testthat::test_path("..", "..")
  loaded <- callr::r(function(path, pkg) {
    pkgload::load_all(pkg, quiet = TRUE, export_all = FALSE)
    before <- loadedNamespaces()
    res <- omicsCore::read_omics(path, omics_type = "rnaseq")
    list(new = setdiff(loadedNamespaces(), before), ok = !is.null(res$input))
  }, args = list(path = path, pkg = pkg))
  expect_true(loaded$ok)
  expect_false(any(c("SummarizedExperiment", "S4Vectors", "DESeq2", "GenomicRanges") %in% loaded$new))
})

test_that("read_summarized_experiment() names a wrong argument", {
  expect_error(read_summarized_experiment(data.frame(a = 1)),
               "`x` must be a SummarizedExperiment or DESeqDataSet, not a data.frame with 1 row(s).",
               fixed = TRUE)
  se <- make_se(list(counts = se_counts()))
  expect_error(read_summarized_experiment(se, omics_type = "metabolomics"), "`omics_type`")
  expect_error(read_summarized_experiment(se, assay = list(1)), "`assay`")
})
