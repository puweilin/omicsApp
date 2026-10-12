# Import defects found in the 2026-10 audit.

skip_if_not_installed("writexl")

ai_genes <- c("TP53", "EGFR", "MYC", "KRAS", "BRCA1", "PTEN", "AKT1", "CDK4",
              "RB1", "APC", paste0("GENE", 1:20))

ai_write <- function(mat, meta, env = parent.frame()) {
  path <- withr::local_tempfile(fileext = ".xlsx", .local_envir = env)
  writexl::write_xlsx(list(expression = mat, samples = meta), path)
  path
}

ai_mat <- function(samples, seed = 2) {
  set.seed(seed)
  data.frame(Gene = ai_genes,
             matrix(round(stats::rnorm(30 * length(samples), 20, 1), 3), 30,
                    dimnames = list(NULL, samples)),
             check.names = FALSE)
}

test_that("metadata rows are put in the matrix's sample order", {
  path <- ai_write(ai_mat(paste0("S", 1:6)),
                   data.frame(sample = paste0("S", 6:1), group = rep(c("B", "A"), each = 3)))
  r <- read_omics(path, omics_type = "proteomics", assay_type = "normalized_intensity")
  expect_identical(rownames(r$input$meta_df), colnames(r$input$expr_mat))
  expect_identical(r$input$meta_df$group, rep(c("A", "B"), each = 3))
})

test_that("an extra or a missing metadata row is reported, not fatal", {
  extra <- ai_write(ai_mat(paste0("S", 1:6)),
                    data.frame(sample = paste0("S", 1:7), group = c(rep(c("A", "B"), each = 3), "B")))
  r <- read_omics(extra, omics_type = "proteomics", assay_type = "normalized_intensity")
  expect_false(is.null(r$input))
  expect_match(r$report$warnings, "not in the matrix", all = FALSE)

  missing <- ai_write(ai_mat(paste0("S", 1:6)),
                      data.frame(sample = paste0("S", 1:5), group = c("A", "A", "A", "B", "B")))
  r <- read_omics(missing, omics_type = "proteomics", assay_type = "normalized_intensity")
  expect_false(is.null(r$input))
  expect_true(is.na(r$input$meta_df["S6", "group"]))
  expect_match(r$report$warnings, "no metadata row", all = FALSE)
})

test_that("metadata is never assigned to samples by row position", {
  # Ids that differ only by "-" vs "_" are matched by name...
  path <- ai_write(ai_mat(c("Ctrl_1", "Ctrl_2", "Trt_1", "Trt_2")),
                   data.frame(sample = c("Trt-1", "Trt-2", "Ctrl-1", "Ctrl-2"),
                              group = c("Treated", "Treated", "Control", "Control")))
  r <- read_omics(path, omics_type = "proteomics", assay_type = "normalized_intensity")
  expect_identical(r$input$meta_df$group, c("Control", "Control", "Treated", "Treated"))
  expect_match(r$report$warnings, "ignoring case and punctuation", all = FALSE)

  # ...and ids that do not match at all attach no metadata, with a reason.
  path2 <- ai_write(ai_mat(c("A1", "A2", "B1", "B2")),
                    data.frame(sample = c("x", "y", "z", "w"), group = c("T", "T", "C", "C")))
  r2 <- read_omics(path2, omics_type = "proteomics", assay_type = "normalized_intensity")
  expect_false("group" %in% names(r2$input$meta_df))
  expect_match(r2$report$warnings, "No metadata column matches", all = FALSE)
})

test_that("numeric sample ids of 100000 and up still match", {
  ids <- c(100000, 100001, 100002, 100003)
  mat <- ai_mat(format(ids, scientific = FALSE, trim = TRUE))
  path <- ai_write(mat, data.frame(sample = ids, group = c("A", "A", "B", "B")))
  r <- read_omics(path, omics_type = "proteomics", assay_type = "normalized_intensity")
  expect_identical(r$input$meta_df$group, c("A", "A", "B", "B"))
})

test_that("duplicated feature ids get the matrix's de-duplicated label, and QC runs", {
  mat <- ai_mat(paste0("S", 1:6))
  mat$Gene[2] <- mat$Gene[1]
  annot <- data.frame(Gene = mat$Gene, description = paste("protein", seq_len(30)))
  path <- withr::local_tempfile(fileext = ".xlsx")
  writexl::write_xlsx(list(expression = mat,
                           samples = data.frame(sample = paste0("S", 1:6),
                                                group = rep(c("A", "B"), each = 3)),
                           features = annot), path)
  r <- read_omics(path, omics_type = "proteomics", assay_type = "normalized_intensity")
  inp <- r$input
  expect_identical(inp$feature_df$feature_id, rownames(inp$expr_mat))
  expect_false(anyDuplicated(inp$feature_df$feature_id) > 0L)
  expect_s3_class(run_qc(inp, outlier_method = "iqr", impute_method = "none"),
                  "analysis_bundle")
})

test_that("omics_input() reorders metadata, and validation refuses what it cannot fix", {
  m <- matrix(1:12 + 10, 3, dimnames = list(paste0("g", 1:3), paste0("s", 1:4)))
  meta <- data.frame(group = c("d", "c", "b", "a"), row.names = paste0("s", 4:1))
  feat <- data.frame(feature_id = paste0("g", 3:1), row.names = paste0("g", 3:1))
  x <- omics_input(m, meta, feat, omics_type = "proteomics",
                   assay_type = "normalized_intensity")
  expect_identical(x$meta_df$group, c("a", "b", "c", "d"))
  expect_identical(x$feature_df$feature_id, paste0("g", 1:3))

  bad <- x
  bad$meta_df <- bad$meta_df[4:1, , drop = FALSE]
  expect_error(validate_omics_input(bad), "not in the order")
  chr <- x
  chr$expr_mat <- matrix(as.character(m), 3, dimnames = dimnames(m))
  expect_error(validate_omics_input(chr), "must be numeric")
  dup <- x
  colnames(dup$expr_mat)[2] <- "s1"
  expect_error(validate_omics_input(dup), "duplicated sample names|not in the order")
  ff <- x
  ff$feature_df$feature_id <- c("zz", "yy", "xx")
  rownames(ff$feature_df) <- c("zz", "yy", "xx")
  expect_error(validate_omics_input(ff), "does not match the rows")
})

test_that("an Olink-style table with samples in rows is read the right way round", {
  set.seed(3)
  prot <- c("IL6", "TNF", "CXCL8", "IL10", "CCL2", "IL1B", "IFNG", "VEGFA",
            "MMP9", "CD40", "IL18", "CSF1", "FGF2", "HGF", "OSM", "TGFA",
            "IL17A", "CCL11", "CXCL10", "IL2")
  samp <- sprintf("S%02d", 1:24)
  npx <- data.frame(SampleID = samp,
                    matrix(round(stats::rnorm(24 * 20, 5), 3), 24,
                           dimnames = list(NULL, prot)),
                    check.names = FALSE)
  path <- withr::local_tempfile(fileext = ".xlsx")
  writexl::write_xlsx(list(npx = npx,
                           samples = data.frame(SampleID = samp,
                                                group = rep(c("A", "B"), 12))), path)
  r <- read_omics(path, omics_type = "proteomics", assay_type = "normalized_intensity")
  expect_identical(colnames(r$input$expr_mat), samp)
  expect_identical(rownames(r$input$expr_mat), prot)

  # And a caller can always say which way round it is.
  r2 <- read_omics(path, omics_type = "proteomics", assay_type = "normalized_intensity",
                   orientation = "samples_in_rows")
  expect_identical(r2$report$suggested_input$orientation_source, "user")
  expect_identical(r2$input$orientation, "samples_in_rows")
})

test_that("a guessed orientation is reported as a guess", {
  m <- ai_mat(c("A", "B", "C", "D", "E", "F"))   # names give nothing away
  m$Gene <- paste0("feat", letters[seq_len(30) %% 26 + 1], seq_len(30) %/% 26)
  m <- m[1:8, ]
  path <- withr::local_tempfile(fileext = ".xlsx")
  writexl::write_xlsx(list(x = m), path)
  r <- read_omics(path, omics_type = "proteomics", assay_type = "normalized_intensity")
  if (r$report$suggested_input$orientation_confidence < 0.6) {
    expect_match(r$report$warnings, "was a guess", all = FALSE)
  } else {
    succeed()
  }
})

test_that("an all-numeric sample sheet is not taken for the matrix", {
  ids <- 101:108
  mat <- ai_mat(as.character(ids))
  meta <- data.frame(id = ids, age = c(30, 41, 52, 33, 44, 55, 36, 47),
                     bmi = round(stats::runif(8, 18, 30), 1),
                     arm = rep(0:1, 4))
  path <- withr::local_tempfile(fileext = ".xlsx")
  writexl::write_xlsx(list(Sheet1 = mat, Sheet2 = meta), path)
  r <- read_omics(path, omics_type = "proteomics", assay_type = "normalized_intensity")
  expect_identical(r$report$suggested_input$matrix_sheet, "Sheet1")
  expect_identical(r$report$suggested_input$metadata_sheet, "Sheet2")
  expect_identical(dim(r$input$expr_mat), c(30L, 8L))
  expect_identical(r$input$meta_df$age, meta$age)
  expect_match(r$report$warnings, "name the samples", all = FALSE)
})

test_that("a workbook that unpacks past the limit is refused before it is read", {
  path <- withr::local_tempfile(fileext = ".xlsx")
  writexl::write_xlsx(list(expression = ai_mat(paste0("S", 1:6))), path)
  withr::local_options(omicsCore.max_unpacked_mb = 0.001)
  expect_error(read_omics(path, omics_type = "proteomics"), "unpacks to")
})

test_that("a saved omics_input that is not valid is reported, not loaded", {
  inp <- omics_input(matrix(1:12 + 10, 3, dimnames = list(paste0("g", 1:3), paste0("s", 1:4))),
                     data.frame(group = c("a", "a", "b", "b"), row.names = paste0("s", 1:4)),
                     data.frame(feature_id = paste0("g", 1:3)),
                     omics_type = "proteomics", assay_type = "normalized_intensity")
  inp$meta_df <- inp$meta_df[4:1, , drop = FALSE]
  path <- withr::local_tempfile(fileext = ".rds")
  saveRDS(inp, path)
  r <- read_omics(path)
  expect_null(r$input)
  expect_match(r$report$warnings, "not valid", all = FALSE)
})
