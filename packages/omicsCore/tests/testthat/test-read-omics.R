# Tests for read_omics() + ImportReport + classifier (slice 1G).
#
# We synthesise tiny workbooks/CSVs/RDS files in tempdir(); never write
# beyond `tempdir()` so the test suite stays self-contained.

make_demo_matrix <- function(n_features = 12L, n_samples = 6L) {
  set.seed(2030)
  expr <- matrix(
    rnorm(n_features * n_samples, mean = 5, sd = 1.0),
    nrow = n_features,
    dimnames = list(paste0("GENE", seq_len(n_features)),
                    paste0("S", seq_len(n_samples)))
  )
  expr
}

write_demo_workbook <- function(path, expr,
                                include_meta = TRUE,
                                include_feat = TRUE,
                                meta_name = "metadata",
                                feat_name = "feature_annot") {
  expr_df <- data.frame(
    feature_id = rownames(expr),
    expr,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  sheets <- list(expression = expr_df)
  if (include_meta) {
    meta_df <- data.frame(
      sample_id = colnames(expr),
      group = rep(c("ctrl", "case"), length.out = ncol(expr)),
      age = seq.int(30L, 30L + ncol(expr) - 1L),
      stringsAsFactors = FALSE
    )
    sheets[[meta_name]] <- meta_df
  }
  if (include_feat) {
    feat_df <- data.frame(
      feature_id = rownames(expr),
      feature_symbol = rownames(expr),
      description = paste("desc for", rownames(expr)),
      stringsAsFactors = FALSE
    )
    sheets[[feat_name]] <- feat_df
  }
  openxlsx::write.xlsx(sheets, file = path, overwrite = TRUE)
  path
}

# ---- ImportReport class -----------------------------------------------

test_that("new_import_report builds the expected shape", {
  r <- new_import_report()
  expect_s3_class(r, "ImportReport")
  expect_true(is_import_report(r))
  expect_s3_class(r$sheets, "data.frame")
  expect_equal(nrow(r$sheets), 0L)
  expect_identical(r$warnings, character(0))
  expect_identical(r$suggested_input, list())
})

test_that("new_import_report coerces partial sheets to the canonical schema", {
  sheets <- data.frame(name = "expression", role = "matrix",
                       stringsAsFactors = FALSE)
  r <- new_import_report(sheets = sheets)
  expect_setequal(colnames(r$sheets),
                  c("name", "role", "n_rows", "n_cols",
                    "confidence", "orientation", "notes"))
  expect_type(r$sheets$confidence, "double")
})

test_that("print.ImportReport works without erroring", {
  r <- new_import_report(
    sheets = data.frame(
      name = "expression", role = "matrix", n_rows = 10L, n_cols = 5L,
      confidence = 0.9, orientation = "features_in_rows",
      notes = "ok", stringsAsFactors = FALSE
    ),
    warnings = "test warning",
    suggested_input = list(matrix_sheet = "expression",
                           omics_type = "proteomics")
  )
  expect_output(print(r), "ImportReport")
  expect_output(print(r), "expression")
  expect_output(print(r), "test warning")
})

# ---- classify_sheet_role ----------------------------------------------

test_that("classify_sheet_role detects a numeric matrix", {
  expr <- make_demo_matrix()
  df <- data.frame(feature_id = rownames(expr), expr,
                   check.names = FALSE, stringsAsFactors = FALSE)
  cls <- classify_sheet_role(df, name = "expression")
  expect_identical(cls$role, "matrix")
  expect_gte(cls$confidence, 0.7)
})

test_that("classify_sheet_role detects metadata by column hints", {
  meta_df <- data.frame(
    sample_id = paste0("S", 1:6),
    group = rep(c("ctrl", "case"), 3L),
    age = c(30, 40, 50, 60, 70, 80),
    stringsAsFactors = FALSE
  )
  cls <- classify_sheet_role(meta_df, name = "metadata")
  expect_identical(cls$role, "metadata")
})

test_that("classify_sheet_role detects a feature annotation table", {
  feat_df <- data.frame(
    feature_id = paste0("GENE", 1:5),
    feature_symbol = paste0("GENE", 1:5),
    description = paste("desc", 1:5),
    stringsAsFactors = FALSE
  )
  cls <- classify_sheet_role(feat_df, name = "feature_annot")
  expect_identical(cls$role, "feature_annot")
})

test_that("classify_sheet_role returns unknown for empty input", {
  cls <- classify_sheet_role(data.frame())
  expect_identical(cls$role, "unknown")
  expect_equal(cls$confidence, 0)
})

# ---- detect_orientation -----------------------------------------------

test_that("detect_orientation favours feature-like row labels", {
  expr <- make_demo_matrix()
  df <- data.frame(feature_id = rownames(expr), expr,
                   check.names = FALSE, stringsAsFactors = FALSE)
  res <- detect_orientation(df)
  expect_identical(res$orientation, "features_in_rows")
})

# ---- detect_id_columns ------------------------------------------------

test_that("detect_id_columns flags UniProt-shaped IDs", {
  df <- data.frame(
    uniprot = c("P12345", "Q67890", "O11111", "P22222", "Q33333"),
    desc = letters[1:5],
    stringsAsFactors = FALSE
  )
  hits <- detect_id_columns(df)
  expect_true("uniprot" %in% hits$column)
  expect_true(hits$pattern[hits$column == "uniprot"] %in%
              c("uniprot", "hgnc_symbol"))
})

test_that("detect_id_columns flags Ensembl-shaped IDs", {
  df <- data.frame(
    gene = paste0("ENSG0000000", 1000:1006),
    stringsAsFactors = FALSE
  )
  hits <- detect_id_columns(df)
  expect_true("gene" %in% hits$column)
  expect_true("ensembl_gene" %in% hits$pattern)
})

# ---- read_omics: Excel path -------------------------------------------

test_that("read_omics reads an Excel workbook and returns input + report", {
  expr <- make_demo_matrix()
  path <- tempfile(fileext = ".xlsx")
  on.exit(unlink(path), add = TRUE)
  write_demo_workbook(path, expr)

  out <- read_omics(path, omics_type = "proteomics",
                    assay_type = "normalized_intensity")
  expect_named(out, c("input", "report"))
  expect_s3_class(out$report, "ImportReport")
  expect_s3_class(out$input, "omics_input")
  expect_equal(nrow(out$input$expr_mat), nrow(expr))
  expect_equal(ncol(out$input$expr_mat), ncol(expr))
  expect_identical(out$input$omics_type, "proteomics")
})

test_that("read_omics ImportReport carries one row per sheet", {
  expr <- make_demo_matrix()
  path <- tempfile(fileext = ".xlsx")
  on.exit(unlink(path), add = TRUE)
  write_demo_workbook(path, expr)

  out <- read_omics(path, omics_type = "proteomics")
  expect_equal(nrow(out$report$sheets), 3L)
  expect_setequal(out$report$sheets$role,
                  c("matrix", "metadata", "feature_annot"))
  expect_true(all(out$report$sheets$confidence >= 0.4))
  expect_identical(out$report$suggested_input$matrix_sheet, "expression")
})

test_that("read_omics returns input = NULL when omics_type is missing", {
  expr <- make_demo_matrix()
  path <- tempfile(fileext = ".xlsx")
  on.exit(unlink(path), add = TRUE)
  write_demo_workbook(path, expr)

  out <- read_omics(path)
  expect_null(out$input)
  expect_true(any(grepl("omics_type", out$report$warnings)))
})

test_that("read_omics still produces a report with no metadata sheet", {
  expr <- make_demo_matrix()
  path <- tempfile(fileext = ".xlsx")
  on.exit(unlink(path), add = TRUE)
  write_demo_workbook(path, expr, include_meta = FALSE)

  out <- read_omics(path, omics_type = "proteomics")
  expect_s3_class(out$input, "omics_input")
  # Synthetic sample_id metadata column should fill in.
  expect_true("sample_id" %in% colnames(out$input$meta_df))
})

# ---- read_omics: CSV path ---------------------------------------------

test_that("read_omics reads a TSV matrix file", {
  expr <- make_demo_matrix()
  df <- data.frame(feature_id = rownames(expr), expr,
                   check.names = FALSE, stringsAsFactors = FALSE)
  path <- tempfile(fileext = ".tsv")
  on.exit(unlink(path), add = TRUE)
  utils::write.table(df, file = path, sep = "\t",
                     quote = FALSE, row.names = FALSE)

  out <- read_omics(path, omics_type = "rnaseq")
  expect_s3_class(out$input, "omics_input")
  expect_equal(nrow(out$report$sheets), 1L)
  expect_identical(out$report$sheets$role, "matrix")
})

# ---- read_omics: RDS path ---------------------------------------------

test_that("read_omics passes an omics_input through unchanged", {
  expr <- make_demo_matrix()
  meta <- data.frame(group = rep(c("ctrl", "case"), 3L),
                     row.names = colnames(expr), stringsAsFactors = FALSE)
  feat <- data.frame(feature_id = rownames(expr),
                     row.names = rownames(expr), stringsAsFactors = FALSE)
  input <- omics_input(expr, meta, feat, omics_type = "proteomics")
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(input, path)

  out <- read_omics(path)
  expect_s3_class(out$input, "omics_input")
  expect_identical(out$input$omics_type, "proteomics")
  expect_equal(out$report$sheets$confidence, 1.0)
})

test_that("read_omics RDS path tolerates a bare data.frame", {
  expr <- make_demo_matrix()
  df <- data.frame(feature_id = rownames(expr), expr,
                   check.names = FALSE, stringsAsFactors = FALSE)
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(df, path)

  out <- read_omics(path, omics_type = "proteomics")
  expect_s3_class(out$report, "ImportReport")
  expect_s3_class(out$input, "omics_input")
})

# ---- read_omics: failure modes ----------------------------------------

test_that("read_omics rejects nonexistent paths", {
  expect_error(read_omics(tempfile(fileext = ".xlsx")),
               "does not exist")
})

test_that("read_omics rejects unsupported extensions", {
  path <- tempfile(fileext = ".bin")
  file.create(path)
  on.exit(unlink(path), add = TRUE)
  expect_error(read_omics(path), "auto-detect")
})

test_that("read_omics flags an empty workbook gracefully", {
  path <- tempfile(fileext = ".xlsx")
  on.exit(unlink(path), add = TRUE)
  openxlsx::write.xlsx(list(empty = data.frame()), file = path, overwrite = TRUE)

  out <- read_omics(path, omics_type = "proteomics")
  expect_null(out$input)
  expect_true(length(out$report$warnings) > 0L)
})
