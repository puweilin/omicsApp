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

# ---- extension vs content --------------------------------------------
# RNA-seq pipelines routinely write tab-separated text and name it .xls.
# Trusting the extension sends it to readxl, which fails with a message
# about the workbook rather than about the format.

test_that("a tab-separated file named .xls is read as text, not Excel", {
  path <- withr::local_tempfile(fileext = ".xls")
  writeLines(c("gene_id\tS1\tS2", "G1\t1\t2", "G2\t3\t4"), path)

  expect_false(is_excel_file(path))
  expect_identical(detect_file_type(path), "csv")
  expect_identical(detect_delimiter(path), "\t")
})

test_that("a real xlsx is still detected as Excel", {
  skip_if_not_installed("writexl")
  path <- withr::local_tempfile(fileext = ".xlsx")
  writexl::write_xlsx(data.frame(gene_id = "G1", S1 = 1), path)

  expect_true(is_excel_file(path))
  expect_identical(detect_file_type(path), "excel")
})

test_that("detect_delimiter picks the delimiter the header actually uses", {
  tsv <- withr::local_tempfile(); writeLines("a\tb\tc", tsv)
  csv <- withr::local_tempfile(); writeLines("a,b,c", csv)
  ssv <- withr::local_tempfile(); writeLines("a;b;c", ssv)
  one <- withr::local_tempfile(); writeLines("only_one_column", one)

  expect_identical(detect_delimiter(tsv), "\t")
  expect_identical(detect_delimiter(csv), ",")
  expect_identical(detect_delimiter(ssv), ";")
  expect_identical(detect_delimiter(one), ",")   # nothing to go on
})

test_that("apostrophes and # in annotation columns do not eat rows", {
  # read.table's defaults treat ' as an opening quote and # as a comment.
  # A gene called 5'-nucleotidase opens a quote that stays open until the
  # next apostrophe, and every row in between is swallowed -- surfacing
  # thousands of lines later as `line N did not have K elements`.
  #
  # Asserted on the parse, not on the classified matrix: how the sheet is
  # then labelled depends on the shape of the fixture, and this is about
  # whether every row survived reading.
  path <- withr::local_tempfile(fileext = ".xls")
  samples <- paste0("S", 1:6)
  descs <- rep(c("5'-nucleotidase", 'some "quoted" name', "contains # hash",
                 "apostrophe's here", "plain"), 6)
  rows <- vapply(seq_along(descs), function(i) {
    paste(c(paste0("G", i), descs[i], as.character(seq_along(samples) * i)),
          collapse = "\t")
  }, character(1L))
  writeLines(c(paste(c("gene_id", "description", samples), collapse = "\t"), rows),
             path)

  res <- read_omics(path, omics_type = "rnaseq", assay_type = "fpkm")
  expect_equal(res$report$sheets$n_rows[[1L]], length(descs))
  expect_equal(res$report$sheets$n_cols[[1L]], 2L + length(samples))
})

test_that("default read.table quoting would have lost those rows", {
  # The counter-test: without the fix the same file parses short or
  # errors outright, which is what makes the one above worth keeping.
  path <- withr::local_tempfile(fileext = ".tsv")
  writeLines(c("gene_id\tdescription\tS1",
               "G1\t5'-nucleotidase\t1",
               "G2\tplain\t2",
               "G3\tanother'\t3"), path)

  lenient <- utils::read.table(path, header = TRUE, sep = "\t",
                               quote = "", comment.char = "",
                               check.names = FALSE, stringsAsFactors = FALSE)
  default <- tryCatch(
    utils::read.table(path, header = TRUE, sep = "\t",
                      check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) NULL)

  expect_equal(nrow(lenient), 3L)
  expect_true(is.null(default) || nrow(default) < 3L)
})

# ---- one measurement per sample --------------------------------------
# Vendor reports carry every sample twice, once per unit. Read whole,
# half the columns are on a different scale from the other half and
# nothing downstream can tell -- PCA, clustering and differential
# testing all run and all mean nothing.

mat_named <- function(nms) {
  m <- matrix(seq_len(4 * length(nms)), nrow = 4,
              dimnames = list(paste0("G", 1:4), nms))
  m
}

test_that("counts win when a table carries both units", {
  res <- select_measurement_columns(
    mat_named(c("A_FPKM", "B_FPKM", "A_count", "B_count")))

  expect_identical(colnames(res$mat), c("A", "B"))
  expect_match(res$note, "kept the 2 count")
  expect_match(res$note, "dropped 2")
})

test_that("the suffix is dropped so both report layouts name samples alike", {
  # The counts-only file from the same vendor calls the sample
  # `RD001_Folli`; metadata written against one must match the other.
  res <- select_measurement_columns(mat_named(c("RD001_Folli_count",
                                                "RD002_Folli_count")))
  expect_identical(colnames(res$mat), c("RD001_Folli", "RD002_Folli"))
})

test_that("a single non-count unit is kept rather than discarded", {
  res <- select_measurement_columns(mat_named(c("A_TPM", "B_TPM")))
  expect_identical(colnames(res$mat), c("A", "B"))
})

test_that("columns without a recognised unit are left alone", {
  nms <- c("RD001_Folli", "RD002_Folli")
  res <- select_measurement_columns(mat_named(nms))
  expect_identical(colnames(res$mat), nms)
  expect_null(res$note)
})

test_that("a partly-suffixed table is not touched", {
  # Acting on a partial match would drop real samples whose names happen
  # to end in a word that looks like a unit.
  nms <- c("A_count", "B_count", "SomeSample")
  res <- select_measurement_columns(mat_named(nms))
  expect_identical(colnames(res$mat), nms)
  expect_null(res$note)
})

test_that("the unit suffix is read as the last word, not across separators", {
  # A greedy match reads RD001_Folli_FPKM as the unit `Folli_FPKM`,
  # which matches nothing -- the table then looks unit-free and both
  # units survive into one matrix.
  expect_identical(measure_suffix("RD001_Folli_FPKM"), "fpkm")
  expect_identical(measure_suffix("RD001_Folli_count"), "count")
  expect_true(is.na(measure_suffix("RD001_Folli")))
})

test_that("names are left intact when stripping would collide", {
  res <- select_measurement_columns(mat_named(c("A_count", "A_counts")))
  expect_identical(colnames(res$mat), c("A_count", "A_counts"))
})
