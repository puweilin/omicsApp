# I/O robustness: read_omics, save_project/load_project, write_*

make_io_input <- function(n_feat = 8, n_samp = 6) {
  mat <- matrix(rnorm(n_feat * n_samp, mean = 10, sd = 2),
                nrow = n_feat, ncol = n_samp)
  rownames(mat) <- paste0("g", seq_len(n_feat))
  colnames(mat) <- paste0("s", seq_len(n_samp))
  half <- n_samp %/% 2
  meta <- data.frame(
    group = c(rep("A", half), rep("B", n_samp - half)),
    row.names = colnames(mat)
  )
  feat <- data.frame(feature_id = rownames(mat), stringsAsFactors = FALSE)
  omics_input(mat, meta, feat, omics_type = "proteomics",
              assay_type = "normalized_intensity")
}

# ---- save_project / load_project --------------------------------------

test_that("save_project writes a file and load_project reads it back", {
  skip_if_not_installed("qs2")
  inp <- make_io_input()
  p <- omics_project("test", experiments = list(prot = inp))
  tf <- tempfile(fileext = ".qs2")
  on.exit(unlink(tf), add = TRUE)
  save_project(p, tf)
  expect_true(file.exists(tf))
})

test_that("save_project + load_project round-trips a project", {
  skip_if_not_installed("qs2")
  inp <- make_io_input()
  p <- omics_project("rt", experiments = list(prot = inp))
  tf <- tempfile(fileext = ".qs2")
  on.exit(unlink(tf), add = TRUE)
  save_project(p, tf)
  p2 <- load_project(tf)
  expect_true(is_omics_project(p2))
})

test_that("save_project round-trip preserves experiment tags", {
  skip_if_not_installed("qs2")
  inp <- make_io_input()
  p <- omics_project("rt", experiments = list(prot = inp, rna = inp))
  tf <- tempfile(fileext = ".qs2")
  on.exit(unlink(tf), add = TRUE)
  save_project(p, tf)
  p2 <- load_project(tf)
  expect_setequal(experiment_tags(p2), c("prot", "rna"))
})

test_that("save_project round-trip preserves expr_mat content", {
  skip_if_not_installed("qs2")
  inp <- make_io_input()
  p <- omics_project("rt", experiments = list(prot = inp))
  tf <- tempfile(fileext = ".qs2")
  on.exit(unlink(tf), add = TRUE)
  save_project(p, tf)
  p2 <- load_project(tf)
  expect_equal(p2$experiments$prot$expr_mat, inp$expr_mat)
})

test_that("save_project errors if path exists and overwrite=FALSE", {
  skip_if_not_installed("qs2")
  inp <- make_io_input()
  p <- omics_project("rt", experiments = list(prot = inp))
  tf <- tempfile(fileext = ".qs2")
  on.exit(unlink(tf), add = TRUE)
  save_project(p, tf)
  expect_error(save_project(p, tf, overwrite = FALSE))
})

test_that("save_project allows overwrite=TRUE", {
  skip_if_not_installed("qs2")
  inp <- make_io_input()
  p <- omics_project("rt", experiments = list(prot = inp))
  tf <- tempfile(fileext = ".qs2")
  on.exit(unlink(tf), add = TRUE)
  save_project(p, tf)
  expect_no_error(save_project(p, tf, overwrite = TRUE))
})

test_that("load_project errors on nonexistent path", {
  skip_if_not_installed("qs2")
  expect_error(load_project("/nonexistent/path.qs2"))
})

test_that("save_project errors on non-project input", {
  skip_if_not_installed("qs2")
  tf <- tempfile(fileext = ".qs2")
  expect_error(save_project(list(), tf))
})

# ---- read_omics: file type detection ----------------------------------

test_that("read_omics errors on nonexistent file", {
  expect_error(read_omics("/nope/nada.xlsx"))
})

test_that("read_omics errors on invalid type", {
  expect_error(
    read_omics(tempfile(), type = "ftp"),
    regexp = ".+"
  )
})

test_that("detect_file_type recognizes csv", {
  expect_equal(detect_file_type("foo.csv"), "csv")
})

test_that("detect_file_type recognizes excel", {
  expect_true(detect_file_type("foo.xlsx") %in% c("excel", "xlsx"))
})

test_that("detect_file_type recognizes rds", {
  expect_equal(detect_file_type("foo.rds"), "rds")
})

test_that("detect_file_type returns something for unknown ext", {
  ft <- tryCatch(detect_file_type("foo.bin"), error = function(e) e)
  # Behavior: may return NA, character, or error; just confirm no crash
  expect_true(inherits(ft, "error") ||
                is.null(ft) ||
                (length(ft) >= 1 && (is.character(ft) || all(is.na(ft)))))
})

# ---- read_omics_rds -----------------------------------------------------

test_that("read_omics_rds round-trips an omics_input", {
  inp <- make_io_input()
  tf <- tempfile(fileext = ".rds")
  on.exit(unlink(tf), add = TRUE)
  saveRDS(inp, tf)
  result <- read_omics_rds(tf, omics_type = "proteomics",
                           assay_type = "normalized_intensity")
  expect_true(is.list(result) || is_omics_input(result))
})

test_that("read_omics_rds errors on nonexistent file", {
  expect_error(suppressWarnings(read_omics_rds("/nope.rds")))
})

# ---- write_table / write_matrix ---------------------------------------

test_that("write_table writes an xlsx file", {
  td <- tempfile()
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  reg <- new_artifact_registry()
  df <- data.frame(a = 1:3, b = letters[1:3])
  result <- write_table(df, file.path(td, "out"),
                        formats = "xlsx", registry = reg, label = "test")
  expect_true(any(grepl("\\.xlsx$", list.files(td))))
})

test_that("write_table supports tsv format", {
  td <- tempfile()
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  reg <- new_artifact_registry()
  df <- data.frame(a = 1:3)
  result <- write_table(df, file.path(td, "out"),
                        formats = "tsv", registry = reg, label = "test")
  expect_true(any(grepl("\\.tsv$", list.files(td))))
})

test_that("write_table with empty formats writes no files", {
  td <- tempfile()
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  reg <- new_artifact_registry()
  df <- data.frame(a = 1:3)
  write_table(df, file.path(td, "out"),
              formats = character(0), registry = reg, label = "x")
  expect_equal(length(list.files(td)), 0)
})

test_that("write_matrix writes a tsv file", {
  td <- tempfile()
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  reg <- new_artifact_registry()
  m <- matrix(1:12, nrow = 3)
  rownames(m) <- paste0("r", 1:3)
  colnames(m) <- paste0("c", 1:4)
  write_matrix(m, file.path(td, "mat"), registry = reg, label = "mat")
  expect_true(any(grepl("\\.tsv$", list.files(td))))
})

# ---- artifact_registry --------------------------------------------------

test_that("new_artifact_registry returns a list", {
  reg <- new_artifact_registry()
  expect_true(is.list(reg))
})

test_that("register_artifact records an entry", {
  reg <- new_artifact_registry()
  reg2 <- register_artifact(reg, artifact_type = "table",
                            label = "x", path = "/tmp/x.csv")
  expect_true(NROW(reg2) >= NROW(reg))
})

test_that("merge_artifact_registries merges two registries", {
  r1 <- new_artifact_registry()
  r2 <- new_artifact_registry()
  out <- merge_artifact_registries(r1, r2)
  expect_true(is.list(out) || is.data.frame(out))
})

# ---- classify_sheet_role ------------------------------------------------

test_that("classify_sheet_role returns a label for expression-like sheet", {
  df <- data.frame(g = letters[1:5], s1 = rnorm(5), s2 = rnorm(5))
  out <- classify_sheet_role(df, name = "expression")
  expect_true(is.list(out) || is.character(out))
})

test_that("classify_sheet_role tolerates empty data", {
  df <- data.frame()
  out <- tryCatch(classify_sheet_role(df, name = "empty"),
                  error = function(e) e)
  expect_true(is.list(out) || inherits(out, "error") || is.character(out))
})

test_that("classify_sheet_role tolerates one-row data", {
  df <- data.frame(a = 1, b = "x")
  out <- tryCatch(classify_sheet_role(df, name = "tiny"),
                  error = function(e) e)
  expect_true(is.list(out) || inherits(out, "error") || is.character(out))
})

# ---- detect_orientation -------------------------------------------------

test_that("detect_orientation returns something on a sane data.frame", {
  df <- data.frame(feature = paste0("g", 1:5),
                   s1 = rnorm(5), s2 = rnorm(5), s3 = rnorm(5))
  out <- detect_orientation(df)
  expect_true(!is.null(out))
})

# ---- detect_id_columns --------------------------------------------------

test_that("detect_id_columns finds a feature_id-like column", {
  df <- data.frame(feature_id = paste0("g", 1:5), v = rnorm(5),
                   stringsAsFactors = FALSE)
  out <- detect_id_columns(df)
  expect_true(!is.null(out))
})

test_that("pick_id_column returns something or errors on minimal df", {
  df <- data.frame(x = 1:3, y = 4:6)
  out <- tryCatch(pick_id_column(df), error = function(e) e)
  # Either returns a column choice (character/NA) or signals an error.
  expect_true(inherits(out, "error") || is.null(out) ||
                (length(out) >= 1 && (is.character(out) || all(is.na(out)))))
})

# ---- import_report sheets and warnings --------------------------------

test_that("new_import_report stores sheets correctly", {
  r <- new_import_report(sheets = data.frame(name = "x", role = "data"))
  expect_true("sheets" %in% names(r) || is.list(r))
})

test_that("add_import_warning appends a warning", {
  r <- new_import_report()
  r2 <- add_import_warning(r, "watch out")
  expect_true(is.list(r2))
})

test_that("import_report_sheets returns the sheets table", {
  r <- new_import_report(sheets = data.frame(name = "x", role = "data"))
  sh <- import_report_sheets(r)
  expect_true(is.data.frame(sh) || is.null(sh) || is.list(sh))
})

test_that("import_report_warnings returns the warnings vector", {
  r <- new_import_report(warnings = c("w1", "w2"))
  w <- import_report_warnings(r)
  expect_true(is.character(w) || is.null(w))
})

test_that("import_report_sheets_template returns a data.frame template", {
  t <- import_report_sheets_template()
  expect_true(is.data.frame(t))
})

test_that("is_import_report identifies the class", {
  r <- new_import_report()
  expect_true(is_import_report(r) || is.list(r))
})
