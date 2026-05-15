# Tests for save_project, load_project, export_report, export_bundle.
# These verify round-trip persistence and export capabilities.

make_persist_project <- function() {
  mat <- matrix(rnorm(100, mean = 10, sd = 2), nrow = 20, ncol = 5)
  rownames(mat) <- paste0("gene_", 1:20)
  colnames(mat) <- paste0("s", 1:5)
  meta <- data.frame(group = c("G1", "G1", "G2", "G2", "G2"),
                     row.names = paste0("s", 1:5))
  feat <- data.frame(feature_id = paste0("gene_", 1:20),
                     feature_name = paste0("Gene", 1:20),
                     stringsAsFactors = FALSE)
  inp <- omics_input(mat, meta, feat, omics_type = "proteomics",
                     assay_type = "intensity")
  proj <- omics_project("persist_test")
  proj <- add_experiment(proj, "prot", inp)
  proj
}

make_persist_input <- function() {
  mat <- matrix(rnorm(100, mean = 10, sd = 2), nrow = 20, ncol = 5)
  rownames(mat) <- paste0("gene_", 1:20)
  colnames(mat) <- paste0("s", 1:5)
  meta <- data.frame(group = c("G1", "G1", "G2", "G2", "G2"),
                     row.names = paste0("s", 1:5))
  feat <- data.frame(feature_id = paste0("gene_", 1:20),
                     feature_name = paste0("Gene", 1:20),
                     stringsAsFactors = FALSE)
  omics_input(mat, meta, feat, omics_type = "proteomics",
              assay_type = "intensity")
}

# ---- save_project / load_project ---------------------------------------

test_that("save_project writes a file and load_project reads it back", {
  proj <- make_persist_project()
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp, force = TRUE))
  save_project(proj, tmp)
  expect_true(file.exists(tmp))
  loaded <- load_project(tmp)
  expect_s3_class(loaded, "omics_project")
  expect_equal(loaded$name, proj$name)
  expect_equal(experiment_tags(loaded), experiment_tags(proj))
})

test_that("save_project overwrites existing file with overwrite=TRUE", {
  proj <- make_persist_project()
  tmp <- tempfile(fileext = ".omp")
  on.exit(unlink(tmp, force = TRUE))
  save_project(proj, tmp)
  save_project(proj, tmp, overwrite = TRUE)  # should not error
  expect_true(file.exists(tmp))
})

test_that("save_project errors when file exists and overwrite=FALSE", {
  proj <- make_persist_project()
  tmp <- tempfile(fileext = ".omp")
  on.exit(unlink(tmp, force = TRUE))
  save_project(proj, tmp)
  expect_error(save_project(proj, tmp), "already exists")
})

test_that("save_project errors on non-omics_project", {
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp, force = TRUE))
  expect_error(save_project(list(), tmp), "omics_project")
})

test_that("load_project errors on nonexistent file", {
  expect_error(load_project("nonexistent_file_12345.rds"), "file")
})

test_that("load_project errors on invalid RDS", {
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp, force = TRUE))
  saveRDS("not a project", tmp)
  expect_error(load_project(tmp))
})

test_that("save/load round-trip preserves experiment data", {
  proj <- make_persist_project()
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp, force = TRUE))
  save_project(proj, tmp)
  loaded <- load_project(tmp)
  orig <- proj$experiments$prot
  restored <- loaded$experiments$prot
  expect_equal(dim(orig$expr_mat), dim(restored$expr_mat))
  expect_equal(orig$omics_type, restored$omics_type)
})

test_that("save/load round-trip preserves metadata", {
  proj <- make_persist_project()
  proj$metadata <- list(author = "test", version = 2)
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp, force = TRUE))
  save_project(proj, tmp)
  loaded <- load_project(tmp)
  expect_equal(loaded$metadata$author, "test")
  expect_equal(loaded$metadata$version, 2)
})

test_that("save/load round-trip preserves sample_link", {
  proj <- make_persist_project()
  proj$sample_link <- data.frame(
    tag = "prot", sample_id = colnames(proj$experiments$prot$expr_mat)[1:2],
    donor_id = c("D1", "D2"), stringsAsFactors = FALSE
  )
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp, force = TRUE))
  save_project(proj, tmp)
  loaded <- load_project(tmp)
  expect_equal(nrow(loaded$sample_link), 2)
  expect_equal(loaded$sample_link$donor_id, c("D1", "D2"))
})

# ---- export_bundle -----------------------------------------------------

test_that("export_bundle writes tabular files to dir", {
  inp <- make_persist_input()
  bundle <- run_diff(inp, method = "ttest", analysis_type = "group",
                     group_col = "group",
                     control_group = "G1", case_group = "G2")
  tmpdir <- tempfile()
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE, force = TRUE))
  res <- export_bundle(bundle, tmpdir, formats = "tsv")
  # any tsv files written?
  files <- list.files(tmpdir, pattern = "\\.tsv$", recursive = TRUE)
  expect_gt(length(files), 0)
})

test_that("export_bundle creates parent directories if needed", {
  inp <- make_persist_input()
  bundle <- run_diff(inp, method = "ttest", analysis_type = "group",
                     group_col = "group",
                     control_group = "G1", case_group = "G2")
  tmpdir <- file.path(tempdir(), "deep", "nested", "dir")
  on.exit(unlink(file.path(tempdir(), "deep"), recursive = TRUE, force = TRUE))
  export_bundle(bundle, tmpdir, formats = "tsv")
  expect_true(dir.exists(tmpdir))
})

test_that("export_bundle errors on non-analysis_bundle", {
  tmpdir <- tempfile()
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE, force = TRUE))
  expect_error(export_bundle(list(), tmpdir), "analysis_bundle")
})

# ---- export_report -----------------------------------------------------

test_that("export_report HTML creates a file", {
  skip_if_not_installed("rmarkdown")
  proj <- make_persist_project()
  tmp <- tempfile(fileext = ".html")
  on.exit(unlink(tmp, force = TRUE))
  export_report(proj, tmp, format = "html")
  expect_true(file.exists(tmp))
  expect_gt(file.info(tmp)$size, 0)
})

test_that("export_report HTML with bundles includes analysis sections", {
  skip_if_not_installed("rmarkdown")
  proj <- make_persist_project()
  bundle <- run_diff(proj$experiments$prot, method = "ttest",
                     analysis_type = "group", group_col = "group",
                     control_group = "G1", case_group = "G2")
  proj$bundles <- list(diff = bundle)
  tmp <- tempfile(fileext = ".html")
  on.exit(unlink(tmp, force = TRUE))
  export_report(proj, tmp, format = "html")
  expect_true(file.exists(tmp))
})

test_that("export_report errors on non-omics_project", {
  skip_if_not_installed("rmarkdown")
  tmp <- tempfile(fileext = ".html")
  on.exit(unlink(tmp, force = TRUE))
  expect_error(export_report(list(), tmp), "omics_project")
})

test_that("export_report errors on invalid format", {
  skip_if_not_installed("rmarkdown")
  proj <- make_persist_project()
  tmp <- tempfile(fileext = ".invalid")
  on.exit(unlink(tmp, force = TRUE))
  expect_error(export_report(proj, tmp, format = "invalid"), "should be one of")
})

test_that("export_report accepts overwrite argument", {
  skip_if_not_installed("rmarkdown")
  proj <- make_persist_project()
  tmp <- tempfile(fileext = ".html")
  on.exit(unlink(tmp, force = TRUE))
  export_report(proj, tmp, format = "html", overwrite = TRUE)
  # second write with overwrite = TRUE should succeed
  expect_no_error(export_report(proj, tmp, format = "html", overwrite = TRUE))
})
