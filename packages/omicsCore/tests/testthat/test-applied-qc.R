qc_project_fixture <- function() {
  set.seed(62)
  x <- matrix(rnorm(120, 20), 20, 6,
              dimnames = list(paste0("G", 1:20), paste0("S", 1:6)))
  x[1, 1:4] <- NA_real_
  x[2, 1] <- NA_real_
  inp <- omics_input(x, data.frame(group = rep(c("A", "B"), each = 3), row.names = colnames(x)),
                     data.frame(feature_id = rownames(x), row.names = rownames(x)),
                     omics_type = "proteomics", assay_type = "normalized_intensity")
  omics_project("QC application", list(protein = inp))
}

test_that("applied QC reaches analyses, survives reload and can be relaxed", {
  p <- qc_project_fixture()
  base <- p$experiments$protein
  q <- run_qc(base, missing_threshold = 0.5, impute_method = "min", outlier_method = "none")
  p$bundles <- list(diff = "old result")
  applied <- apply_qc(p, "protein", q)
  expect_equal(nrow(applied$experiments$protein$expr_mat), 19)
  expect_false(anyNA(applied$experiments$protein$expr_mat))
  expect_identical(qc_baseline_input(applied, "protein"), base)
  expect_null(applied$bundles$diff)
  diff <- run_diff(applied$experiments$protein, method = "ttest", group_col = "group",
                   control_group = "A", case_group = "B")
  expect_equal(nrow(diff$results$diff_result_df), 19)
  path <- withr::local_tempfile(fileext = ".omp")
  save_project(applied, path)
  restored <- load_project(path)
  expect_identical(restored$experiments, applied$experiments)
  loose <- run_qc(qc_baseline_input(restored, "protein"), missing_threshold = 0.9,
                  impute_method = "min", outlier_method = "none")
  expect_equal(nrow(apply_qc(restored, "protein", loose)$experiments$protein$expr_mat), 20)
  # Replacing a layer must not resurrect its previous import.
  restored$experiments$protein <- base
  restored$experiments$protein$expr_mat[] <- 1
  expect_identical(qc_baseline_input(restored, "protein"), restored$experiments$protein)
  q$input_info$layer <- "other"
  expect_error(apply_qc(p, "protein", q), "different layer")
})

test_that("legacy desktop preprocessing migrates only for the matching import", {
  p <- qc_project_fixture()
  p$experiments$protein$version <- "old-version"
  q <- run_qc(p$experiments$protein, missing_threshold = 0.5,
              impute_method = "min", outlier_method = "none")
  p$preprocessing <- list(protein = list(base_input_version = "old-version",
    cleaned = qc_cleaned_input(q, p$experiments$protein), params = q$params))
  path <- withr::local_tempfile(fileext = ".omp")
  save_project(p, path)
  migrated <- load_project(path)
  expect_equal(nrow(migrated$experiments$protein$expr_mat), 19)
  expect_equal(nrow(qc_baseline_input(migrated, "protein")$expr_mat), 20)
  expect_identical(migrate_desktop_qc(migrated), migrated)
  p$preprocessing$protein$base_input_version <- "stale"
  expect_identical(migrate_desktop_qc(p)$experiments, p$experiments)
})

test_that("applied cleaning preserves upload identity and excluded samples", {
  p <- qc_project_fixture()
  p$experiments$protein$source_fingerprint <- "123456789abc:file"
  q <- run_qc(p$experiments$protein, impute_method = "min", outlier_method = "none")
  p <- apply_qc(p, "protein", q)
  expect_match(p$experiments$protein$source_fingerprint, "^123456789abc:file:qc=")
  p$experiments$protein <- subset_omics_samples(p$experiments$protein, paste0("S", 2:6))
  p$experiments$protein$excluded_samples <- "S1"
  expect_false("S1" %in% colnames(qc_baseline_input(p, "protein")$expr_mat))
  expect_false(is.null(p$experiments$protein$applied_qc))
})
