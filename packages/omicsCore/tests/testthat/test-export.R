# export_bundle() + export_report() tests. Reuses the diff test helper.

make_export_input <- function(n_features = 30L, n_per_group = 4L) {
  set.seed(2028)
  n_samples <- 2L * n_per_group
  expr <- matrix(
    rnorm(n_features * n_samples, mean = 5, sd = 1.0),
    nrow = n_features,
    dimnames = list(paste0("f", seq_len(n_features)),
                    paste0("s", seq_len(n_samples)))
  )
  expr[1:5, seq_len(n_per_group)] <- expr[1:5, seq_len(n_per_group)] - 2
  meta <- data.frame(group = rep(c("ctrl", "case"), each = n_per_group),
                     row.names = paste0("s", seq_len(n_samples)),
                     stringsAsFactors = FALSE)
  feat <- data.frame(feature_id = paste0("f", seq_len(n_features)),
                     feature_symbol = paste0("GENE", seq_len(n_features)),
                     row.names = paste0("f", seq_len(n_features)),
                     stringsAsFactors = FALSE)
  omics_input(expr, meta, feat, omics_type = "proteomics",
              assay_type = "normalized_intensity")
}

make_export_diff_bundle <- function() {
  x <- make_export_input()
  run_diff(x, method = "ttest", analysis_type = "group",
           group_col = "group", control_group = "ctrl", case_group = "case")
}

# ---- export_bundle ----------------------------------------------------

test_that("export_bundle writes tables + params for a diff bundle", {
  b <- make_export_diff_bundle()
  dir <- tempfile()
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)

  reg <- export_bundle(b, dir = dir, formats = c("xlsx", "tsv"))
  expect_s3_class(reg, "data.frame")
  # At minimum: bundle_params + the diff_result_df.
  expect_true("params" %in% reg$artifact_type)
  expect_true("table" %in% reg$artifact_type)

  expect_true(file.exists(file.path(dir, "run_diff_diff_result_df.tsv")))
  expect_true(file.exists(file.path(dir, "run_diff_diff_result_df.xlsx")))
  expect_true(file.exists(file.path(dir, "run_diff_params.json")))

  json <- jsonlite::fromJSON(file.path(dir, "run_diff_params.json"))
  expect_identical(json$analysis_name, "run_diff")
})

test_that("export_bundle writes plots when supplied", {
  b <- make_export_diff_bundle()
  dir <- tempfile()
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)

  reg <- export_bundle(
    b, dir = dir, formats = c("tsv", "pdf"),
    plots = list(volcano = plot_volcano(b, top_n = 5))
  )
  expect_true(any(reg$artifact_type == "plot"))
  expect_true(file.exists(file.path(dir, "run_diff_volcano.pdf")))
})

test_that("export_bundle honors prefix", {
  b <- make_export_diff_bundle()
  dir <- tempfile()
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)

  export_bundle(b, dir = dir, formats = "tsv", prefix = "demo_")
  expect_true(file.exists(file.path(dir, "demo_run_diff_diff_result_df.tsv")))
})

test_that("export_bundle errors on non-bundle inputs", {
  expect_error(export_bundle(list(), tempfile()), "analysis_bundle")
})

# ---- export_report ----------------------------------------------------

test_that("export_report renders an HTML file", {
  skip_if_not_installed("rmarkdown")
  skip_on_cran()
  if (!nzchar(Sys.which("pandoc")) && !rmarkdown::pandoc_available()) {
    skip("pandoc unavailable")
  }
  x <- make_export_input()
  p <- omics_project("export_demo", experiments = list(proteo = x))
  p$bundles <- list(diff = run_diff(x, method = "ttest",
                                    analysis_type = "group",
                                    group_col = "group",
                                    control_group = "ctrl",
                                    case_group = "case"))
  tmp <- tempfile(fileext = ".html")
  on.exit(unlink(tmp), add = TRUE)
  out <- export_report(p, tmp)
  expect_identical(out, tmp)
  expect_true(file.exists(tmp))
  expect_gt(file.info(tmp)$size, 0L)
})

test_that("export_report errors on non-projects", {
  expect_error(export_report(list(), tempfile(fileext = ".html")),
               "omics_project")
})

test_that("export_report refuses to overwrite by default", {
  skip_if_not_installed("rmarkdown")
  skip_on_cran()
  if (!nzchar(Sys.which("pandoc")) && !rmarkdown::pandoc_available()) {
    skip("pandoc unavailable")
  }
  x <- make_export_input()
  p <- omics_project("export_demo", experiments = list(proteo = x))
  tmp <- tempfile(fileext = ".html")
  on.exit(unlink(tmp), add = TRUE)
  export_report(p, tmp)
  expect_error(export_report(p, tmp), "already exists")
})
