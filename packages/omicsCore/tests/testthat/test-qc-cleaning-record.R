# The QC bundle keeps a record of how the input was cleaned instead of a
# second copy of it, and MinProb falls back to MinDet when it cannot be
# estimated.

cr_input <- function(n_feat = 400, n_samp = 12, missing = 0.15,
                     assay_type = "normalized_intensity", seed = 11) {
  set.seed(seed)
  mu <- stats::rnorm(n_feat, 22, 2)
  m <- matrix(stats::rnorm(n_feat * n_samp, mu, 0.6), n_feat, n_samp,
              dimnames = list(paste0("P", seq_len(n_feat)), paste0("S", seq_len(n_samp))))
  m[sample(length(m), round(missing * length(m)))] <- NA
  if (n_samp >= 10L) m[1:5, 1:10] <- NA   # features the default filter drops
  if (!assay_type %in% LOG_SCALE_ASSAY_TYPES) m <- 2^m
  meta <- data.frame(group = rep(c("A", "B"), length.out = n_samp),
                     row.names = colnames(m))
  feat <- data.frame(feature_id = rownames(m), symbol = paste0("G", seq_len(n_feat)),
                     row.names = rownames(m))
  omics_input(m, meta, feat, omics_type = "proteomics", assay_type = assay_type)
}

# What run_qc() used to store as results$cleaned_input, built step by step.
cr_cleaned_by_hand <- function(input, b, method) {
  keep_s <- b$results$cleaning$kept_samples
  keep_f <- setdiff(rownames(input$expr_mat), b$results$qc_summary$recommended_filters$remove_features)
  cl <- subset_omics(input, samples = keep_s, features = keep_f)
  cl$raw_mat <- cl$raw_mat %||% cl$expr_mat
  cl$expr_mat <- impute_matrix(cl$expr_mat, method = method)
  cl$assay_type <- "imputed_intensity"
  cl
}

test_that("the bundle carries a record, not a copy, and the copy is rebuilt exactly", {
  inp <- cr_input()
  b <- run_qc(inp, impute_method = "min", outlier_method = "pca")
  expect_null(b$results$cleaned_input)
  expect_named(b$results, c("qc_summary", "cleaning", "plot_data"))
  rec <- b$results$cleaning
  expect_identical(rec$dropped_features, paste0("P", 1:5))
  expect_identical(rec$kept_samples, colnames(inp$expr_mat))
  expect_length(rec$imputed_values, b$results$qc_summary$imputation$n_imputed)

  cl <- qc_cleaned_input(b, inp)
  expect_true(is_omics_input(cl))
  expect_identical(cl, cr_cleaned_by_hand(inp, b, "min"))
  expect_identical(dim(cl$expr_mat), c(b$input_info$n_features_out, b$input_info$n_samples_out))
})

test_that("the rebuilt input follows the linear-scale path too", {
  skip_if_not_installed("imputeLCMD")
  inp <- cr_input(assay_type = "raw_intensity")
  na <- is.na(inp$expr_mat[-(1:5), ])
  b <- run_qc(inp, impute_method = "MinProb", outlier_method = "none")
  cl <- qc_cleaned_input(b, inp)
  expect_identical(cl$assay_type, "raw_intensity")
  expect_false(anyNA(cl$expr_mat))
  expect_true(all(cl$expr_mat[na] > 0))
  # Observed values untouched; raw_mat is the matrix before imputation
  expect_identical(cl$expr_mat[!na], inp$expr_mat[-(1:5), ][!na])
  expect_identical(cl$raw_mat, inp$expr_mat[-(1:5), ])
})

test_that("a removed outlier stays removed in the rebuilt input", {
  inp <- cr_input(missing = 0)
  inp$expr_mat[, 1] <- inp$expr_mat[, 1] + 10
  b <- run_qc(inp, outlier_method = "iqr", outlier_sd_threshold = 1.5,
              impute_method = "none", remove_outliers = TRUE)
  cl <- qc_cleaned_input(b, inp)
  expect_false("S1" %in% colnames(cl$expr_mat))
  expect_identical(b$results$cleaning$kept_samples, colnames(cl$expr_mat))
  expect_null(b$results$qc_summary$imputation)
})

test_that("qc_cleaned_input refuses an input QC did not run on, and says what it needs", {
  inp <- cr_input()
  b <- run_qc(inp, impute_method = "min", outlier_method = "none")
  expect_error(qc_cleaned_input(b), "pass it as `input`", fixed = TRUE)
  other <- inp
  other$expr_mat[10, 3] <- other$expr_mat[10, 3] + 1
  expect_error(qc_cleaned_input(b, other), "not the data this QC result was computed on")
  expect_error(qc_cleaned_input(list(), inp), "analysis_bundle from run_qc")
  # The design and annotation may change; the numbers may not.
  relabelled <- set_study_design(inp, "group")
  expect_no_error(qc_cleaned_input(b, relabelled))
})

test_that("the record survives saving and opening a project", {
  skip_if_not_installed("qs2")
  inp <- cr_input()
  b <- run_qc(inp, impute_method = "min", outlier_method = "none")
  b$input_info$layer <- "prot"
  proj <- omics_project("rec", experiments = list(prot = inp))
  proj$bundles <- list(qc = b)
  path <- tempfile(fileext = ".omp")
  on.exit(unlink(path), add = TRUE)
  save_project(proj, path)
  back <- load_project(path)
  expect_identical(qc_cleaned_input(back$bundles$qc, back$experiments$prot),
                   qc_cleaned_input(b, inp))
  expect_s3_class(plot_qc(back$bundles$qc, view = "pca"), "ggplot")
})

test_that("the bundle is a fraction of the size it was", {
  inp <- suppressMessages(normalize_omics(cr_input(n_feat = 2000, n_samp = 30,
                                                   assay_type = "raw_intensity"),
                                          method = "log2"))
  b <- run_qc(inp, impute_method = "min", outlier_method = "pca")
  old <- b
  old$results$cleaned_input <- qc_cleaned_input(b, inp)
  new_size <- as.numeric(utils::object.size(b))
  old_size <- as.numeric(utils::object.size(old))
  # A quarter at most; on 8,000 x 60 with 15% missing it is a tenth.
  expect_lt(new_size, old_size / 4)
})

test_that("plots drawn from the record match those drawn from a stored copy", {
  inp <- cr_input()
  b <- run_qc(inp, impute_method = "min", outlier_method = "pca")
  old <- b
  old$results$plot_data <- NULL
  old$results$cleaning <- NULL
  old$results$cleaned_input <- qc_cleaned_input(b, inp)

  p_new <- plot_qc(b, view = "pca", color_by = "group")
  p_old <- plot_qc(old, view = "pca", color_by = "group")
  expect_equal(p_new$data, p_old$data)
  expect_identical(p_new$labels$x, p_old$labels$x)

  c_new <- plot_qc(b, view = "connectivity")
  c_old <- plot_qc(old, view = "connectivity")
  expect_equal(c_new$data, c_old$data)

  i_new <- plot_qc(b, view = "imputation")
  i_old <- plot_qc(old, view = "imputation")
  expect_setequal(unique(i_new$data$type), c("observed", "imputed"))
  # The imputed curve is drawn from the same values either way
  expect_equal(i_new$data[i_new$data$type == "imputed", ],
               i_old$data[i_old$data$type == "imputed", ], ignore_attr = TRUE)
  # The observed one from a summary, close to the full set
  o_new <- i_new$data[i_new$data$type == "observed", ]
  o_old <- i_old$data[i_old$data$type == "observed", ]
  expect_lt(max(abs(o_new$density - o_old$density)), 0.05 * max(o_old$density))

  expect_error(plot_qc(b, view = "pca", color_by = "nope"), "not found")
})

test_that("a PCA that cannot be drawn says why when it is asked for", {
  inp <- cr_input(n_samp = 1, missing = 0)
  b <- run_qc(inp, outlier_method = "none", impute_method = "none",
              missing_threshold = 1)
  expect_s3_class(plot_qc(b, view = "missing"), "ggplot")
  expect_error(plot_qc(b, view = "pca"), "at least 2 samples")
  expect_error(plot_qc(b, view = "imputation"), "no imputation step")
})

test_that("a result saved with its cleaned input still opens, plots and rebuilds", {
  skip_if_not_installed("qs2")
  f <- testthat::test_path("fixtures", "omp",
                           "schema-1.0.0__omicsCore-0.2.0__qc-full-cleaned-input.omp")
  skip_if_not(file.exists(f))
  proj <- load_project(f)
  qc <- proj$bundles$qc
  expect_true(is_omics_input(qc$results$cleaned_input))
  expect_null(qc$results$cleaning)
  for (v in c("missing", "depth", "pca", "connectivity", "imputation")) {
    expect_s3_class(plot_qc(qc, view = v), "ggplot")
  }
  expect_identical(qc_cleaned_input(qc), qc$results$cleaned_input)
  expect_identical(qc_cleaned_input(qc, proj$experiments$proteomics),
                   qc$results$cleaned_input)
  # Run again today, the same data gives the same cleaned input
  again <- run_qc(proj$experiments$proteomics)
  expect_identical(qc_cleaned_input(again, proj$experiments$proteomics),
                   qc$results$cleaned_input)
  expect_no_error(parse(text = export_script(proj)))
})

# ---- MinProb that cannot be estimated -----------------------------------

test_that("MinProb with fewer than two well-measured features falls back to MinDet", {
  skip_if_not_installed("imputeLCMD")
  # One feature seen in more than half the samples: imputeLCMD stops
  m1 <- rbind(a = c(20, 21, 22, 23), b = c(NA, NA, NA, 19), c = c(NA, NA, 18, NA))
  colnames(m1) <- paste0("s", 1:4)
  # None: imputeLCMD returns every gap still empty
  m0 <- rbind(a = c(20, NA, NA, NA), b = c(NA, NA, 19, NA), c = c(NA, 18, NA, 21))
  colnames(m0) <- paste0("s", 1:4)
  for (m in list(m1, m0)) {
    expect_warning(out <- impute_matrix(m, method = "MinProb"),
                   class = "omics_impute_fallback")
    expect_false(anyNA(out))
    expect_identical(out, impute_matrix(m, method = "MinDet"))
    expect_identical(dimnames(out), dimnames(m))
  }
  w <- tryCatch(impute_matrix(m1, method = "MinProb"), warning = function(w) w)
  expect_match(conditionMessage(w), "MinDet")
  expect_match(conditionMessage(w), "this data has 1")
})

test_that("MinProb on two or more well-measured features is MinProb, unchanged", {
  skip_if_not_installed("imputeLCMD")
  m <- rbind(a = c(20, 21, 22, 23), d = c(19, 20, NA, 21), b = c(NA, NA, NA, 19))
  colnames(m) <- paste0("s", 1:4)
  expect_no_warning(out <- impute_matrix(m, method = "MinProb"))
  expect_false(anyNA(out))
  expect_false(identical(out, impute_matrix(m, method = "MinDet")))
})

test_that("run_qc notes the fallback and records which method ran", {
  skip_if_not_installed("imputeLCMD")
  inp <- cr_input(n_feat = 6, n_samp = 6, missing = 0, seed = 3)
  m <- inp$expr_mat
  m[2:6, 1:3] <- NA     # only P1 seen in more than half the samples
  inp$expr_mat <- m
  expect_no_warning(b <- run_qc(inp, missing_threshold = 1, outlier_method = "none"))
  expect_identical(b$params$impute_method, "MinProb")
  imp <- b$results$qc_summary$imputation
  expect_identical(imp$method, "MinDet")
  expect_identical(imp$requested_method, "MinProb")
  expect_match(b$warnings, "MinProb could not be used", all = FALSE)
  cl <- qc_cleaned_input(b, inp)
  expect_false(anyNA(cl$expr_mat))

  # Enough features: MinProb runs and nothing is said
  b2 <- run_qc(cr_input(), outlier_method = "none")
  expect_identical(b2$results$qc_summary$imputation$method, "MinProb")
  expect_null(b2$results$qc_summary$imputation$requested_method)
  expect_false(any(grepl("MinProb could not", b2$warnings)))
})

test_that("the report draws QC from the record and names the method that ran", {
  skip_if_not_installed("rmarkdown")
  skip_if_not(rmarkdown::pandoc_available(), "pandoc unavailable")
  skip_if_not_installed("imputeLCMD")
  inp <- cr_input(n_feat = 6, n_samp = 6, missing = 0, seed = 3)
  m <- inp$expr_mat
  m[2:6, 1:3] <- NA
  inp$expr_mat <- m
  b <- run_qc(inp, missing_threshold = 1, outlier_method = "none")
  b$input_info$layer <- "prot"
  proj <- omics_project("r", experiments = list(prot = inp))
  proj$bundles <- list(qc = b)
  path <- tempfile(fileext = ".html")
  on.exit(unlink(path), add = TRUE)
  suppressMessages(export_report(proj, path, format = "html"))
  html <- paste(readLines(path, warn = FALSE), collapse = "\n")
  # Pandoc wraps lines, so words are separated by any whitespace.
  expect_match(html, "imputed\\s+with\\s+MinDet")
  expect_match(html, "MinProb\\s+could\\s+not\\s+be\\s+estimated")
  expect_match(html, "Quality control", fixed = TRUE)
  expect_match(html, "<img", fixed = TRUE)
})
