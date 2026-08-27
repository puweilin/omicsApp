# =============================================================================
# assay_type vocabulary, scale checking, and proteomics normalization
# =============================================================================
# `assay_type` is the only record of what scale `expr_mat` is on. limma reads
# `expr_mat` untouched; `log2(x + 1)` is applied only for "raw_count". So a
# wrong label does not fail anywhere -- it just changes the answer. These tests
# cover the three things that make the label trustworthy: a vocabulary, a
# numeric cross-check, and a normalization step that moves an input from one
# scale to the other.

# Linear intensities, the scale a proteomics Excel sheet actually carries.
make_linear_input <- function(n_features = 50, n_samples = 4,
                              assay_type = "raw_intensity") {
  set.seed(1)
  expr <- matrix(
    2^rnorm(n_features * n_samples, mean = 20, sd = 2),
    nrow = n_features,
    dimnames = list(paste0("P", seq_len(n_features)),
                    paste0("s", seq_len(n_samples)))
  )
  meta <- data.frame(
    group = rep(c("A", "B"), each = n_samples / 2),
    row.names = colnames(expr)
  )
  feat <- data.frame(
    feature_id = rownames(expr), row.names = rownames(expr)
  )
  # Deliberately not suppressing warnings: several tests below assert on the
  # warning this raises for a mislabelled or superseded assay_type
  omics_input(expr, meta, feat, omics_type = "proteomics",
              assay_type = assay_type)
}

# ---- vocabulary -------------------------------------------------------------

test_that("SUPPORTED_ASSAY_TYPES covers every supported omics type", {
  expect_type(SUPPORTED_ASSAY_TYPES, "list")
  expect_setequal(names(SUPPORTED_ASSAY_TYPES), SUPPORTED_OMICS_TYPES)
  expect_true(all(vapply(SUPPORTED_ASSAY_TYPES, is.character, logical(1))))
  expect_true("raw_intensity" %in% SUPPORTED_ASSAY_TYPES$proteomics)
  expect_true("normalized_intensity" %in% SUPPORTED_ASSAY_TYPES$proteomics)
  expect_true("raw_count" %in% SUPPORTED_ASSAY_TYPES$rnaseq)
})

test_that("every log-scale assay type is itself a known assay type", {
  known <- unlist(SUPPORTED_ASSAY_TYPES, use.names = FALSE)
  expect_true(all(LOG_SCALE_ASSAY_TYPES %in% known))
})

test_that("deprecated aliases point at values in the vocabulary", {
  known <- unlist(SUPPORTED_ASSAY_TYPES, use.names = FALSE)
  expect_true(all(unname(DEPRECATED_ASSAY_TYPE_ALIASES) %in% known))
  # The import view stamped "intensity" on every proteomics upload
  expect_identical(unname(DEPRECATED_ASSAY_TYPE_ALIASES[["intensity"]]),
                   "raw_intensity")
  # The docs said "raw_counts" while the code only matched "raw_count"
  expect_identical(unname(DEPRECATED_ASSAY_TYPE_ALIASES[["raw_counts"]]),
                   "raw_count")
})

test_that("unknown assay types warn but are still accepted", {
  # CyTOF drives omicsCore with arcsinh-transformed values; it must keep working
  expect_warning(
    input <- make_linear_input(assay_type = "arcsinh_intensity"),
    "Unrecognised `assay_type`"
  )
  expect_identical(input$assay_type, "arcsinh_intensity")
})

test_that("a missing assay type warns", {
  set.seed(1)
  expr <- matrix(rnorm(20), nrow = 5, dimnames = list(
    paste0("P", 1:5), paste0("s", 1:4)))
  meta <- data.frame(group = rep(c("A", "B"), each = 2), row.names = colnames(expr))
  feat <- data.frame(feature_id = rownames(expr), row.names = rownames(expr))
  expect_warning(
    omics_input(expr, meta, feat, omics_type = "proteomics", assay_type = NULL),
    "`assay_type` is missing"
  )
})

# ---- deprecated aliases -----------------------------------------------------

test_that("omics_input rewrites a deprecated alias", {
  expect_warning(
    input <- make_linear_input(assay_type = "intensity"),
    "superseded"
  )
  expect_identical(input$assay_type, "raw_intensity")
})

test_that("validate_omics_input accepts a stored alias without calling it unknown", {
  # Projects saved before the vocabulary existed carry assay_type = "intensity"
  # and never pass back through the constructor
  input <- make_linear_input()
  input$assay_type <- "intensity"

  expect_warning(validate_omics_input(input), "superseded")
  expect_no_warning(
    suppressWarnings(validate_omics_input(input)) # not "Unrecognised"
  )
  w <- tryCatch(validate_omics_input(input), warning = function(w) conditionMessage(w))
  expect_false(grepl("Unrecognised", w))
})

# ---- scale checking ---------------------------------------------------------

test_that("check_assay_scale passes when values match the declared scale", {
  expect_silent(check_assay_scale(make_linear_input()))
  expect_true(check_assay_scale(make_linear_input()))
})

test_that("check_assay_scale catches log-scale values labelled linear", {
  input <- make_linear_input()
  input$expr_mat <- log2(input$expr_mat)   # now ~20, labelled raw_intensity

  expect_warning(res <- check_assay_scale(input), "looks .*already log")
  expect_false(res)
})

test_that("check_assay_scale catches linear values labelled normalized", {
  # Values are ~1e6 but the label promises glog2 output, so construction
  # already warns; the point here is that a direct call warns too
  input <- suppressWarnings(make_linear_input(assay_type = "normalized_intensity"))

  expect_warning(res <- check_assay_scale(input), "implies log-scale")
  expect_false(res)
})

test_that("check_assay_scale stays quiet for assay types it cannot place", {
  input <- suppressWarnings(make_linear_input(assay_type = "arcsinh_intensity"))
  expect_silent(check_assay_scale(input))
})

test_that("omics_input runs the scale check once at construction", {
  set.seed(1)
  expr <- matrix(rnorm(200, mean = 20, sd = 2), nrow = 50, dimnames = list(
    paste0("P", 1:50), paste0("s", 1:4)))
  meta <- data.frame(group = rep(c("A", "B"), each = 2), row.names = colnames(expr))
  feat <- data.frame(feature_id = rownames(expr), row.names = rownames(expr))

  expect_warning(
    omics_input(expr, meta, feat, omics_type = "proteomics",
                assay_type = "raw_intensity"),
    "already log"
  )
})

# ---- normalization ----------------------------------------------------------

test_that("normalize_omics with log2 moves the input onto a log scale", {
  input <- make_linear_input()
  out <- suppressMessages(normalize_omics(input, method = "log2"))

  expect_identical(out$assay_type, "normalized_intensity")
  expect_lt(max(out$expr_mat, na.rm = TRUE), MAX_PLAUSIBLE_LOG_SCALE_VALUE)
  expect_identical(dimnames(out$expr_mat), dimnames(input$expr_mat))
  # The declared scale and the values now agree
  expect_silent(check_assay_scale(out))
})

test_that("normalize_omics keeps the pre-normalization matrix in raw_mat", {
  input <- make_linear_input()
  out <- suppressMessages(normalize_omics(input, method = "log2"))

  expect_identical(out$raw_mat, input$expr_mat)
  expect_identical(out$normalized_mat, out$expr_mat)
})

test_that("normalize_omics with vsn returns glog2-scale values", {
  skip_if_not_installed("vsn")

  input <- make_linear_input(n_features = 200, n_samples = 8)
  out <- suppressMessages(normalize_omics(input, method = "vsn"))

  expect_identical(out$assay_type, "normalized_intensity")
  expect_identical(dim(out$expr_mat), dim(input$expr_mat))
  expect_silent(check_assay_scale(out))
  # glog2 tracks log2 at the high end, so a matrix built around 2^20 lands near 20
  expect_gt(max(out$expr_mat, na.rm = TRUE), 10)
  expect_lt(max(out$expr_mat, na.rm = TRUE), MAX_PLAUSIBLE_LOG_SCALE_VALUE)
})

test_that("normalize_omics refuses to normalize twice", {
  input <- make_linear_input()
  out <- suppressMessages(normalize_omics(input, method = "log2"))

  expect_error(normalize_omics(out, method = "log2"), "already")
  # And every proteomics log-scale label is refused, not just the one just set
  proteomics_log <- intersect(LOG_SCALE_ASSAY_TYPES,
                              SUPPORTED_ASSAY_TYPES$proteomics)
  expect_gt(length(proteomics_log), 1)
  for (at in proteomics_log) {
    marked <- input
    marked$assay_type <- at
    expect_error(normalize_omics(marked, method = "log2"), "already")
  }
})

test_that("normalize_omics treats zero and negative intensities as missing", {
  input <- make_linear_input()
  input$expr_mat[1, 1] <- 0
  input$expr_mat[2, 1] <- -5

  expect_warning(
    out <- suppressMessages(normalize_omics(input, method = "log2")),
    "negative"
  )
  expect_true(is.na(out$expr_mat[1, 1]))
  expect_true(is.na(out$expr_mat[2, 1]))
})

test_that("normalize_omics rejects RNA-seq with an explanation", {
  set.seed(1)
  expr <- matrix(rpois(200, 200), nrow = 50, dimnames = list(
    paste0("G", 1:50), paste0("s", 1:4)))
  meta <- data.frame(group = rep(c("A", "B"), each = 2), row.names = colnames(expr))
  feat <- data.frame(feature_id = rownames(expr), row.names = rownames(expr))
  input <- omics_input(expr, meta, feat, omics_type = "rnaseq",
                       assay_type = "raw_count")

  expect_error(normalize_omics(input), "proteomics only")
})

test_that("normalize_omics output survives validation and subsetting", {
  input <- make_linear_input()
  out <- suppressMessages(normalize_omics(input, method = "log2"))

  expect_silent(validate_omics_input(out))
  subset <- subset_omics(out, samples = colnames(out$expr_mat)[1:2])
  expect_identical(subset$assay_type, "normalized_intensity")
  expect_identical(ncol(subset$expr_mat), 2L)
})
