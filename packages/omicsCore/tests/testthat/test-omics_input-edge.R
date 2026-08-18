# Edge-case tests for omics_input construction and validation.

# ---- omics_input: construction -----------------------------------------

test_that("omics_input creates a valid object from matrix input", {
  mat <- matrix(1:20, nrow = 5, dimnames = list(paste0("g", 1:5), paste0("s", 1:4)))
  meta <- data.frame(group = rep(c("A", "B"), each = 2), row.names = paste0("s", 1:4))
  feat <- data.frame(feature_id = paste0("g", 1:5), name = paste0("Gene", 1:5),
                     stringsAsFactors = FALSE)
  inp <- omics_input(mat, meta, feat, omics_type = "proteomics",
                     assay_type = "intensity")
  expect_s3_class(inp, "omics_input")
  expect_equal(inp$omics_type, "proteomics")
  expect_equal(inp$assay_type, "intensity")
  expect_equal(dim(inp$expr_mat), c(5, 4))
})

test_that("omics_input works with rnaseq type", {
  mat <- matrix(1:20, nrow = 5, dimnames = list(paste0("g", 1:5), paste0("s", 1:4)))
  meta <- data.frame(group = rep(c("A", "B"), each = 2), row.names = paste0("s", 1:4))
  feat <- data.frame(feature_id = paste0("g", 1:5), name = paste0("Gene", 1:5),
                     stringsAsFactors = FALSE)
  inp <- omics_input(mat, meta, feat, omics_type = "rnaseq",
                     assay_type = "raw_count")
  expect_equal(inp$omics_type, "rnaseq")
  expect_equal(inp$assay_type, "raw_count")
})

test_that("omics_input stores raw_mat when provided", {
  mat <- matrix(1:20, nrow = 5, dimnames = list(paste0("g", 1:5), paste0("s", 1:4)))
  meta <- data.frame(group = rep(c("A", "B"), each = 2), row.names = paste0("s", 1:4))
  feat <- data.frame(feature_id = paste0("g", 1:5), name = paste0("Gene", 1:5),
                     stringsAsFactors = FALSE)
  raw <- matrix(1:20, nrow = 5, dimnames = list(paste0("g", 1:5), paste0("s", 1:4)))
  inp <- omics_input(mat, meta, feat, omics_type = "proteomics",
                     raw_mat = raw)
  expect_equal(inp$raw_mat, raw)
})

test_that("omics_input stores normalized_mat when provided", {
  mat <- matrix(1:20, nrow = 5, dimnames = list(paste0("g", 1:5), paste0("s", 1:4)))
  meta <- data.frame(group = rep(c("A", "B"), each = 2), row.names = paste0("s", 1:4))
  feat <- data.frame(feature_id = paste0("g", 1:5), name = paste0("Gene", 1:5),
                     stringsAsFactors = FALSE)
  norm <- mat / 10
  inp <- omics_input(mat, meta, feat, omics_type = "proteomics",
                     normalized_mat = norm)
  expect_equal(inp$normalized_mat, norm)
})

test_that("omics_input stores raw_object", {
  mat <- matrix(1:20, nrow = 5, dimnames = list(paste0("g", 1:5), paste0("s", 1:4)))
  meta <- data.frame(group = rep(c("A", "B"), each = 2), row.names = paste0("s", 1:4))
  feat <- data.frame(feature_id = paste0("g", 1:5), name = paste0("Gene", 1:5),
                     stringsAsFactors = FALSE)
  raw_obj <- list(type = "xlsx", path = "/path/to/file.xlsx")
  inp <- omics_input(mat, meta, feat, omics_type = "proteomics",
                     raw_object = raw_obj)
  expect_equal(inp$raw_object, raw_obj)
})

# ---- omics_input: validation errors ------------------------------------

test_that("omics_input errors when expr_mat has no colnames", {
  mat <- matrix(1:12, nrow = 3, ncol = 4)
  meta <- data.frame(group = 1:4, row.names = paste0("s", 1:4))
  feat <- data.frame(feature_id = paste0("g", 1:3))
  expect_error(omics_input(mat, meta, feat, omics_type = "proteomics",
                           assay_type = "intensity"), "sample column names")
})

test_that("omics_input errors when meta has no rownames", {
  mat <- matrix(1:12, nrow = 3, dimnames = list(paste0("g", 1:3), paste0("s", 1:4)))
  meta <- data.frame(group = 1:4)
  feat <- data.frame(feature_id = paste0("g", 1:3))
  expect_error(omics_input(mat, meta, feat, omics_type = "proteomics",
                           assay_type = "intensity"), "rownames")
})

test_that("omics_input errors when feature_df lacks feature_id column", {
  mat <- matrix(1:12, nrow = 3, dimnames = list(paste0("g", 1:3), paste0("s", 1:4)))
  meta <- data.frame(group = 1:4, row.names = paste0("s", 1:4))
  feat <- data.frame(id = paste0("g", 1:3))
  expect_error(omics_input(mat, meta, feat, omics_type = "proteomics",
                           assay_type = "intensity"), "feature_id")
})

test_that("omics_input errors when samples in expr_mat not in meta", {
  mat <- matrix(1:12, nrow = 3, dimnames = list(paste0("g", 1:3), paste0("s", 1:4)))
  meta <- data.frame(group = 1:3, row.names = paste0("s", 1:3))
  feat <- data.frame(feature_id = paste0("g", 1:3))
  expect_error(omics_input(mat, meta, feat, omics_type = "proteomics",
                           assay_type = "intensity"),
               "ncol\\(expr_mat\\).*nrow\\(meta_df\\)|Not all samples")
})

test_that("omics_input errors when features not in feature_df", {
  mat <- matrix(1:12, nrow = 3, dimnames = list(paste0("g", 1:3), paste0("s", 1:4)))
  meta <- data.frame(group = 1:4, row.names = paste0("s", 1:4))
  feat <- data.frame(feature_id = paste0("g", 1:2))
  expect_error(omics_input(mat, meta, feat, omics_type = "proteomics",
                           assay_type = "intensity"),
               "nrow\\(expr_mat\\).*nrow\\(feature_df\\)|Not all features")
})

test_that("omics_input allows custom omics_type strings but warns", {
  mat <- matrix(1:12, nrow = 3, dimnames = list(paste0("g", 1:3), paste0("s", 1:4)))
  meta <- data.frame(group = 1:4, row.names = paste0("s", 1:4))
  feat <- data.frame(feature_id = paste0("g", 1:3))
  # Non-standard omics_type strings still construct — omicsCore is a
  # general engine — but validate_omics_input() warns so typos surface.
  expect_warning(
    inp <- omics_input(mat, meta, feat, omics_type = "custom_type",
                       assay_type = "custom"),
    "Unrecognised `omics_type`"
  )
  expect_equal(inp$omics_type, "custom_type")
  expect_equal(inp$assay_type, "custom")
})

test_that("omics_input rejects a missing or empty omics_type", {
  mat <- matrix(1:12, nrow = 3, dimnames = list(paste0("g", 1:3), paste0("s", 1:4)))
  meta <- data.frame(group = 1:4, row.names = paste0("s", 1:4))
  feat <- data.frame(feature_id = paste0("g", 1:3))
  expect_error(omics_input(mat, meta, feat, omics_type = ""),
               "non-empty single string")
  expect_error(omics_input(mat, meta, feat, omics_type = NA_character_),
               "non-empty single string")
})

test_that("supported omics types validate without a warning", {
  mat <- matrix(1:12, nrow = 3, dimnames = list(paste0("g", 1:3), paste0("s", 1:4)))
  meta <- data.frame(group = 1:4, row.names = paste0("s", 1:4))
  feat <- data.frame(feature_id = paste0("g", 1:3))
  for (type in c("proteomics", "rnaseq")) {
    expect_no_warning(omics_input(mat, meta, feat, omics_type = type))
  }
})

test_that("omics_input keeps rownames matching column order", {
  mat <- matrix(1:12, nrow = 3, dimnames = list(paste0("g", 1:3), c("s2", "s4", "s1", "s3")))
  meta <- data.frame(group = 1:4, row.names = c("s1", "s2", "s3", "s4"))
  feat <- data.frame(feature_id = paste0("g", 1:3))
  inp <- omics_input(mat, meta, feat, omics_type = "proteomics")
  expect_equal(nrow(inp$meta_df), ncol(mat))
  expect_true(all(colnames(mat) %in% rownames(inp$meta_df)))
})

# ---- is_omics_input ----------------------------------------------------

test_that("is_omics_input returns TRUE for omics_input objects", {
  mat <- matrix(1:12, nrow = 3, dimnames = list(paste0("g", 1:3), paste0("s", 1:4)))
  meta <- data.frame(group = 1:4, row.names = paste0("s", 1:4))
  feat <- data.frame(feature_id = paste0("g", 1:3))
  inp <- omics_input(mat, meta, feat, omics_type = "proteomics",
                     assay_type = "intensity")
  expect_true(is_omics_input(inp))
})

test_that("is_omics_input returns FALSE for other objects", {
  expect_false(is_omics_input(list()))
  expect_false(is_omics_input(NULL))
  expect_false(is_omics_input(data.frame()))
})

# ---- new_omics_input ---------------------------------------------------

test_that("new_omics_input creates without validation", {
  mat <- matrix(1:12, nrow = 3, dimnames = list(paste0("g", 1:3), paste0("s", 1:4)))
  meta <- data.frame(group = 1:4, row.names = paste0("s", 1:4))
  feat <- data.frame(feature_id = paste0("g", 1:3))
  inp <- new_omics_input(omics_type = "proteomics", assay_type = "intensity",
                         expr_mat = mat, meta_df = meta, feature_df = feat)
  expect_s3_class(inp, "omics_input")
})

# ---- print.omics_input -------------------------------------------------

test_that("print.omics_input prints without error", {
  mat <- matrix(1:12, nrow = 3, dimnames = list(paste0("g", 1:3), paste0("s", 1:4)))
  meta <- data.frame(group = 1:4, row.names = paste0("s", 1:4))
  feat <- data.frame(feature_id = paste0("g", 1:3))
  inp <- omics_input(mat, meta, feat, omics_type = "proteomics",
                     assay_type = "intensity")
  expect_output(print(inp), "omics_input")
  expect_invisible(print(inp))
})
