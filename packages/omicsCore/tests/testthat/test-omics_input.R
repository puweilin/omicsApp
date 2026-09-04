make_demo_input <- function(
  n_features = 6, n_samples = 5,
  omics_type = "proteomics"
) {
  set.seed(42)
  feat_ids <- paste0("g", seq_len(n_features))
  samp_ids <- paste0("s", seq_len(n_samples))
  expr <- matrix(
    rnorm(n_features * n_samples),
    nrow = n_features,
    dimnames = list(feat_ids, samp_ids)
  )
  meta <- data.frame(
    group = rep(c("A", "B"), length.out = n_samples),
    age = seq(20, by = 5, length.out = n_samples),
    row.names = samp_ids,
    stringsAsFactors = FALSE
  )
  feat <- data.frame(
    feature_id = feat_ids,
    feature_symbol = feat_ids,
    row.names = feat_ids,
    stringsAsFactors = FALSE
  )
  omics_input(expr, meta, feat, omics_type = omics_type,
              assay_type = "normalized_intensity")
}

test_that("omics_input() constructs and validates", {
  x <- make_demo_input()
  expect_true(is_omics_input(x))
  expect_silent(validate_omics_input(x))
  expect_equal(nrow(x$expr_mat), 6)
  expect_equal(ncol(x$expr_mat), 5)
})

test_that("omics_input() rejects mismatched dimensions", {
  set.seed(1)
  expr <- matrix(rnorm(20), nrow = 4,
                 dimnames = list(paste0("g", 1:4), paste0("s", 1:5)))
  meta <- data.frame(group = c("A", "B"),
                     row.names = c("s1", "s2"))
  feat <- data.frame(feature_id = paste0("g", 1:4),
                     row.names = paste0("g", 1:4))
  expect_error(
    omics_input(expr, meta, feat, omics_type = "proteomics",
                assay_type = "normalized_intensity"),
    "ncol\\(expr_mat\\)"
  )
})

test_that("omics_input() requires feature_id column", {
  set.seed(1)
  expr <- matrix(rnorm(20), nrow = 4,
                 dimnames = list(paste0("g", 1:4), paste0("s", 1:5)))
  meta <- data.frame(group = rep("A", 5),
                     row.names = paste0("s", 1:5))
  feat <- data.frame(symbol = paste0("g", 1:4),
                     row.names = paste0("g", 1:4))
  expect_error(
    omics_input(expr, meta, feat, omics_type = "proteomics",
                assay_type = "normalized_intensity"),
    "feature_id"
  )
})

test_that("subset_omics() handles samples, features, or both", {
  x <- make_demo_input()

  by_samples <- subset_omics(x, samples = c("s1", "s3"))
  expect_equal(ncol(by_samples$expr_mat), 2)
  expect_equal(colnames(by_samples$expr_mat), c("s1", "s3"))

  by_features <- subset_omics(x, features = c("g2", "g4"))
  expect_equal(nrow(by_features$expr_mat), 2)
  expect_equal(rownames(by_features$expr_mat), c("g2", "g4"))

  both <- subset_omics(x, samples = "s2", features = c("g1", "g3"))
  expect_equal(dim(both$expr_mat), c(2L, 1L))
})

test_that("subset_omics() with both NULL returns equivalent input", {
  x <- make_demo_input()
  y <- subset_omics(x)
  expect_equal(dim(y$expr_mat), dim(x$expr_mat))
})

test_that("subset_omics_samples / features error on no overlap", {
  x <- make_demo_input()
  expect_error(subset_omics_samples(x, "missing"), "No requested sample")
  expect_error(subset_omics_features(x, "missing"), "No requested feature")
})

test_that("drop_meta_na removes rows with NA in target column", {
  x <- make_demo_input()
  x$meta_df$batch <- c("b1", NA, "b1", "b2", NA)
  y <- drop_meta_na(x, "batch")
  expect_equal(ncol(y$expr_mat), 3)
  expect_false(any(is.na(y$meta_df$batch)))
})

test_that("select_complete_cases applies missingness cutoff", {
  x <- make_demo_input()
  x$expr_mat[1, ] <- NA          # 100% missing
  x$expr_mat[2, 1:3] <- NA       # 60% missing
  y <- select_complete_cases(x, feature_missing_cutoff = 0.5)
  expect_equal(nrow(y$expr_mat), 4)
  expect_false("g1" %in% rownames(y$expr_mat))
  expect_false("g2" %in% rownames(y$expr_mat))
})
