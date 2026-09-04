# Edge-case tests for subset_omics, subset_omics_samples,
# subset_omics_features, drop_meta_na, and select_complete_cases.

make_subset_input <- function() {
  mat <- matrix(rnorm(60, mean = 10, sd = 2), nrow = 10, ncol = 6)
  rownames(mat) <- paste0("gene_", 1:10)
  colnames(mat) <- paste0("s", 1:6)
  meta <- data.frame(group = c("A", "A", "A", "B", "B", "B"),
                     row.names = paste0("s", 1:6))
  feat <- data.frame(feature_id = paste0("gene_", 1:10),
                     feature_name = paste0("Gene", 1:10),
                     stringsAsFactors = FALSE)
  omics_input(mat, meta, feat, omics_type = "proteomics",
              assay_type = "normalized_intensity")
}

# ---- subset_omics ------------------------------------------------------

test_that("subset_omics subsets both axes at once", {
  inp <- make_subset_input()
  samps <- colnames(inp$expr_mat)[1:3]
  feats <- rownames(inp$expr_mat)[1:5]
  sub <- subset_omics(inp, samples = samps, features = feats)
  expect_equal(colnames(sub$expr_mat), samps)
  expect_equal(rownames(sub$expr_mat), feats)
  expect_equal(nrow(sub$meta_df), 3)
  expect_equal(nrow(sub$feature_df), 5)
})

test_that("subset_omics with NULL samples retains all samples", {
  inp <- make_subset_input()
  sub <- subset_omics(inp, samples = NULL, features = rownames(inp$expr_mat)[1:3])
  expect_equal(ncol(sub$expr_mat), ncol(inp$expr_mat))
  expect_equal(nrow(sub$expr_mat), 3)
})

test_that("subset_omics with NULL features retains all features", {
  inp <- make_subset_input()
  sub <- subset_omics(inp, samples = colnames(inp$expr_mat)[1:3], features = NULL)
  expect_equal(nrow(sub$expr_mat), nrow(inp$expr_mat))
  expect_equal(ncol(sub$expr_mat), 3)
})

test_that("subset_omics with both NULL returns identical copy", {
  inp <- make_subset_input()
  sub <- subset_omics(inp, samples = NULL, features = NULL)
  expect_equal(sub$expr_mat, inp$expr_mat)
  expect_equal(nrow(sub$meta_df), nrow(inp$meta_df))
})

test_that("subset_omics validates input first", {
  expect_error(subset_omics(list(), samples = "x"), "omics_input")
})

# ---- subset_omics_samples ----------------------------------------------

test_that("subset_omics_samples errors when no samples match", {
  inp <- make_subset_input()
  expect_error(subset_omics_samples(inp, c("nonexistent_1", "nonexistent_2")),
               "No requested sample")
})

test_that("subset_omics_samples intersects with existing samples", {
  inp <- make_subset_input()
  existing <- colnames(inp$expr_mat)[1:2]
  sub <- subset_omics_samples(inp, c(existing, "fake_sample"))
  expect_equal(colnames(sub$expr_mat), existing)
  expect_equal(ncol(sub$expr_mat), 2)
})

test_that("subset_omics_samples preserves raw_mat when present", {
  inp <- make_subset_input()
  inp$raw_mat <- inp$expr_mat + 1
  rownames(inp$raw_mat) <- rownames(inp$expr_mat)
  colnames(inp$raw_mat) <- colnames(inp$expr_mat)
  samp <- colnames(inp$expr_mat)[1:2]
  sub <- subset_omics_samples(inp, samp)
  expect_equal(colnames(sub$raw_mat), samp)
  expect_equal(sub$raw_mat[1, 1], inp$raw_mat[1, 1])
})

test_that("subset_omics_samples preserves normalized_mat when present", {
  inp <- make_subset_input()
  inp$normalized_mat <- inp$expr_mat * 2
  rownames(inp$normalized_mat) <- rownames(inp$expr_mat)
  colnames(inp$normalized_mat) <- colnames(inp$expr_mat)
  samp <- colnames(inp$expr_mat)[2:4]
  sub <- subset_omics_samples(inp, samp)
  expect_equal(colnames(sub$normalized_mat), samp)
})

test_that("subset_omics_samples preserves feature_df unchanged", {
  inp <- make_subset_input()
  samp <- colnames(inp$expr_mat)[1:3]
  sub <- subset_omics_samples(inp, samp)
  expect_equal(nrow(sub$feature_df), nrow(inp$feature_df))
})

# ---- subset_omics_features ---------------------------------------------

test_that("subset_omics_features errors when no features match", {
  inp <- make_subset_input()
  expect_error(subset_omics_features(inp, c("nonexistent1", "nonexistent2")),
               "No requested feature")
})

test_that("subset_omics_features aligns feature_df rows", {
  inp <- make_subset_input()
  feats <- rownames(inp$expr_mat)[2:4]
  sub <- subset_omics_features(inp, feats)
  expect_equal(rownames(sub$expr_mat), feats)
  expect_equal(sub$feature_df$feature_id, feats)
})

test_that("subset_omics_features intersects with existing features", {
  inp <- make_subset_input()
  existing <- rownames(inp$expr_mat)[1:3]
  sub <- subset_omics_features(inp, c(existing, "fake"))
  expect_equal(rownames(sub$expr_mat), existing)
  expect_equal(nrow(sub$expr_mat), 3)
})

test_that("subset_omics_features preserves sample metadata unchanged", {
  inp <- make_subset_input()
  feats <- rownames(inp$expr_mat)[1:5]
  sub <- subset_omics_features(inp, feats)
  expect_equal(nrow(sub$meta_df), nrow(inp$meta_df))
})

# ---- drop_meta_na ------------------------------------------------------

test_that("drop_meta_na removes samples with NA in specified columns", {
  inp <- make_subset_input()
  inp$meta_df$group[1] <- NA
  inp$meta_df$group[3] <- NA
  sub <- drop_meta_na(inp, "group")
  expect_lt(nrow(sub$meta_df), nrow(inp$meta_df))
  expect_false(anyNA(sub$meta_df$group))
})

test_that("drop_meta_na errors when column missing", {
  inp <- make_subset_input()
  expect_error(drop_meta_na(inp, "nonexistent_col"), "required columns")
})

test_that("drop_meta_na on complete data retains all samples", {
  inp <- make_subset_input()
  sub <- drop_meta_na(inp, "group")
  expect_equal(nrow(sub$meta_df), nrow(inp$meta_df))
})

test_that("drop_meta_na respects multiple columns", {
  inp <- make_subset_input()
  inp$meta_df$extra <- "x"
  inp$meta_df$extra[c(2, 5)] <- NA
  sub <- drop_meta_na(inp, c("group", "extra"))
  expect_equal(nrow(sub$meta_df), sum(stats::complete.cases(inp$meta_df[, c("group", "extra")])))
})

# ---- select_complete_cases ----------------------------------------------

test_that("select_complete_cases filters features by missingness", {
  inp <- make_subset_input()
  n_before <- nrow(inp$expr_mat)
  # inject NAs into some features
  inp$expr_mat[1, 1:3] <- NA_real_   # 3/6 = 0.5 missing
  inp$expr_mat[2, 1:6] <- NA_real_   # all missing
  sub <- select_complete_cases(inp, feature_missing_cutoff = 0.4)
  expect_lt(nrow(sub$expr_mat), n_before)
  expect_false(rownames(inp$expr_mat)[1] %in% rownames(sub$expr_mat))
  expect_false(rownames(inp$expr_mat)[2] %in% rownames(sub$expr_mat))
})

test_that("select_complete_cases with cutoff=1 keeps all features", {
  inp <- make_subset_input()
  sub <- select_complete_cases(inp, feature_missing_cutoff = 1)
  expect_equal(nrow(sub$expr_mat), nrow(inp$expr_mat))
})

test_that("select_complete_cases with cutoff=0 keeps only complete features", {
  inp <- make_subset_input()
  inp$expr_mat[1, 1] <- NA_real_
  sub <- select_complete_cases(inp, feature_missing_cutoff = 0)
  expect_false(anyNA(sub$expr_mat))
})

test_that("select_complete_cases with cutoff=0 removes all features if all incomplete", {
  inp <- make_subset_input()
  # make every feature have at least one NA
  for (i in seq_len(nrow(inp$expr_mat))) {
    inp$expr_mat[i, 1] <- NA_real_
  }
  # Implementation errors when zero features survive the filter
  expect_error(select_complete_cases(inp, feature_missing_cutoff = 0),
               "No requested feature")
})
