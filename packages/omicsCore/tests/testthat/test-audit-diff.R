# Differential-analysis defects found in the 2026-10 audit, one test each.

ad_counts <- function(n_genes = 300L, per_group = 4L, groups = c("A", "B"),
                      lib = NULL, seed = 1L) {
  set.seed(seed)
  g <- rep(groups, each = per_group)
  ids <- paste0("S", seq_along(g))
  mu <- exp(stats::rnorm(n_genes, log(200), 1))
  m <- matrix(stats::rnbinom(n_genes * length(g), mu = mu, size = 20), n_genes,
              dimnames = list(paste0("G", seq_len(n_genes)), ids))
  if (!is.null(lib)) m <- round(sweep(m, 2L, lib, "*"))
  storage.mode(m) <- "integer"
  omics_input(m, data.frame(group = g, row.names = ids, stringsAsFactors = FALSE),
              data.frame(feature_id = rownames(m), feature_symbol = rownames(m),
                         row.names = rownames(m)),
              omics_type = "rnaseq", assay_type = "raw_count")
}

ad_intensity <- function(log = FALSE, seed = 2L) {
  set.seed(seed)
  g <- rep(c("A", "B"), each = 5)
  ids <- paste0("S", 1:10)
  m <- matrix(stats::rnorm(100 * 10, 20, 0.3), 100, dimnames = list(paste0("P", 1:100), ids))
  m[1:20, g == "B"] <- m[1:20, g == "B"] + 1
  if (!log) m <- 2^m
  omics_input(m, data.frame(group = g, row.names = ids, stringsAsFactors = FALSE),
              data.frame(feature_id = rownames(m), row.names = rownames(m)),
              omics_type = "proteomics",
              assay_type = if (log) "normalized_intensity" else "raw_intensity")
}

test_that("linear intensities are log2-transformed before limma, so effects are log2 fold changes", {
  inp <- suppressWarnings(ad_intensity(log = FALSE))
  b <- suppressWarnings(run_diff(inp, group_col = "group", control_group = "A",
                                 case_group = "B"))
  expect_identical(b$params$method, "limma")
  expect_equal(stats::median(b$results$diff_result_df$effect[1:20]), 1, tolerance = 0.15)
  expect_match(b$warnings, "log2-transformed", all = FALSE)
  expect_identical(b$input_info$assay_type, "raw_intensity")
  expect_identical(b$input_info$analysed_scale, "normalized_intensity")
})

test_that("the t-test on counts normalizes library size, so deeper libraries are not 'up'", {
  inp <- ad_counts(per_group = 5L, lib = rep(c(1, 2), each = 5))
  b <- suppressWarnings(run_diff(inp, method = "ttest", group_col = "group",
                                 control_group = "A", case_group = "B"))
  hits <- sum(b$results$diff_result_df$adj_p_value < 0.05, na.rm = TRUE)
  expect_lt(hits, 15L)   # was 347 of 500 before
  expect_match(b$warnings, "log2-CPM", all = FALSE)
})

test_that("edgeR sets aside genes too low to test, and says how many", {
  skip_if_not_installed("edgeR")
  inp <- ad_counts()
  m <- inp$expr_mat
  m[1:100, ] <- 0L
  m[1:100, 1] <- 1L
  inp$expr_mat <- m
  b <- run_diff(inp, method = "edger", group_col = "group",
                control_group = "A", case_group = "B")
  df <- b$results$diff_result_df
  expect_equal(nrow(df), 300L)
  expect_true(all(is.na(df$p_value[1:100])))
  expect_match(b$warnings, "filterByExpr", all = FALSE)
})

test_that("the edgeR global test does not test a covariate whose name starts with the group column", {
  skip_if_not_installed("edgeR")
  inp <- ad_counts(groups = c("A", "B", "C"), per_group = 4L)
  inp$meta_df$group_batch <- rep(c("x", "y"), 6)
  one <- run_diff(inp, method = "edger", analysis_type = "anova",
                  group_col = "group", covariates = "group_batch")
  inp$meta_df$batch2 <- inp$meta_df$group_batch
  two <- run_diff(inp, method = "edger", analysis_type = "anova",
                  group_col = "group", covariates = "batch2")
  expect_equal(one$results$diff_result_df$p_value, two$results$diff_result_df$p_value)
})

test_that("the limma spline fit is reachable", {
  set.seed(5)
  ids <- paste0("S", 1:16)
  age <- seq(20, 80, length.out = 16)
  m <- matrix(stats::rnorm(50 * 16, 10), 50, dimnames = list(paste0("P", 1:50), ids))
  m[1:5, ] <- m[1:5, ] + rep(((age - 50) / 15)^2, each = 5)
  inp <- omics_input(m, data.frame(age = age, row.names = ids),
                     data.frame(feature_id = rownames(m), row.names = rownames(m)),
                     omics_type = "proteomics", assay_type = "normalized_intensity")
  b <- run_diff_continuous(inp, method = "limma", continuous_col = "age",
                           model = "spline", df = 3)
  df <- b$results$diff_result_df
  expect_identical(unique(df$analysis_type), "continuous_spline")
  expect_true(all(df$adj_p_value[1:5] < 0.05))
})

test_that("a factor covariate level used only outside the comparison is not 'confounded'", {
  inp <- ad_intensity(log = TRUE)
  inp$expr_mat <- cbind(inp$expr_mat, inp$expr_mat[, 1:4] + 0.1)
  colnames(inp$expr_mat)[11:14] <- paste0("C", 1:4)
  meta <- inp$meta_df
  meta <- rbind(meta, data.frame(group = rep("C", 4), row.names = paste0("C", 1:4)))
  meta$batch <- factor(c(rep(c("b1", "b2"), 5), rep("b3", 4)))
  inp$meta_df <- meta
  b <- run_diff(inp, method = "limma", group_col = "group", control_group = "A",
                case_group = "B", covariates = "batch")
  expect_identical(b$params$covariates, "batch")
})

test_that("the global test checks its design like the comparisons do", {
  inp <- ad_intensity(log = TRUE)
  inp$meta_df$group <- rep(c("A", "B", "C", "A", "B"), 2)
  inp$meta_df$age <- c(NA, 31:39)
  expect_error(run_diff(inp, method = "limma", analysis_type = "anova",
                        group_col = "group", covariates = "age"),
               "missing values")
  expect_error(run_diff(inp, method = "limma", analysis_type = "anova",
                        group_col = "group", selected_groups = c("A", "B", "Z")),
               "selected_groups")
})

test_that("a covariate the pairing already absorbs is refused clearly for DESeq2", {
  skip_if_not_installed("DESeq2")
  inp <- ad_counts(per_group = 4L)
  inp$meta_df$subject <- rep(paste0("s", 1:4), 2)
  inp$meta_df$sex <- rep(c("F", "M", "F", "M"), 2)
  expect_error(run_diff(inp, method = "deseq2", group_col = "group",
                        control_group = "A", case_group = "B",
                        paired_col = "subject", covariates = "sex"),
               "does not vary within the blocks")
})

test_that("numeric pair ids are a block, not a slope, in the limma continuous fit", {
  set.seed(8)
  ids <- paste0("S", 1:12)
  m <- matrix(stats::rnorm(20 * 12, 10), 20, dimnames = list(paste0("P", 1:20), ids))
  meta <- data.frame(age = rep(c(30, 50), 6), subj = rep(1:6, each = 2),
                     `Body mass` = stats::rnorm(12, 70, 5),
                     row.names = ids, check.names = FALSE)
  inp <- omics_input(m, meta, data.frame(feature_id = rownames(m), row.names = rownames(m)),
                     omics_type = "proteomics", assay_type = "normalized_intensity")
  b <- run_diff_continuous(inp, method = "limma", continuous_col = "age",
                           paired_col = "subj", covariates = "Body mass")
  y <- m[1, ]
  ref <- summary(stats::lm(y ~ meta$age + factor(meta$subj) + meta$`Body mass`))$adj.r.squared
  expect_equal(b$results$diff_result_df$model_fit[1], ref, tolerance = 1e-8)
})

test_that("a group column that is not an R name works in every engine", {
  inp <- ad_counts()
  names(inp$meta_df) <- "Treatment Group"
  for (m in c("edger", "deseq2")) {
    skip_if_not_installed(if (m == "edger") "edgeR" else "DESeq2")
    b <- suppressWarnings(suppressMessages(
      run_diff(inp, method = m, group_col = "Treatment Group",
               control_group = "A", case_group = "B")))
    expect_identical(b$params$comparison, "B_vs_A", info = m)
  }
  li <- ad_intensity(log = TRUE)
  names(li$meta_df) <- "Treatment Group"
  li$meta_df$`Body mass` <- stats::rnorm(10, 70, 3)
  b <- run_diff(li, method = "lm", group_col = "Treatment Group",
                control_group = "A", case_group = "B", covariates = "Body mass")
  expect_true(all(b$results$diff_result_df$effect[1:20] > 0.5))
})

test_that("figures refuse a bundle that holds several comparisons", {
  inp <- ad_intensity(log = TRUE)
  inp$meta_df$group <- rep(c("A", "B", "C", "A", "B"), 2)
  b <- run_diff(inp, method = "limma", group_col = "group", control_group = "A",
                case_group = c("B", "C"))
  expect_error(plot_volcano(b), "select_comparison")
  expect_error(plot_ma(b), "select_comparison")
  expect_s3_class(plot_volcano(select_comparison(b, "B_vs_A")), "ggplot")
})
