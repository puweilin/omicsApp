# Values a matrix can carry that no engine was written for.
#
# An intensity matrix reaches run_diff() with -Inf where a log of zero
# was taken upstream; a counts matrix arrives with a sample that
# sequenced nothing, a value past the integer range because a column
# was summed twice, or a fraction because it was scaled. What the
# engines said before: an infinite effect with no p-value, "missing
# value where TRUE/FALSE needed", "every gene contains at least one
# zero", and an invalid-object error after a coercion warning. What
# they say now is below, and each names the value or the sample.

skip_if_not_installed("limma")

prot <- realistic_input("proteomics")
rna <- realistic_input("rnaseq")

group_diff <- function(input, method) {
  run_diff(input, method = method, analysis_type = "group",
           group_col = "group", control_group = "G1", case_group = "G2")
}

expect_sane_table <- function(bundle) {
  df <- bundle$results$diff_result_df
  p <- df$p_value[!is.na(df$p_value)]
  expect_true(all(is.finite(p)))
  expect_true(all(p >= 0 & p <= 1))
  expect_true(all(is.finite(df$effect) | is.na(df$effect)))
  invisible(df)
}

# ---- intensities --------------------------------------------------------

test_that("an infinite intensity is treated as missing, and the bundle says so", {
  i <- prot
  i$expr_mat[1, 1] <- Inf
  i$expr_mat[2, 2] <- -Inf
  for (m in c("ttest", "limma")) {
    expect_warning(group_diff(i, m), "2 infinite value\\(s\\) in 2 feature\\(s\\)")
    b <- suppressWarnings(group_diff(i, m))
    expect_match(b$warnings, "treated as missing")
    df <- expect_sane_table(b)
    # The features are still tested, on the values that remain.
    expect_false(is.na(df$p_value[df$feature_id == rownames(i$expr_mat)[1]]))
  }
})

test_that("a clean matrix carries no such warning", {
  b <- group_diff(prot, "ttest")
  expect_length(b$warnings, 0L)
})

test_that("NaN is missing already, and needs no warning", {
  i <- prot
  i$expr_mat[3, 3] <- NaN
  expect_no_warning(b <- group_diff(i, "ttest"))
  expect_sane_table(b)
})

test_that("a feature with no values, a constant feature and a zero row still yield a table", {
  i <- prot
  i$expr_mat[4, ] <- NA
  i$expr_mat[6, ] <- 7
  i$expr_mat[9, ] <- 0
  for (m in c("ttest", "limma")) {
    b <- suppressWarnings(group_diff(i, m))
    df <- expect_sane_table(b)
    expect_true(is.na(df$p_value[df$feature_id == rownames(i$expr_mat)[4]]))
    expect_identical(nrow(df), nrow(i$expr_mat))
  }
})

test_that("a sample with no values does not stop QC, the PCA or the test", {
  i <- prot
  i$expr_mat[, 5] <- NA
  expect_s3_class(suppressWarnings(run_qc(i)), "analysis_bundle")
  expect_no_error(ggplot2::ggplot_build(plot_pca(i, color_by = "group")))
  expect_sane_table(group_diff(i, "ttest"))
})

test_that("a value at the edge of double precision is carried, not crashed on", {
  i <- prot
  i$expr_mat[7, 1] <- 1e300
  expect_sane_table(suppressWarnings(group_diff(i, "limma")))
  expect_no_error(ggplot2::ggplot_build(plot_volcano(group_diff(i, "ttest"))))
})

test_that("a feature observed in one group only gets no estimate, not an error", {
  i <- prot
  i$expr_mat[10, i$meta_df$group == "G1"] <- NA
  for (m in c("ttest", "limma")) {
    df <- expect_sane_table(suppressWarnings(group_diff(i, m)))
    expect_true(is.na(df$effect[df$feature_id == rownames(i$expr_mat)[10]]))
  }
})

# ---- counts ---------------------------------------------------------------

count_engines <- c(if (requireNamespace("DESeq2", quietly = TRUE)) "deseq2",
                   if (requireNamespace("edgeR", quietly = TRUE)) "edger")

test_that("infinite, negative and missing counts are refused by count", {
  skip_if(length(count_engines) == 0L, "no count engine installed")
  cases <- list(
    infinite = function(i) { i$expr_mat[1, 1] <- Inf; i },
    negative = function(i) { i$expr_mat[2, 2] <- -5; i },
    missing  = function(i) { i$expr_mat[3, 3] <- NA; i }
  )
  for (m in count_engines) for (nm in names(cases)) {
    i <- cases[[nm]](rna)
    expect_error(group_diff(i, m), "Counts must be finite and non-negative", info = paste(m, nm))
    expect_error(group_diff(i, m), sprintf("1 %s", nm), info = paste(m, nm))
  }
})

test_that("a sample that sequenced nothing is named", {
  skip_if(length(count_engines) == 0L, "no count engine installed")
  i <- rna
  i$expr_mat[, 6] <- 0L
  for (m in count_engines) {
    expect_error(group_diff(i, m),
                 sprintf("1 sample\\(s\\) have no counts at all: %s", colnames(i$expr_mat)[6]))
  }
})

test_that("an empty sample outside the contrast does not block it", {
  skip_if(length(count_engines) == 0L, "no count engine installed")
  i <- rna
  i$meta_df$group <- as.character(i$meta_df$group)
  i$meta_df$group[6] <- "G3"
  i$expr_mat[, 6] <- 0L
  for (m in count_engines) expect_sane_table(suppressMessages(group_diff(i, m)))
})

test_that("a count past the integer range is refused for DESeq2 and carried by edgeR", {
  i <- rna
  i$expr_mat[5, 5] <- 3e9
  if ("deseq2" %in% count_engines) {
    expect_error(group_diff(i, "deseq2"), "exceed 2147483647")
  }
  if ("edger" %in% count_engines) {
    expect_sane_table(group_diff(i, "edger"))
  }
})

test_that("a fractional count is refused by DESeq2 in words and accepted by edgeR", {
  i <- rna
  i$expr_mat[4, 4] <- 2.5
  if ("deseq2" %in% count_engines) {
    expect_error(group_diff(i, "deseq2"), "integer-like")
  }
  if ("edger" %in% count_engines) {
    expect_sane_table(group_diff(i, "edger"))
  }
})

test_that("a gene with no counts, or none in one group, gets a row and no error", {
  skip_if(length(count_engines) == 0L, "no count engine installed")
  i <- rna
  i$expr_mat[7, ] <- 0L
  i$expr_mat[8, i$meta_df$group == "G1"] <- 0L
  for (m in count_engines) {
    df <- expect_sane_table(suppressMessages(group_diff(i, m)))
    expect_identical(nrow(df), nrow(i$expr_mat))
  }
})

test_that("QC and depth still run on every one of these matrices", {
  i <- rna
  i$expr_mat[1, 1] <- Inf
  i$expr_mat[, 6] <- 0L
  i$expr_mat[5, 5] <- 3e9
  expect_s3_class(suppressWarnings(run_qc(i)), "analysis_bundle")
  d <- qc_depth(i)
  expect_identical(nrow(d), ncol(i$expr_mat))
  expect_true(colnames(i$expr_mat)[6] %in% qc_depth_outliers(d))
})
