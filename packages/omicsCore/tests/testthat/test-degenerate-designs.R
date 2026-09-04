# Designs that cannot answer the question, refused before an engine
# turns them into a table.
#
# A probe of what each engine did with these found the same pattern
# again: limma and lm dropped a covariate confounded with the group and
# returned the *unadjusted* result under the adjusted label; lm returned
# a table of NA for a constant covariate; a continuous variable with no
# variation produced all-NA tables from both; a case level absent from
# the column reached limma as "trying to take contrast of non-estimable
# coefficient". Only DESeq2 and edgeR refused. run_diff() now refuses
# for every engine, naming the column.

base <- function() realistic_input(n_per_group = 4L)

group_run <- function(inp, method = "limma", control = "G1", case = "G2", ...) {
  run_diff(inp, method = method, analysis_type = "group", group_col = "group",
           control_group = control, case_group = case, ...)
}

engines <- function() c("ttest", "lm", if (requireNamespace("limma", quietly = TRUE)) "limma")

test_that("a covariate identical to the group is refused, not silently dropped", {
  inp <- base()
  inp$meta_df$batch <- inp$meta_df$group
  for (m in setdiff(engines(), "ttest")) {
    expect_error(group_run(inp, m, covariates = "batch"), "confounded with `group`",
                 info = m)
  }
  # And through the count engines, with the same words
  skip_if_not_installed("DESeq2")
  rna <- realistic_input("rnaseq", n_per_group = 3L)
  rna$meta_df$batch <- rna$meta_df$group
  expect_error(group_run(rna, "deseq2", covariates = "batch"), "confounded with `group`")
})

test_that("a covariate that merely correlates with the group is allowed", {
  inp <- base()
  # Balanced batch: two of each group in each batch
  inp$meta_df$batch <- rep(c("A", "B"), times = 4L)
  for (m in setdiff(engines(), "ttest")) {
    b <- group_run(inp, m, covariates = "batch")
    expect_gt(nrow(b$results$diff_result_df), 0L)
    expect_false(anyNA(b$results$diff_result_df$effect))
  }
})

test_that("a constant covariate is refused", {
  inp <- base()
  inp$meta_df$site <- "one site"
  for (m in setdiff(engines(), "ttest")) {
    expect_error(group_run(inp, m, covariates = "site"), "constant", info = m)
  }
})

test_that("a covariate with missing values names the samples", {
  inp <- base()
  inp$meta_df$age[c(2L, 5L)] <- NA
  for (m in setdiff(engines(), "ttest")) {
    expect_error(group_run(inp, m, covariates = "age"), "missing values in 2", info = m)
    expect_error(group_run(inp, m, covariates = "age"), rownames(inp$meta_df)[2L],
                 fixed = TRUE, info = m)
  }
})

test_that("a missing value outside the contrast does not block it", {
  inp <- base()
  inp$meta_df$group[1L] <- "G3"
  inp$meta_df$age[1L] <- NA   # the sample that is not in G1 vs G2
  for (m in setdiff(engines(), "ttest")) {
    b <- group_run(inp, m, covariates = "age")
    expect_gt(nrow(b$results$diff_result_df), 0L)
  }
})

test_that("a group level that is not in the column is named, with what is", {
  inp <- base()
  for (m in engines()) {
    err <- tryCatch(group_run(inp, m, case = "G9"),
                    error = function(e) conditionMessage(e))
    expect_match(err, "'G9' is not a level of `group`", fixed = TRUE, info = m)
    expect_match(err, "'G1', 'G2'", fixed = TRUE, info = m)
  }
})

test_that("control and case must differ", {
  inp <- base()
  expect_error(group_run(inp, "ttest", control = "G1", case = "G1"), "distinct")
})

test_that("a continuous variable with no variation is refused", {
  inp <- base()
  inp$meta_df$age <- 40
  for (m in setdiff(engines(), "ttest")) {
    expect_error(run_diff(inp, method = m, analysis_type = "continuous",
                          continuous_col = "age"),
                 "no variation", info = m)
  }
})

test_that("a continuous variable with missing values names the samples", {
  inp <- base()
  inp$meta_df$age[3L] <- NA
  for (m in setdiff(engines(), "ttest")) {
    expect_error(run_diff(inp, method = m, analysis_type = "continuous",
                          continuous_col = "age"),
                 rownames(inp$meta_df)[3L], fixed = TRUE, info = m)
  }
})

test_that("a covariate confounded with the continuous variable is refused", {
  inp <- base()
  inp$meta_df$age2 <- inp$meta_df$age * 2
  for (m in setdiff(engines(), "ttest")) {
    expect_error(run_diff(inp, method = m, analysis_type = "continuous",
                          continuous_col = "age", covariates = "age2"),
                 "confounded with `age`", info = m)
  }
})

test_that("the messages the backends already used are unchanged", {
  inp <- base()
  expect_error(run_diff(inp, method = "ttest", analysis_type = "group",
                        group_col = "nope", control_group = "G1", case_group = "G2"),
               "`group_col` not found in `meta_df`: nope", fixed = TRUE)
  expect_error(group_run(inp, "lm", covariates = "nope"),
               "Missing covariates: nope", fixed = TRUE)
  expect_error(run_diff(inp, method = "lm", analysis_type = "continuous",
                        continuous_col = "nope"),
               "`continuous_col` not found", fixed = TRUE)
})

test_that("a valid design still runs on every engine", {
  inp <- base()
  for (m in engines()) {
    b <- group_run(inp, m)
    expect_identical(nrow(b$results$diff_result_df), nrow(inp$expr_mat))
  }
  for (m in setdiff(engines(), "ttest")) {
    b <- run_diff(inp, method = m, analysis_type = "continuous",
                  continuous_col = "age", covariates = "sex")
    expect_identical(nrow(b$results$diff_result_df), nrow(inp$expr_mat))
  }
})
