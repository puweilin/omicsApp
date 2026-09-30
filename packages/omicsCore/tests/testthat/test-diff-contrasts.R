# Several treatment groups against one control.

mc_input <- function(levels = c("ctrl", "A", "B"), n = 4L, omics = "proteomics") {
  set.seed(7)
  groups <- rep(levels, each = n)
  samp <- paste0("S", seq_along(groups))
  n_feat <- 40L
  ids <- paste0("F", seq_len(n_feat))
  if (omics == "proteomics") {
    m <- matrix(rnorm(n_feat * length(groups), 10, 0.4), n_feat,
                dimnames = list(ids, samp))
    m[1:5, groups == levels[2]] <- m[1:5, groups == levels[2]] + 3
    m[6:10, groups == levels[3]] <- m[6:10, groups == levels[3]] - 3
    assay <- "normalized_intensity"
  } else {
    m <- matrix(rpois(n_feat * length(groups), 200), n_feat,
                dimnames = list(ids, samp))
    m[1:5, groups == levels[2]] <- m[1:5, groups == levels[2]] * 6L
    m[6:10, groups == levels[3]] <- pmax(m[6:10, groups == levels[3]] %/% 6L, 1L)
    assay <- "raw_count"
  }
  omics_input(m,
              data.frame(group = groups, donor = rep(paste0("D", 1:n), length(levels)),
                         row.names = samp, stringsAsFactors = FALSE),
              data.frame(feature_id = ids, feature_symbol = paste0("G", seq_len(n_feat)),
                         row.names = ids, stringsAsFactors = FALSE),
              omics_type = omics, assay_type = assay)
}

test_that("several case groups come back stacked, one comparison each", {
  inp <- mc_input()
  b <- run_diff(inp, method = "limma", group_col = "group",
                control_group = "ctrl", case_group = c("A", "B"))
  expect_identical(diff_comparisons(b), c("A_vs_ctrl", "B_vs_ctrl"))
  df <- b$results$diff_result_df
  check_diff_result_schema(df)
  expect_equal(nrow(df), 80L)
  expect_identical(b$params$case_group, c("A", "B"))
})

test_that("a single case group gives the same result as before", {
  inp <- mc_input()
  one <- run_diff(inp, method = "limma", group_col = "group",
                  control_group = "ctrl", case_group = "A")
  expect_identical(one$params$comparison, "A_vs_ctrl")
  expect_false("comparison" %in% names(one$results$diff_raw_df))
  # Two-group subset: identical to fitting only those samples.
  sub <- subset_omics_samples(inp, rownames(inp$meta_df)[inp$meta_df$group != "B"])
  alone <- run_diff(sub, method = "limma", group_col = "group",
                    control_group = "ctrl", case_group = "A")
  expect_equal(one$results$diff_result_df$p_value,
               alone$results$diff_result_df$p_value)
})

test_that("the shared limma fit pools variance across all groups", {
  inp <- mc_input()
  both <- run_diff(inp, method = "limma", group_col = "group",
                   control_group = "ctrl", case_group = c("A", "B"))
  one <- run_diff(inp, method = "limma", group_col = "group",
                  control_group = "ctrl", case_group = "A")
  a_shared <- select_comparison(both, "A_vs_ctrl")
  # Same log fold changes, different (pooled) standard errors.
  expect_equal(a_shared$results$diff_result_df$effect,
               one$results$diff_result_df$effect)
  expect_false(isTRUE(all.equal(a_shared$results$diff_result_df$p_value,
                                one$results$diff_result_df$p_value)))
  expect_gt(both$results$diff_object$df.residual[1],
            one$results$diff_object$df.residual[1])
})

test_that("group labels that are not R names work in limma", {
  inp <- mc_input(levels = c("ctrl", "Drug A", "24h"))
  b <- run_diff(inp, method = "limma", group_col = "group",
                control_group = "ctrl", case_group = c("Drug A", "24h"))
  expect_identical(diff_comparisons(b), c("Drug A_vs_ctrl", "24h_vs_ctrl"))
  hits <- select_comparison(b, "Drug A_vs_ctrl")$results$diff_result_df
  expect_true(all(hits$effect[1:5] > 2))
  single <- run_diff(inp, method = "limma", group_col = "group",
                     control_group = "ctrl", case_group = "24h")
  expect_true(all(single$results$diff_result_df$effect[6:10] < -2))
})

test_that("t-test and lm run one pair at a time and stack the results", {
  inp <- mc_input()
  for (m in c("ttest", "lm")) {
    b <- run_diff(inp, method = m, group_col = "group",
                  control_group = "ctrl", case_group = c("A", "B"))
    expect_identical(diff_comparisons(b), c("A_vs_ctrl", "B_vs_ctrl"))
    s <- select_comparison(b, "B_vs_ctrl")
    expect_equal(nrow(s$results$diff_result_df), 40L)
  }
})

test_that("edgeR and DESeq2 fit once and read every contrast off it", {
  skip_if_not_installed("edgeR")
  skip_if_not_installed("DESeq2")
  inp <- mc_input(omics = "rnaseq", levels = c("ctrl", "A", "AB"))
  for (m in c("edger", "deseq2")) {
    b <- suppressWarnings(suppressMessages(
      run_diff(inp, method = m, group_col = "group",
               control_group = "ctrl", case_group = c("A", "AB"))))
    expect_identical(diff_comparisons(b), c("A_vs_ctrl", "AB_vs_ctrl"))
    a <- select_comparison(b, "A_vs_ctrl")$results$diff_result_df
    ab <- select_comparison(b, "AB_vs_ctrl")$results$diff_result_df
    # "AB" must not be mistaken for "A" (edgeR used a substring match).
    expect_true(all(a$effect[1:5] > 1.5))
    expect_true(all(ab$effect[6:10] < -1.5))
    expect_true(all(abs(ab$effect[1:5]) < 1.5))
  }
})

test_that("select_comparison returns an ordinary single-contrast bundle", {
  inp <- mc_input()
  b <- run_diff(inp, method = "limma", group_col = "group",
                control_group = "ctrl", case_group = c("A", "B"))
  s <- select_comparison(b, "B_vs_ctrl")
  expect_identical(s$params$comparison, "B_vs_ctrl")
  expect_identical(s$params$case_group, "B")
  expect_identical(s$params$all_case_groups, c("A", "B"))
  expect_identical(unique(s$results$diff_result_df$comparison), "B_vs_ctrl")
  expect_false("comparison" %in% names(s$results$diff_raw_df))
  expect_s3_class(plot_volcano(s), "ggplot")
  expect_error(select_comparison(b), "name one")
  expect_error(select_comparison(b, "C_vs_ctrl"), "not in the bundle")
  expect_error(select_comparison(list(), "x"), "run_diff")
})

test_that("summarize_diff_contrasts counts hits per comparison", {
  inp <- mc_input()
  b <- run_diff(inp, method = "limma", group_col = "group",
                control_group = "ctrl", case_group = c("A", "B"))
  s <- summarize_diff_contrasts(b, p_cutoff = 0.05, effect_cutoff = 1)
  expect_identical(s$comparison, c("A_vs_ctrl", "B_vs_ctrl"))
  expect_gte(s$n_up[1], 5L)
  expect_gte(s$n_down[2], 5L)
  expect_equal(s$n_hits, s$n_up + s$n_down)
  expect_s3_class(plot_diff_contrasts(b), "ggplot")
})

test_that("a paired design with several treatments per donor is accepted", {
  inp <- mc_input()
  b <- run_diff(inp, method = "limma", group_col = "group",
                control_group = "ctrl", case_group = c("A", "B"),
                paired_col = "donor")
  expect_identical(diff_comparisons(b), c("A_vs_ctrl", "B_vs_ctrl"))
  bad <- inp
  bad$meta_df$donor[bad$meta_df$group == "ctrl"] <- paste0("X", 1:4)
  expect_error(run_diff(bad, method = "limma", group_col = "group",
                        control_group = "ctrl", case_group = c("A", "B"),
                        paired_col = "donor"),
               "paired design")
})

test_that("case_group and control_group must not overlap", {
  inp <- mc_input()
  expect_error(run_diff(inp, method = "limma", group_col = "group",
                        control_group = "ctrl", case_group = c("A", "ctrl")),
               "distinct")
  expect_error(run_diff(inp, method = "limma", group_col = "group",
                        control_group = "ctrl", case_group = c("A", "A")),
               "case_group")
})

test_that("the ANOVA F-test asks whether groups differ, not whether means are zero", {
  set.seed(3)
  groups <- rep(c("a", "b", "c"), each = 5)
  m <- matrix(rnorm(200 * 15, mean = 20, sd = 1), 200,
              dimnames = list(paste0("F", 1:200), paste0("S", 1:15)))
  m[1:10, groups == "c"] <- m[1:10, groups == "c"] + 4
  inp <- omics_input(m, data.frame(group = groups, row.names = colnames(m)),
                     data.frame(feature_id = rownames(m), row.names = rownames(m)),
                     omics_type = "proteomics", assay_type = "normalized_intensity")
  b <- run_diff(inp, method = "limma", analysis_type = "anova", group_col = "group")
  df <- b$results$diff_result_df
  # Pure noise elsewhere: most features must not be significant. The
  # cell-means F-test made all 200 significant (means of ~20 are not 0).
  expect_lt(sum(df$adj_p_value[11:200] < 0.05), 20)
  expect_true(all(df$adj_p_value[1:10] < 0.05))
  expect_true(all(is.na(df$is_significant)))
})

test_that("log-scale RNA-seq is analysed with limma under auto, so covariates are kept", {
  skip_if_not_installed("limma")
  inp <- mc_input()
  inp$omics_type <- "rnaseq"
  inp$assay_type <- "logcpm"
  expect_identical(auto_select_diff_method(inp, "group"), "limma")
  b <- run_diff(inp, group_col = "group", control_group = "ctrl",
                case_group = "A", covariates = NULL)
  expect_identical(b$params$method, "limma")
})

test_that("an argument the chosen engine cannot use is reported, not silently dropped", {
  inp <- mc_input()
  inp$meta_df$age <- seq_len(nrow(inp$meta_df))
  expect_warning(
    b <- run_diff(inp, method = "ttest", group_col = "group",
                  control_group = "ctrl", case_group = "A", covariates = "age"),
    "does not support `covariates`")
  expect_match(b$warnings, "ignored", all = FALSE)
})

test_that("the exported script reproduces a shared fit, then takes the same contrast", {
  inp <- mc_input()
  b <- run_diff(inp, method = "limma", group_col = "group",
                control_group = "ctrl", case_group = c("A", "B"))
  s <- select_comparison(b, "B_vs_ctrl")
  proj <- omics_project("p", list(proteomics = inp))
  proj$bundles <- list(diff = s)
  code <- export_script(proj, include_plots = FALSE)
  txt <- paste(code, collapse = "\n")
  expect_match(txt, 'case_group    = c("A", "B")', fixed = TRUE)
  expect_match(txt, 'diff <- select_comparison(diff, "B_vs_ctrl")', fixed = TRUE)
})
