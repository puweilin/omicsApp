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

# ---- contrasts beyond "each vs control" ----------------------------------

test_that("pairwise_contrasts lists every pair, later against earlier", {
  expect_identical(pairwise_contrasts(c("ctrl", "A", "B")),
                   c("A - ctrl", "B - ctrl", "B - A"))
  expect_identical(pairwise_contrasts(c("Control", "Drug A")),
                   "`Drug A` - Control")
  expect_error(pairwise_contrasts("A"), "two groups")
})

test_that("contrast expressions become zero-sum weights, and nothing else is accepted", {
  expect_identical(contrast_labels(c("B - A", "(B + C)/2 - A"), c("A", "B", "C")),
                   c("B_vs_A", "(B + C)/2 - A"))
  w <- contrast_weights("(B + C)/2 - A", c("A", "B", "C"))
  expect_equal(unname(w), c(-1, 0.5, 0.5))
  expect_identical(contrast_labels("`Drug A` - `24h`", c("24h", "Drug A")), "Drug A_vs_24h")
  expect_error(contrast_labels("B + A", c("A", "B")), "sum to")
  expect_error(contrast_labels("D - A", c("A", "B")), "not a group")
  expect_error(contrast_labels("B * A", c("A", "B")), "multiplies")
  expect_error(contrast_labels("B - A + 1", c("A", "B")), "adds a number")
  expect_error(contrast_labels("B / 0 - A", c("A", "B")), "divides")
  expect_error(contrast_labels("unlink('x') - A", c("A", "B")), "something other")
  expect_error(contrast_labels(c("B - A", "B - A"), c("A", "B")), "twice")
})

test_that("run_diff(contrasts = 'pairwise') fits every pair in one model, control first", {
  inp <- mc_input()
  b <- run_diff(inp, method = "limma", group_col = "group",
                control_group = "ctrl", contrasts = "pairwise")
  expect_identical(diff_comparisons(b), c("A_vs_ctrl", "B_vs_ctrl", "B_vs_A"))
  # The two "vs control" contrasts are exactly the each-vs-control run.
  each <- run_diff(inp, method = "limma", group_col = "group",
                   control_group = "ctrl", case_group = c("A", "B"))
  expect_equal(select_comparison(b, "B_vs_ctrl")$results$diff_result_df$p_value,
               select_comparison(each, "B_vs_ctrl")$results$diff_result_df$p_value)
  s <- select_comparison(b, "B_vs_A")
  expect_identical(s$params$case_group, "B")
  expect_identical(s$params$control_group, "A")
  expect_identical(s$params$contrasts, "B - A")
})

test_that("a weighted contrast is read off the same fit in limma, edgeR and DESeq2", {
  inp <- mc_input()
  b <- run_diff(inp, method = "limma", group_col = "group",
                contrasts = c("(A + B)/2 - ctrl"))
  df <- b$results$diff_result_df
  expect_identical(unique(df$comparison), "(A + B)/2 - ctrl")
  # Features 1-5 are +3 in A only, so half of that on the average.
  expect_equal(mean(df$effect[1:5]), 1.5, tolerance = 0.3)
  expect_null(b$params$case_group)

  skip_if_not_installed("edgeR")
  skip_if_not_installed("DESeq2")
  cnt <- mc_input(omics = "rnaseq")
  for (m in c("edger", "deseq2")) {
    r <- suppressWarnings(suppressMessages(
      run_diff(cnt, method = m, group_col = "group",
               contrasts = c("(A + B)/2 - ctrl", "B - A"))))
    expect_identical(diff_comparisons(r), c("(A + B)/2 - ctrl", "B_vs_A"))
    avg <- select_comparison(r, "(A + B)/2 - ctrl")$results$diff_result_df
    # log2(6) / 2 ~ 1.3 on the A-only genes.
    expect_true(all(avg$effect[1:5] > 0.8 & avg$effect[1:5] < 1.8), info = m)
    ba <- select_comparison(r, "B_vs_A")$results$diff_result_df
    expect_true(all(ba$effect[1:5] < -1.5), info = m)
  }
})

test_that("the per-pair engines take pairs but refuse weighted contrasts", {
  inp <- mc_input()
  b <- run_diff(inp, method = "ttest", group_col = "group", contrasts = "B - A")
  expect_identical(diff_comparisons(b), "B_vs_A")
  expect_error(run_diff(inp, method = "ttest", group_col = "group",
                        contrasts = "(A + B)/2 - ctrl"),
               "needs limma")
})

test_that("the exported script repeats a contrast run and takes the same comparison", {
  inp <- mc_input()
  b <- run_diff(inp, method = "limma", group_col = "group",
                control_group = "ctrl", contrasts = "pairwise")
  proj <- omics_project("p", list(proteomics = inp))
  proj$bundles <- list(diff = select_comparison(b, "B_vs_A"))
  txt <- paste(export_script(proj, include_plots = FALSE), collapse = "\n")
  expect_match(txt, 'contrasts     = c("A - ctrl", "B - ctrl", "B - A")', fixed = TRUE)
  expect_match(txt, 'diff <- select_comparison(diff, "B_vs_A")', fixed = TRUE)
})

# ---- overlap -----------------------------------------------------------

test_that("hit sets and the overlap plot follow the thresholds", {
  inp <- mc_input(levels = c("ctrl", "A", "B"))
  # Features 1-5 up in A, 6-10 down in B; add 11-15 up in both.
  inp$expr_mat[11:15, inp$meta_df$group != "ctrl"] <-
    inp$expr_mat[11:15, inp$meta_df$group != "ctrl"] + 3
  b <- run_diff(inp, method = "limma", group_col = "group",
                control_group = "ctrl", case_group = c("A", "B"))
  sets <- diff_hit_sets(b, effect_cutoff = 1)
  expect_identical(names(sets), c("A_vs_ctrl", "B_vs_ctrl"))
  expect_true(all(paste0("F", 11:15) %in% intersect(sets[[1]], sets[[2]])))
  expect_true(all(paste0("F", 1:5) %in% setdiff(sets[[1]], sets[[2]])))
  up <- diff_hit_sets(b, effect_cutoff = 1, direction = "up")
  expect_false(any(paste0("F", 6:10) %in% up$B_vs_ctrl))
  expect_s3_class(plot_diff_overlap(b, effect_cutoff = 1), "ggplot")
  # One comparison, or no hits: a message, not an error.
  one <- run_diff(inp, method = "limma", group_col = "group",
                  control_group = "ctrl", case_group = "A")
  expect_s3_class(plot_diff_overlap(one), "ggplot")
  expect_s3_class(plot_diff_overlap(b, p_cutoff = 1e-300), "ggplot")
})

test_that("a multi-contrast run is exported one table per comparison", {
  inp <- mc_input()
  b <- run_diff(inp, method = "limma", group_col = "group",
                contrasts = c("A - ctrl", "(A + B)/2 - ctrl"))
  dir <- withr::local_tempdir()
  reg <- export_bundle(b, dir, formats = "tsv")
  files <- basename(reg$path)
  expect_true("run_diff_A_vs_ctrl_diff_result_df.tsv" %in% files)
  expect_true("run_diff_A_B_2_ctrl_diff_result_df.tsv" %in% files)
  expect_true("run_diff_contrast_summary.tsv" %in% files)
  one <- utils::read.delim(file.path(dir, "run_diff_A_vs_ctrl_diff_result_df.tsv"))
  expect_identical(unique(one$comparison), "A_vs_ctrl")
})

test_that("the report covers every comparison, with its own top hits", {
  skip_if_not_installed("rmarkdown")
  skip_if_not(rmarkdown::pandoc_available(), "pandoc is not available")
  inp <- mc_input()
  b <- run_diff(inp, method = "limma", group_col = "group",
                control_group = "ctrl", case_group = c("A", "B"))
  proj <- omics_project("Multi", list(proteomics = inp))
  proj$bundles <- list(diff = b)
  out <- withr::local_tempfile(fileext = ".html")
  export_report(proj, out)
  html <- paste(readLines(out, warn = FALSE), collapse = "\n")
  expect_match(html, "A vs ctrl", fixed = TRUE)
  expect_match(html, "B vs ctrl", fixed = TRUE)
  expect_match(html, "fitted together in one model", fixed = TRUE)
})

test_that("the script takes the comparison enrichment ran on out of a full run", {
  inp <- mc_input()
  b <- run_diff(inp, method = "limma", group_col = "group",
                control_group = "ctrl", case_group = c("A", "B"))
  proj <- omics_project("p", list(proteomics = inp))
  enrich <- new_analysis_bundle("run_enrichment",
                                params = list(type = "ora", database = "hallmark",
                                              comparison = "B_vs_ctrl"))
  proj$bundles <- list(diff = b, enrich = enrich)
  txt <- paste(export_script(proj, include_plots = FALSE), collapse = "\n")
  expect_match(txt, 'diff_shown <- select_comparison(diff, "B_vs_ctrl")', fixed = TRUE)
  expect_match(txt, "run_enrichment(\n  diff_shown", fixed = TRUE)
})

# ---- global test across groups ------------------------------------------

test_that("counts get a global test from edgeR and DESeq2 too", {
  skip_if_not_installed("edgeR")
  skip_if_not_installed("DESeq2")
  cnt <- mc_input(omics = "rnaseq")
  for (m in c("edger", "deseq2")) {
    b <- suppressWarnings(suppressMessages(
      run_diff(cnt, method = m, analysis_type = "anova", group_col = "group")))
    df <- b$results$diff_result_df
    check_diff_result_schema(df)
    expect_identical(unique(df$analysis_type), "anova")
    # The seeded genes (1-10) differ between groups; most others do not.
    expect_true(all(df$adj_p_value[1:10] < 0.05), info = m)
    expect_lt(mean(df$adj_p_value[11:40] < 0.05), 0.2)
    expect_true(all(is.na(df$is_significant)))
  }
  expect_error(run_diff(cnt, method = "ttest", analysis_type = "anova",
                        group_col = "group"), "ANOVA")
})

# ---- paired t-test with several treatments ------------------------------

# Six patients, all sampled under Control and TreatA; only P1-P4 under
# TreatB. Samples in a scrambled order, so pairs are found by patient and
# not by position.
pt_input <- function() {
  set.seed(31)
  pts <- paste0("P", 1:6)
  meta <- data.frame(patient = c(pts, pts, pts[1:4]),
                     group = rep(c("Control", "TreatA", "TreatB"), c(6, 6, 4)),
                     stringsAsFactors = FALSE)
  rownames(meta) <- paste0(meta$patient, "_", meta$group)
  meta <- meta[sample(nrow(meta)), ]
  n_feat <- 25L
  patient_effect <- stats::rnorm(6, 0, 1)
  m <- matrix(stats::rnorm(n_feat * nrow(meta), 10, 0.3), n_feat,
              dimnames = list(paste0("F", seq_len(n_feat)), rownames(meta)))
  m <- m + rep(patient_effect[match(meta$patient, pts)], each = n_feat)
  m[1:5, meta$group == "TreatA"] <- m[1:5, meta$group == "TreatA"] + 1
  m[6:10, meta$group == "TreatB"] <- m[6:10, meta$group == "TreatB"] - 1
  omics_input(m, meta, data.frame(feature_id = rownames(m)),
              omics_type = "proteomics", assay_type = "normalized_intensity")
}

pt_manual <- function(inp, case, control, patients) {
  m <- inp$expr_mat
  t(vapply(rownames(m), function(f) {
    x <- m[f, paste0(patients, "_", case)]
    y <- m[f, paste0(patients, "_", control)]
    tt <- stats::t.test(x, y, paired = TRUE)
    c(t = unname(tt$statistic), p = tt$p.value, d = mean(x - y))
  }, numeric(3)))
}

test_that("a paired t-test with several treatments pairs each comparison on its own", {
  inp <- pt_input()
  b <- run_diff(inp, method = "ttest", group_col = "group", control_group = "Control",
                case_group = c("TreatA", "TreatB"), paired_col = "patient")
  expect_identical(diff_comparisons(b), c("TreatA_vs_Control", "TreatB_vs_Control"))
  df <- b$results$diff_result_df
  for (cmp in list(list("TreatA", paste0("P", 1:6)), list("TreatB", paste0("P", 1:4)))) {
    got <- df[df$comparison == paste0(cmp[[1]], "_vs_Control"), ]
    want <- pt_manual(inp, cmp[[1]], "Control", cmp[[2]])
    got <- got[match(rownames(want), got$feature_id), ]
    expect_equal(got$statistic, unname(want[, "t"]), info = cmp[[1]])
    expect_equal(got$p_value, unname(want[, "p"]), info = cmp[[1]])
    expect_equal(got$effect, unname(want[, "d"]), info = cmp[[1]])
  }
  # Said in words, with the count for each comparison and who was left out.
  note <- grep("Paired t-test", b$warnings, value = TRUE)
  expect_length(note, 1L)
  expect_match(note, "TreatA vs Control used 6 pairs", fixed = TRUE)
  expect_match(note, "TreatB vs Control used 4 pairs (2 left out", fixed = TRUE)
  expect_match(note, "P5, P6", fixed = TRUE)
})

test_that("pairwise paired comparisons use the patients both groups share", {
  inp <- pt_input()
  b <- run_diff(inp, method = "ttest", group_col = "group", control_group = "Control",
                contrasts = "pairwise", paired_col = "patient")
  df <- b$results$diff_result_df
  got <- df[df$comparison == "TreatB_vs_TreatA", ]
  want <- pt_manual(inp, "TreatB", "TreatA", paste0("P", 1:4))
  got <- got[match(rownames(want), got$feature_id), ]
  expect_equal(got$statistic, unname(want[, "t"]))
  expect_match(b$warnings, "TreatB vs TreatA used 4 pairs", fixed = TRUE, all = FALSE)
})

test_that("a comparison with fewer than two complete pairs is refused", {
  inp <- pt_input()
  keep <- !(inp$meta_df$group == "TreatB" & inp$meta_df$patient %in% paste0("P", 2:4))
  inp <- subset_omics(inp, samples = rownames(inp$meta_df)[keep])
  expect_error(run_diff(inp, method = "ttest", group_col = "group",
                        control_group = "Control", case_group = c("TreatA", "TreatB"),
                        paired_col = "patient"),
               "TreatB against Control needs at least 2 complete pairs")
})

test_that("a subject sampled twice under one group is left out of that comparison", {
  inp <- pt_input()
  inp$meta_df$patient[rownames(inp$meta_df) == "P6_TreatA"] <- "P1"
  b <- run_diff(inp, method = "ttest", group_col = "group", control_group = "Control",
                case_group = c("TreatA", "TreatB"), paired_col = "patient")
  expect_match(b$warnings,
               "TreatA vs Control used 4 pairs (2 left out, without one sample in each group: P1, P6)",
               fixed = TRUE, all = FALSE)
})

test_that("a balanced paired design needs no note, and one comparison pairs the same way", {
  inp <- pt_input()
  bal <- subset_omics(inp, samples = rownames(inp$meta_df)[inp$meta_df$patient %in% paste0("P", 1:4)])
  b <- run_diff(bal, method = "ttest", group_col = "group", control_group = "Control",
                case_group = c("TreatA", "TreatB"), paired_col = "patient")
  expect_false(any(grepl("Paired t-test", b$warnings)))
  # A single comparison with a patient missing its partner uses the
  # complete pairs, as the same comparison does inside a multi-group run.
  one <- run_diff(inp, method = "ttest", group_col = "group",
                  control_group = "Control", case_group = "TreatB",
                  paired_col = "patient")
  multi <- run_diff(inp, method = "ttest", group_col = "group",
                    control_group = "Control", case_group = c("TreatA", "TreatB"),
                    paired_col = "patient")
  m <- multi$results$diff_result_df
  m <- m[m$comparison == "TreatB_vs_Control", , drop = FALSE]
  expect_equal(one$results$diff_result_df$p_value, m$p_value)
  expect_match(one$warnings, "TreatB vs Control used", fixed = TRUE, all = FALSE)
})
