# The group-aware missing-value filter: a feature is kept when it is
# observed well enough within a group, rather than over the whole study.

# Six samples, A A A B B B, plus one sample with no group. Rows built so
# each rule gives a different answer at a cutoff of 0.4.
mf_input <- function(design = FALSE) {
  o <- 20
  x <- rbind(
    on_off   = c(o,  o,  o,  NA, NA, NA, o),   # all A, no B
    complete = c(o,  o,  o,  o,  o,  o,  o),
    one_each = c(NA, o,  o,  NA, o,  o,  o),   # 1/3 missing in each group
    two_in_a = c(NA, NA, o,  o,  o,  o,  o),   # 2/3 of A, none of B
    sparse   = c(NA, NA, o,  NA, NA, o,  NA)   # 2/3 in both groups
  )
  x <- x + matrix(seq_len(length(x)) / 100, nrow(x))
  colnames(x) <- paste0("s", 1:7)
  meta <- data.frame(group = c("A", "A", "A", "B", "B", "B", NA),
                     batch = rep(c("x", "y"), length.out = 7),
                     row.names = colnames(x))
  inp <- omics_input(x, meta, data.frame(feature_id = rownames(x), row.names = rownames(x)),
                     omics_type = "proteomics", assay_type = "normalized_intensity")
  if (design) inp <- set_study_design(inp, "group", "A")
  inp
}

test_that("each rule flags the right features on a hand-built matrix", {
  inp <- mf_input()
  global <- qc_missingness(inp, feature_missing_cutoff = 0.4)
  # The ungrouped sample counts here: on_off is 3/7 missing.
  expect_setequal(global$flagged_features, c("on_off", "sparse"))

  any_g <- qc_missingness(inp, feature_missing_cutoff = 0.4,
                          missing_filter = "any_group", group_col = "group")
  expect_setequal(any_g$flagged_features, "sparse")
  expect_equal(any_g$feature_metrics$filter_missing_rate,
               c(0, 0, 1 / 3, 0, 2 / 3))

  all_g <- qc_missingness(inp, feature_missing_cutoff = 0.4,
                          missing_filter = "all_groups", group_col = "group")
  expect_setequal(all_g$flagged_features, c("on_off", "two_in_a", "sparse"))
  expect_equal(all_g$feature_metrics$filter_missing_rate,
               c(1, 0, 1 / 3, 2 / 3, 2 / 3))
  expect_identical(all_g$settings$missing_filter, "all_groups")
  expect_identical(all_g$settings$group_col, "group")

  # The cutoff keeps its meaning: at or below is kept.
  expect_setequal(qc_missingness(inp, feature_missing_cutoff = 0.7,
                                 missing_filter = "all_groups",
                                 group_col = "group")$flagged_features,
                  "on_off")
})

test_that("run_qc filters by group and records the resolved group column", {
  inp <- mf_input(design = TRUE)
  b <- run_qc(inp, missing_threshold = 0.4, missing_filter = "any_group",
              outlier_method = "none", impute_method = "none")
  expect_identical(b$params$missing_filter, "any_group")
  # From the recorded study design.
  expect_identical(b$params$group_col, "group")
  expect_setequal(b$results$qc_summary$recommended_filters$remove_features, "sparse")
  expect_false("sparse" %in% rownames(qc_cleaned_input(b, inp)$expr_mat))
  expect_true("on_off" %in% rownames(qc_cleaned_input(b, inp)$expr_mat))

  b2 <- run_qc(inp, missing_threshold = 0.4, missing_filter = "all_groups",
               group_col = "batch", outlier_method = "none", impute_method = "none")
  expect_identical(b2$params$group_col, "batch")
})

test_that("the default stays global and gives the old result", {
  inp <- mf_input(design = TRUE)
  b <- run_qc(inp, missing_threshold = 0.4, outlier_method = "none",
              impute_method = "none")
  expect_identical(b$params$missing_filter, "global")
  expect_null(b$params$group_col)
  expect_setequal(b$results$qc_summary$recommended_filters$remove_features,
                  c("on_off", "sparse"))
  expect_null(b$results$qc_summary$missingness$feature_metrics$filter_missing_rate)
  # A group column given to the global rule changes nothing.
  b2 <- run_qc(inp, missing_threshold = 0.4, outlier_method = "none",
               impute_method = "none", group_col = "batch")
  expect_identical(qc_cleaned_input(b2, inp)$expr_mat, qc_cleaned_input(b, inp)$expr_mat)
  expect_null(b2$params$group_col)
})

test_that("a group rule without a group column says what it needs", {
  inp <- mf_input()
  expect_error(run_qc(inp, missing_filter = "any_group", outlier_method = "none"),
               "needs `group_col`", fixed = TRUE)
  expect_error(qc_missingness(inp, missing_filter = "any_group", group_col = "nope"),
               "not a column", fixed = TRUE)
  expect_error(run_qc(inp, missing_filter = "per_batch"), "should be one of")
})

test_that("the exported script repeats the group filter and reproduces it", {
  inp <- mf_input(design = TRUE)
  inp$source_path <- "raw/mf.xlsx"
  original <- run_qc(inp, missing_threshold = 0.4, missing_filter = "all_groups",
                     outlier_method = "none", impute_method = "none")
  proj <- omics_project("mf", experiments = list(proteomics = inp))
  proj$bundles <- list(qc = original)
  lines <- export_script(proj, include_plots = FALSE)
  call_start <- grep("^qc <- run_qc\\($", lines)
  expect_length(call_start, 1L)
  call_end <- call_start + which(lines[-seq_len(call_start)] == ")")[1L]
  call <- lines[call_start:call_end]
  expect_true(any(grepl('missing_filter *= "all_groups"', call)))
  expect_true(any(grepl('group_col *= "group"', call)))

  # Run the emitted call against the same layer: same features out.
  env <- new.env(parent = asNamespace("omicsCore"))
  env$input <- inp
  eval(parse(text = call), envir = env)
  expect_identical(rownames(qc_cleaned_input(env$qc, inp)$expr_mat),
                   rownames(qc_cleaned_input(original, inp)$expr_mat))
  expect_identical(env$qc$params, original$params)

  # A bundle from before the option existed has no missing_filter, and
  # the script leaves it out: the default is the old behaviour.
  old <- original
  old$params$missing_filter <- NULL
  old$params$group_col <- NULL
  proj$bundles <- list(qc = old)
  expect_false(any(grepl("missing_filter", export_script(proj, include_plots = FALSE))))
})
