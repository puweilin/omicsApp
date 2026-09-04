# What the HTML report says, not whether the file exists.
#
# test-export.R checks that export_report() writes a file of non-zero
# size. A report that renders but leaves out the analyses would pass
# that. These parse the HTML and look for the numbers and names a
# reader would check against the app.

skip_if_no_report <- function() {
  skip_if_not_installed("rmarkdown")
  skip_if_not(rmarkdown::pandoc_available(), "pandoc unavailable")
}

report_html <- function(project) {
  path <- tempfile(fileext = ".html")
  suppressMessages(export_report(project, path, format = "html"))
  paste(readLines(path, warn = FALSE), collapse = "\n")
}

test_that("the report names every layer with its shape and label", {
  skip_if_no_report()
  inp <- realistic_input(n_per_group = 3L)
  rna <- realistic_input("rnaseq", n_per_group = 3L)
  p <- omics_project("Report project",
                     experiments = list(proteomics = inp, rnaseq = rna))
  html <- report_html(p)
  expect_match(html, "Report project", fixed = TRUE)
  for (tag in c("proteomics", "rnaseq")) expect_match(html, tag, fixed = TRUE)
  expect_match(html, "normalized_intensity", fixed = TRUE)
  expect_match(html, "raw_count", fixed = TRUE)
  expect_match(html, as.character(nrow(inp$expr_mat)), fixed = TRUE)
  expect_match(html, as.character(ncol(inp$expr_mat)), fixed = TRUE)
})

test_that("the report carries each analysis and the head of its result table", {
  skip_if_no_report()
  inp <- realistic_input(n_per_group = 3L)
  diff <- run_diff(inp, method = "ttest", analysis_type = "group",
                   group_col = "group", control_group = "G1", case_group = "G2")
  p <- omics_project("With results", experiments = list(proteomics = inp))
  p$bundles <- list(diff = diff)
  html <- report_html(p)
  expect_match(html, "run_diff", fixed = TRUE)
  expect_match(html, "diff_result_df", fixed = TRUE)
  expect_match(html, diff$params$comparison, fixed = TRUE)
  # The first rows of the result table are printed; the first feature
  # and its symbol must be there, verbatim.
  head_rows <- utils::head(diff$results$diff_result_df, 10L)
  expect_match(html, head_rows$feature_id[[1L]], fixed = TRUE)
  expect_match(html, head_rows$feature_symbol[[1L]], fixed = TRUE)
  expect_false(grepl("No analysis_bundle objects", html, fixed = TRUE))
})

test_that("a project with no analyses says so rather than rendering nothing", {
  skip_if_no_report()
  p <- omics_project("Empty", experiments = list(proteomics = realistic_input(n_per_group = 3L)))
  html <- report_html(p)
  expect_match(html, "No analysis_bundle objects", fixed = TRUE)
})

test_that("a restored project reports the same as the live one", {
  skip_if_no_report()
  skip_if_not_installed("qs2")
  inp <- realistic_input(n_per_group = 3L)
  diff <- run_diff(inp, method = "ttest", analysis_type = "group",
                   group_col = "group", control_group = "G1", case_group = "G2")
  p <- omics_project("Round trip", experiments = list(proteomics = inp))
  p$bundles <- list(diff = diff)
  f <- tempfile(fileext = ".omp")
  save_project(p, f)
  live <- report_html(p)
  restored <- report_html(load_project(f))
  # Everything but the timestamp and session info
  strip <- function(html) sub("Session info.*$", "", html)
  strip_dates <- function(html) gsub("[0-9]{4}-[0-9]{2}-[0-9]{2}[^<]*", "", html)
  expect_identical(strip_dates(strip(live)), strip_dates(strip(restored)))
})
