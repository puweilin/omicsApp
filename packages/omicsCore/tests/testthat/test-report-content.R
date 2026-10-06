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
  # In words, not internal labels.
  expect_match(html, "log2 intensities", fixed = TRUE)
  expect_match(html, "raw read counts", fixed = TRUE)
  expect_match(html, as.character(nrow(inp$expr_mat)), fixed = TRUE)
  expect_match(html, as.character(ncol(inp$expr_mat)), fixed = TRUE)
})

test_that("the report describes the analysis and lists its strongest hits", {
  skip_if_no_report()
  inp <- realistic_input(n_per_group = 3L)
  diff <- run_diff(inp, method = "ttest", analysis_type = "group",
                   group_col = "group", control_group = "G1", case_group = "G2")
  p <- omics_project("With results", experiments = list(proteomics = inp))
  p$bundles <- list(diff = diff)
  html <- report_html(p)
  expect_match(html, "Methods", fixed = TRUE)
  expect_match(html, "t-test", fixed = TRUE)
  expect_match(html, "G2 vs G1", fixed = TRUE)
  # The top hits are the most significant features, not the first rows
  # of the table in file order.
  df <- diff$results$diff_result_df
  best <- df$feature_symbol[which.min(df$p_value)]
  expect_match(html, best, fixed = TRUE)
  expect_match(html, "Top features by p-value", fixed = TRUE)
  # And p-values are written out, never rounded to 0.
  expect_false(grepl("<td[^>]*>0</td>", html))
  expect_false(grepl("No analyses have been run yet", html, fixed = TRUE))
})

test_that("a project with no analyses says so rather than rendering nothing", {
  skip_if_no_report()
  p <- omics_project("Empty", experiments = list(proteomics = realistic_input(n_per_group = 3L)))
  html <- report_html(p)
  expect_match(html, "No analyses have been run yet", fixed = TRUE)
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
  strip <- function(html) gsub("<img[^>]*>", "", sub("Software versions.*$", "", html))
  strip_dates <- function(html) gsub("[0-9]{4}-[0-9]{2}-[0-9]{2}[^<]*", "", html)
  expect_identical(strip_dates(strip(live)), strip_dates(strip(restored)))
})

test_that("the report says which gene lists ORA tested and which way pathways went", {
  skip_if_no_report()
  inp <- realistic_input(n_per_group = 3L)
  diff <- run_diff(inp, method = "limma", analysis_type = "group",
                   group_col = "group", control_group = "G1", case_group = "G2")
  enr <- data.frame(
    database = "hallmark", result_type = "ora", comparison = "G2_vs_G1",
    pathway_id = c("HALLMARK_G2M_CHECKPOINT", "HALLMARK_ADIPOGENESIS"),
    pathway_name = c("G2M CHECKPOINT", "ADIPOGENESIS"),
    effect = NA_real_, effect_type = NA_character_, direction = c("up", "down"),
    p_value = c(1e-8, 1e-3), adj_p_value = c(1e-7, 1e-2), q_value = NA_real_,
    gene_set_size = 40, overlap_size = c(30, 8), overlap_features = "A",
    leading_features = "A", source_label = "ora", stringsAsFactors = FALSE)
  eb <- new_analysis_bundle("run_enrichment", input_info = diff$input_info,
                            params = list(type = "ora", database = c("hallmark", "kegg"),
                                          organism = "Homo sapiens",
                                          direction = "separate",
                                          p_adjust_scope = "all"),
                            results = list(enrich_result_df = enr))
  p <- omics_project("Directions", experiments = list(proteomics = inp))
  p$bundles <- list(diff = diff, enrich = eb)
  # Pandoc wraps long lines; a phrase can straddle one.
  html <- gsub("\\s+", " ", report_html(p))
  expect_match(html, "tested as two separate lists", fixed = TRUE)
  expect_match(html, "adjusted across all databases together", fixed = TRUE)
  expect_match(html, "up-regulated genes", fixed = TRUE)
  expect_match(html, "down-regulated genes", fixed = TRUE)

  # An older, pooled result is described as pooled.
  eb$params$direction <- "both"
  eb$params$database <- "hallmark"
  eb$results$enrich_result_df$direction <- NA_character_
  p$bundles$enrich <- eb
  html <- gsub("\\s+", " ", report_html(p))
  expect_match(html, "up and down pooled into one list", fixed = TRUE)

  # A directional ActivePathways result: the way each pathway went.
  idf <- data.frame(
    feature_id = c("P1", "P2"), feature_symbol = c("HALLMARK_ONE", "HALLMARK_TWO"),
    result_type = "active_pathways", experiments = "a vs b", comparison = "x | y",
    effect = c(10, 6), effect_type = "neg_log10_padj", statistic = c(1e-10, 1e-6),
    statistic_type = "adjusted_p_val", p_value = c(1e-10, 1e-6),
    adj_p_value = c(1e-10, 1e-6), direction = c("up", "mixed"), quadrant = "a,b",
    is_significant = TRUE, source_label = "ap", evidence = c("shared", "combined"),
    stringsAsFactors = FALSE)
  p$bundles$integration <- new_analysis_bundle(
    "run_integration",
    params = list(method = "active_pathways", experiments = c("a", "b"),
                  merge_method = "DPM", constraints_vector = c(1, 1)),
    results = list(integration_df = idf))
  html <- gsub("\\s+", " ", report_html(p))
  expect_match(html, "2 pathways were significant: 1 up in both layers, 0 down in both, 1 mixed",
               fixed = TRUE)
  expect_match(html, "expecting the layers to change in the same direction", fixed = TRUE)
  expect_match(html, "Found by", fixed = TRUE)
})
