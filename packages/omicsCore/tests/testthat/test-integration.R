# Tests for the integration slice (correlation, concordance, active_pathways,
# dispatcher, plot_integration). Heavy backends are gated by
# skip_if_not_installed() so the suite stays green on a fresh install.

# Two layers (proteomics + rnaseq) sharing donors via a sample_link, plus an
# overlapping subset of HGNC symbols so cross-omics joins are non-empty.
make_integration_project <- function(n_features = 60L, n_donors = 8L) {
  set.seed(2026)
  symbols <- c(
    "TP53", "EGFR", "MYC", "AKT1", "BRCA1", "BRCA2", "PTEN", "RB1", "KRAS",
    "PIK3CA", "STAT3", "JUN", "FOS", "VEGFA", "TNF", "IL6", "CDKN1A",
    "MAPK1", "MAPK3", "BAX", "BCL2", "CASP3", "CASP8", "ATM", "ATR",
    "CHEK1", "CHEK2", "MDM2", "RAD51", "XRCC1"
  )
  symbols <- c(symbols, paste0("GENE", seq_len(max(0L, n_features - length(symbols)))))
  symbols <- symbols[seq_len(n_features)]

  groups <- rep(c("ctrl", "case"), each = n_donors / 2L)

  # Proteomics layer
  feat_a <- paste0("p", seq_len(n_features))
  samp_a <- paste0("PA_", seq_len(n_donors))
  mat_a <- matrix(rnorm(n_features * n_donors, mean = 5, sd = 1),
                  nrow = n_features,
                  dimnames = list(feat_a, samp_a))
  mat_a[1:10, groups == "case"] <- mat_a[1:10, groups == "case"] + 2.5
  meta_a <- data.frame(group = groups, row.names = samp_a,
                       stringsAsFactors = FALSE)
  fdf_a <- data.frame(
    feature_id = feat_a,
    feature_symbol = symbols,
    row.names = feat_a,
    stringsAsFactors = FALSE
  )
  input_a <- omics_input(mat_a, meta_a, fdf_a, omics_type = "proteomics",
                          assay_type = "normalized_intensity")

  # RNA-seq layer using same donors but distinct sample IDs (so we exercise
  # sample_link). Symbols overlap perfectly so feature_pairs is non-empty;
  # the first 10 are pumped up in the case group to be concordant with
  # proteomics, and the next 5 are pumped DOWN in case to create a
  # discordant slice for quadrant testing.
  feat_b <- paste0("r", seq_len(n_features))
  samp_b <- paste0("RB_", seq_len(n_donors))
  base <- matrix(rnorm(n_features * n_donors, mean = 6, sd = 0.5),
                 nrow = n_features,
                 dimnames = list(feat_b, samp_b))
  base[1:10, groups == "case"] <- base[1:10, groups == "case"] + 2.0
  base[11:15, groups == "case"] <- base[11:15, groups == "case"] - 2.0
  # Force the RNA layer's first 10 features to track the proteomics features
  # across donors so per-feature correlation is high.
  base[1:10, ] <- mat_a[1:10, ] + matrix(
    rnorm(10 * n_donors, sd = 0.2), nrow = 10
  )
  meta_b <- data.frame(group = groups, row.names = samp_b,
                       stringsAsFactors = FALSE)
  fdf_b <- data.frame(
    feature_id = feat_b,
    feature_symbol = symbols,
    row.names = feat_b,
    stringsAsFactors = FALSE
  )
  input_b <- omics_input(base, meta_b, fdf_b, omics_type = "rnaseq",
                          assay_type = "normalized_intensity")

  donor_ids <- paste0("D", seq_len(n_donors))
  link <- data.frame(
    tag = rep(c("proteo", "rna"), each = n_donors),
    sample_id = c(samp_a, samp_b),
    donor_id = rep(donor_ids, times = 2L),
    stringsAsFactors = FALSE
  )

  omics_project(
    "integ_demo",
    experiments = list(proteo = input_a, rna = input_b),
    sample_link = link
  )
}

build_integration_diff_bundles <- function(project) {
  list(
    proteo = run_diff(project$experiments$proteo, method = "ttest",
                      analysis_type = "group", group_col = "group",
                      control_group = "ctrl", case_group = "case"),
    rna    = run_diff(project$experiments$rna, method = "ttest",
                      analysis_type = "group", group_col = "group",
                      control_group = "ctrl", case_group = "case")
  )
}

# ---- helpers / schema --------------------------------------------------

test_that("integration result template has the standard columns", {
  out <- new_integration_result_template()
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 0L)
  expect_true(all(INTEGRATION_RESULT_REQUIRED_COLS %in% colnames(out)))
})

test_that("resolve_experiment_pair returns both tags when project has two", {
  p <- make_integration_project()
  expect_equal(resolve_experiment_pair(p, NULL), c("proteo", "rna"))
  expect_equal(resolve_experiment_pair(p, c("rna", "proteo")), c("rna", "proteo"))
  expect_error(resolve_experiment_pair(p, c("rna", "foo")), "not found")
  expect_error(resolve_experiment_pair(p, "rna"), "length-2")
})

test_that("build_sample_pairs uses sample_link donor map", {
  p <- make_integration_project()
  pairs <- build_sample_pairs(p, "proteo", "rna")
  expect_setequal(colnames(pairs), c("donor_id", "proteo", "rna"))
  expect_equal(nrow(pairs), 8L)
  expect_true(all(startsWith(pairs$proteo, "PA_")))
  expect_true(all(startsWith(pairs$rna, "RB_")))
})

test_that("build_feature_pairs joins by feature_symbol", {
  p <- make_integration_project()
  fp <- build_feature_pairs(p, "proteo", "rna")
  expect_setequal(colnames(fp),
                  c("feature_id", "feature_symbol", "feature_a", "feature_b"))
  expect_gt(nrow(fp), 30L)
})

test_that("classify_concordance_quadrant labels all four quadrants", {
  q <- classify_concordance_quadrant(
    c("up", "down", "up", "down", "ns"),
    c("up", "down", "down", "up", "up")
  )
  expect_equal(q, c("up_up", "down_down", "up_down", "down_up", NA))
})

# ---- correlation backend ----------------------------------------------

test_that("run_integration correlation produces a standard bundle", {
  p <- make_integration_project()
  b <- run_integration(p, method = "correlation",
                       experiments = c("proteo", "rna"))
  expect_true(is_analysis_bundle(b))
  expect_identical(b$analysis_name, "run_integration")
  expect_identical(b$params$method, "correlation")
  df <- b$results$integration_df
  expect_s3_class(df, "data.frame")
  check_integration_result_schema(df)
  expect_gt(nrow(df), 0L)
  expect_true(all(df$result_type == "correlation"))
  # Direction is positive/negative or NA, never up/down.
  expect_true(all(is.na(df$direction) | df$direction %in% c("positive", "negative")))
  # The first 10 features were forced to track across layers, so several
  # of them should land with a high positive correlation.
  expect_true(any(df$effect > 0.7, na.rm = TRUE))
})

test_that("correlation requires the configured minimum sample count", {
  p <- make_integration_project(n_donors = 4L)
  expect_error(
    run_integration(p, method = "correlation", min_samples = 20L),
    "paired samples"
  )
})

# ---- concordance backend ----------------------------------------------

test_that("run_integration concordance produces a standard bundle", {
  p <- make_integration_project()
  diffs <- build_integration_diff_bundles(p)
  b <- run_integration(p, method = "concordance",
                       experiments = c("proteo", "rna"),
                       diff_bundles = diffs)
  expect_true(is_analysis_bundle(b))
  expect_identical(b$params$method, "concordance")
  df <- b$results$integration_df
  check_integration_result_schema(df)
  expect_gt(nrow(df), 0L)
  expect_true(all(df$result_type == "concordance"))
  # Quadrant values come from the classifier (or NA for ns features).
  qok <- df$quadrant %in% c("up_up", "down_down", "up_down", "down_up") |
         is.na(df$quadrant)
  expect_true(all(qok))
  # Direction is concordant/discordant or NA.
  expect_true(all(is.na(df$direction) |
                  df$direction %in% c("concordant", "discordant")))
})

test_that("concordance requires diff_bundles", {
  p <- make_integration_project()
  expect_error(
    run_integration(p, method = "concordance"),
    "diff_bundles"
  )
})

test_that("concordance rejects non-run_diff bundles", {
  p <- make_integration_project()
  diffs <- build_integration_diff_bundles(p)
  diffs$rna$analysis_name <- "not_run_diff"
  expect_error(
    run_integration(p, method = "concordance",
                    experiments = c("proteo", "rna"),
                    diff_bundles = diffs),
    "run_diff"
  )
})

# ---- active_pathways backend (Suggests-gated) -------------------------

test_that("run_integration active_pathways requires ActivePathways", {
  p <- make_integration_project()
  diffs <- build_integration_diff_bundles(p)
  if (requireNamespace("ActivePathways", quietly = TRUE)) {
    skip("ActivePathways is installed; skipping the gate test")
  }
  expect_error(
    run_integration(p, method = "active_pathways",
                    experiments = c("proteo", "rna"),
                    diff_bundles = diffs),
    "ActivePathways"
  )
})

test_that("active_pathways returns the integration schema when installed", {
  skip_if_not_installed("ActivePathways")
  skip_if_not_installed("msigdbr")
  skip_if_not_installed("clusterProfiler")
  p <- make_integration_project()
  diffs <- build_integration_diff_bundles(p)
  b <- run_integration(p, method = "active_pathways",
                       experiments = c("proteo", "rna"),
                       diff_bundles = diffs,
                       database = "hallmark", organism = "Hs")
  expect_true(is_analysis_bundle(b))
  df <- b$results$integration_df
  check_integration_result_schema(df)
  # If any pathways come back, direction should be shared/unique.
  if (nrow(df) > 0L) {
    expect_true(all(df$direction %in% c("shared", "unique")))
  }
})

# ---- plot_integration -------------------------------------------------

test_that("plot_integration scatter view works for correlation", {
  p <- make_integration_project()
  b <- run_integration(p, method = "correlation")
  gg <- plot_integration(b, view = "scatter")
  expect_s3_class(gg, "ggplot")
})

test_that("plot_integration dual_volcano + quadrant work for concordance", {
  p <- make_integration_project()
  diffs <- build_integration_diff_bundles(p)
  b <- run_integration(p, method = "concordance",
                       experiments = c("proteo", "rna"),
                       diff_bundles = diffs)
  expect_s3_class(plot_integration(b, view = "dual_volcano"), "ggplot")
  expect_s3_class(plot_integration(b, view = "quadrant"), "ggplot")
})

test_that("plot_integration rejects view/method mismatches", {
  p <- make_integration_project()
  b <- run_integration(p, method = "correlation")
  expect_error(plot_integration(b, view = "dual_volcano"), "concordance")
  expect_error(plot_integration(b, view = "quadrant"), "concordance")
  expect_error(plot_integration(b, view = "dotplot"), "active_pathways")
})
