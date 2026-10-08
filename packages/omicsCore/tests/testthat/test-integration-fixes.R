# Integration defects found in the 2026-09 review, one test each. The
# fixtures come from test-integration.R and test-integration-sample-link.R,
# which testthat sources first (alphabetical order), but are rebuilt here
# so this file also runs on its own.

fx_concordance_project <- function() {
  set.seed(11)
  sym <- c("TP53", "EGFR", "MYC", "AKT1", "PTEN", "KRAS", "STAT3", "JUN",
           paste0("G", 1:22))
  groups <- rep(c("ctrl", "case"), each = 4)
  mk <- function(prefix, omics, assay, shift) {
    ids <- paste0(prefix, seq_along(sym))
    samp <- paste0(prefix, "_S", 1:8)
    m <- matrix(rnorm(length(sym) * 8, 8, 0.3), nrow = length(sym),
                dimnames = list(ids, samp))
    m[1:5, groups == "case"] <- m[1:5, groups == "case"] + shift
    omics_input(m,
                data.frame(group = groups, donor = paste0("D", 1:8),
                           row.names = samp, stringsAsFactors = FALSE),
                data.frame(feature_id = ids, feature_symbol = sym,
                           row.names = ids, stringsAsFactors = FALSE),
                omics_type = omics, assay_type = assay)
  }
  omics_project("fx", list(prot = mk("p", "proteomics", "normalized_intensity", 3),
                           rna = mk("r", "rnaseq", "logcpm", 2)))
}

fx_diffs <- function(p) {
  lapply(p$experiments, function(e) {
    run_diff(e, method = "ttest", group_col = "group",
             control_group = "ctrl", case_group = "case")
  })
}

test_that("concordance keeps each layer's own effect, so the effect-pair plot has coordinates", {
  p <- fx_concordance_project()
  d <- fx_diffs(p)
  b <- run_integration(p, "concordance", c("prot", "rna"), diff_bundles = d)
  df <- b$results$integration_df
  expect_true(all(c("effect_a", "effect_b", "p_value_a", "p_value_b",
                    "significant_a", "significant_b") %in% names(df)))
  ra <- d$prot$results$diff_result_df
  expect_equal(df$effect_a[df$feature_symbol == "TP53"],
               ra$effect[ra$feature_symbol == "TP53"])
  expect_equal(df$effect, df$effect_a - df$effect_b)

  # What the plot used to do: y = effect_a - (effect_a - effect_b) went
  # through a recovery that put every point at y = 0.
  g <- ggplot2::ggplot_build(plot_integration(b, view = "effect_pair"))
  is_pts <- vapply(g$data, function(l) nrow(l) == nrow(df), logical(1))
  pts <- g$data[[which(is_pts)[[1L]]]]
  expect_false(all(pts$y == 0))
})

test_that("an older concordance result without per-layer effects draws a message, not fake points", {
  p <- fx_concordance_project()
  b <- run_integration(p, "concordance", c("prot", "rna"), diff_bundles = fx_diffs(p))
  b$results$integration_df$effect_a <- NULL
  b$results$integration_df$effect_b <- NULL
  g <- plot_integration(b, view = "effect_pair")
  expect_true(any(vapply(g$layers, function(l) inherits(l$geom, "GeomText"),
                         logical(1))))
})

test_that("the combined p-value is Fisher on the raw p-values, not on adjusted ones", {
  p <- fx_concordance_project()
  d <- fx_diffs(p)
  b <- run_integration(p, "concordance", c("prot", "rna"), diff_bundles = d,
                       p_preference = "adjusted")
  df <- b$results$integration_df
  expected <- stats::pchisq(-2 * (log(df$p_value_a) + log(df$p_value_b)),
                            df = 4, lower.tail = FALSE)
  expect_equal(df$p_value, expected)
  expect_equal(df$adj_p_value, stats::p.adjust(expected, "BH"))
})

test_that("a p-value that underflowed to zero is kept, not dropped", {
  p <- fx_concordance_project()
  d <- fx_diffs(p)
  d$prot$results$diff_result_df$p_value[1] <- 0
  d$prot$results$diff_result_df$adj_p_value[1] <- 0
  b <- run_integration(p, "concordance", c("prot", "rna"), diff_bundles = d)
  df <- b$results$integration_df
  sym <- d$prot$results$diff_result_df$feature_symbol[1]
  expect_false(is.na(df$p_value[df$feature_symbol == sym]))
  expect_lt(df$p_value[df$feature_symbol == sym], 1e-100)
})

test_that("effect_cutoff and p_preference decide the per-layer hits", {
  p <- fx_concordance_project()
  d <- fx_diffs(p)
  loose <- run_integration(p, "concordance", c("prot", "rna"), diff_bundles = d,
                           p_preference = "raw", effect_cutoff = 0)
  strict <- run_integration(p, "concordance", c("prot", "rna"), diff_bundles = d,
                            p_preference = "raw", effect_cutoff = 2.5)
  expect_gt(sum(loose$results$integration_df$significant_b),
            sum(strict$results$integration_df$significant_b))
  expect_identical(strict$params$effect_cutoff, 2.5)
  expect_error(run_integration(p, "concordance", c("prot", "rna"),
                               diff_bundles = d, effect_cutoff = -1),
               "effect_cutoff")
})

test_that("continuous diffs (positive / negative) still fall into quadrants", {
  p <- fx_concordance_project()
  d <- fx_diffs(p)
  for (tag in names(d)) {
    df <- d[[tag]]$results$diff_result_df
    df$direction <- ifelse(df$effect > 0, "positive", "negative")
    d[[tag]]$results$diff_result_df <- df
  }
  b <- run_integration(p, "concordance", c("prot", "rna"), diff_bundles = d)
  expect_false(all(is.na(b$results$integration_df$quadrant)))
})

test_that("symbols are matched ignoring case, and a second feature of a gene is kept as its own pair", {
  p <- fx_concordance_project()
  d <- fx_diffs(p)
  rb <- d$rna$results$diff_result_df
  rb$feature_symbol <- tolower(rb$feature_symbol)
  # A second, low-abundance row for TP53 in the protein layer.
  ra <- d$prot$results$diff_result_df
  dup <- ra[ra$feature_symbol == "TP53", ]
  dup$feature_id <- "p_dup"
  dup$effect <- -99
  dup$base_mean <- -1
  ra$base_mean <- 10
  d$prot$results$diff_result_df <- rbind(dup, ra)
  d$rna$results$diff_result_df <- rb
  b <- run_integration(p, "concordance", c("prot", "rna"), diff_bundles = d)
  df <- b$results$integration_df
  # Both TP53 proteins are paired with the one TP53 transcript. Keeping
  # only the most abundant (as before feature links) hid the second.
  expect_equal(nrow(df), 31L)
  tp53 <- df[df$feature_symbol == "TP53", ]
  expect_setequal(tp53$feature_id_a, c("p1", "p_dup"))
  expect_identical(unique(tp53$feature_id_b), "r1")
  expect_setequal(tp53$feature_id, c("TP53 (p1)", "TP53 (p_dup)"))
  expect_equal(tp53$effect_a[tp53$feature_id_a == "p_dup"], -99)
  # Lower-cased RNA symbols still matched: every other gene has one row.
  expect_equal(sum(df$feature_symbol != "TP53"), 29L)
  expect_identical(b$params$method_info$feature_pairing$n_a_sharing, 2L)
})

test_that("a multi-contrast diff must be narrowed before it is integrated", {
  p <- fx_concordance_project()
  d <- fx_diffs(p)
  two <- d$prot$results$diff_result_df
  two2 <- two
  two2$comparison <- "other_vs_ctrl"
  d$prot$results$diff_result_df <- rbind(two, two2)
  expect_error(run_integration(p, "concordance", c("prot", "rna"), diff_bundles = d),
               "select_comparison")
})

# ---- sample pairing -----------------------------------------------------

test_that("correlation runs on a donor-column pairing, as the preview promised", {
  p <- fx_concordance_project()
  # Distinct ids, shared donor column: the preview said "from the donor
  # column" while the run failed with "No shared sample IDs".
  expect_identical(sample_pairing_preview(p, "prot", "rna")$source, "donor")
  b <- run_integration(p, "correlation", c("prot", "rna"))
  expect_identical(b$params$method_info$pairing_source, "donor")
  expect_identical(b$params$method_info$n_samples, 8L)
})

test_that("a saved link that does not cover this pair of layers falls through to the donor column", {
  p <- fx_concordance_project()
  p$sample_link <- data.frame(tag = c("other", "other"), sample_id = c("x", "y"),
                              donor_id = c("D1", "D2"), stringsAsFactors = FALSE)
  res <- sample_pairing_preview(p, "prot", "rna")
  expect_identical(res$source, "donor")
  expect_equal(nrow(res$pairs), 8L)
})

test_that("a guessed pairing is never used by a run until it is accepted", {
  ids_a <- paste0("RD00", 1:5, "-C")
  ids_b <- paste0("RD00", 1:5, "_F")
  mk <- function(ids) {
    m <- matrix(rnorm(20), 4, dimnames = list(paste0("F", 1:4), ids))
    omics_input(m, data.frame(g = rep("a", 5), row.names = ids),
                data.frame(feature_id = paste0("F", 1:4),
                           feature_symbol = paste0("S", 1:4),
                           row.names = paste0("F", 1:4)),
                omics_type = "proteomics", assay_type = "normalized_intensity")
  }
  p <- omics_project("g", list(a = mk(ids_a), b = mk(ids_b)))
  expect_identical(sample_pairing_preview(p, "a", "b")$source, "suggested")
  expect_error(run_integration(p, "correlation", c("a", "b"), min_samples = 3L),
               "guess")
})

test_that("a donor with two samples in one layer is paired once and reported", {
  p <- fx_concordance_project()
  meta <- p$experiments$prot$meta_df
  meta$donor[2] <- "D1"  # two protein samples for donor D1
  p$experiments$prot$meta_df <- meta
  b <- run_integration(p, "correlation", c("prot", "rna"))
  expect_identical(b$params$method_info$n_ambiguous_samples, 1L)
  expect_match(b$warnings, "share a donor")
  expect_identical(b$params$method_info$n_samples, 7L)
})

test_that("raw counts are correlated as log-CPM, not as log counts", {
  m <- matrix(c(10, 20, 40, 80,
                100, 200, 400, 800), nrow = 2, byrow = TRUE)
  out <- coerce_to_continuous(m, "raw_count", lib_size = c(1e6, 2e6, 4e6, 8e6))
  # Library size doubles with the counts, so CPM is constant per gene.
  expect_equal(unname(out[1, ]), rep(log2(10 + 1), 4))
  expect_equal(coerce_to_continuous(m, "raw_intensity"), log2(m + 1))
  expect_identical(coerce_to_continuous(m, "logcpm"), m)
})

test_that("a strong Spearman correlation on a dozen samples has a finite, non-zero p-value", {
  set.seed(4)
  n <- 12L
  ids <- paste0("S", seq_len(n))
  x <- matrix(rnorm(5 * n), 5, dimnames = list(paste0("f", 1:5), ids))
  y <- x + matrix(rnorm(5 * n, sd = 0.2), 5)
  y[5, ] <- x[5, ] + 1  # identical ranking: rho = 1
  mk <- function(m, omics) {
    omics_input(m, data.frame(g = rep("a", n), row.names = ids),
                data.frame(feature_id = rownames(m), feature_symbol = paste0("G", 1:5),
                           row.names = rownames(m)),
                omics_type = omics, assay_type = if (omics == "rnaseq") "logcpm"
                                                 else "normalized_intensity")
  }
  p <- omics_project("s", list(a = mk(x, "proteomics"), b = mk(y, "rnaseq")))
  df <- run_integration(p, "correlation", c("a", "b"))$results$integration_df
  expect_true(all(df$effect > 0.9))
  # cor.test()'s Edgeworth "exact" p for n > 9 truncated these to 0.
  expect_true(all(df$p_value > 0))
  expect_true(all(df$p_value < 1e-3))
  expect_equal(df$p_value[df$feature_symbol == "G5"], 2 / factorial(12))
})

# ---- the effect-pair and top-hits views (figure review, 2026-10) -------

fx_hits_bundle <- function() {
  p <- fx_concordance_project()
  run_integration(p, "concordance", c("prot", "rna"), diff_bundles = fx_diffs(p))
}

test_that("the top-hits view lists the hits in both layers, in the table's order", {
  b <- fx_hits_bundle()
  df <- b$results$integration_df
  both <- df$significant_a & df$significant_b
  expect_gt(sum(both), 1L)
  expected <- df$feature_symbol[both][order(df$adj_p_value[both], df$p_value[both])]

  g <- plot_integration(b, view = "top_hits", top_n = 3L)
  # Top row first: the y axis runs bottom-up, so the levels are reversed.
  expect_identical(rev(levels(g$data$.label)), utils::head(expected, 3L))
  # One dot per layer per feature, at that layer's own effect.
  expect_identical(nrow(g$data), 2L * 3L)
  first <- g$data[g$data$.label == expected[[1L]], ]
  row <- df[df$feature_symbol == expected[[1L]], ]
  expect_setequal(first$effect, c(row$effect_a, row$effect_b))
  expect_identical(levels(g$data$layer), c("prot", "rna"))
})

test_that("a feature the layers disagree on is still a top hit, with its line across zero", {
  b <- fx_hits_bundle()
  df <- b$results$integration_df
  i <- which(df$significant_a & df$significant_b)[[1L]]
  df$effect_b[i] <- -abs(df$effect_a[i])
  df$quadrant[i] <- if (df$effect_a[i] > 0) "up_down" else "down_up"
  b$results$integration_df <- df
  g <- plot_integration(b, view = "top_hits")
  expect_true(df$feature_symbol[i] %in% levels(g$data$.label))
})

test_that("the top-hits view says so when no feature is a hit in both layers", {
  b <- fx_hits_bundle()
  b$results$integration_df$significant_b <- FALSE
  g <- plot_integration(b, view = "top_hits")
  txt <- unlist(lapply(g$layers, function(l) l$aes_params$label %||% l$data$label))
  expect_match(paste(txt, collapse = " "), "No feature is a hit in both layers")
})

test_that("the top-hits view of a result without per-layer effects draws the message", {
  b <- fx_hits_bundle()
  b$results$integration_df$effect_a <- NULL
  b$results$integration_df$effect_b <- NULL
  g <- plot_integration(b, view = "top_hits")
  expect_true(any(vapply(g$layers, function(l) inherits(l$geom, "GeomText"),
                         logical(1))))
})

test_that("the effect-pair legend counts each class, and only the top hits are ringed", {
  b <- fx_hits_bundle()
  df <- b$results$integration_df
  both <- df$significant_a & df$significant_b
  g <- plot_integration(b, view = "effect_pair", top_n = 2L)
  labels <- g$scales$get_scales("colour")$labels
  expect_true(any(grepl(sprintf("^up in both \\(%d\\)$",
                                sum(both & df$quadrant == "up_up")), labels)))
  expect_true(any(grepl(sprintf("^not a hit in both \\(%d\\)$", sum(!both)), labels)))
  # Equal axes, symmetric about zero: the diagonal is the line of agreement.
  expect_identical(g$coordinates$ratio, 1)
  expect_equal(g$coordinates$limits$x, -rev(g$coordinates$limits$x))

  rings <- Filter(function(l) inherits(l$geom, "GeomPoint") &&
                    identical(l$aes_params$shape, 21), g$layers)
  expect_length(rings, 1L)
  top <- df$feature_id[both][order(df$adj_p_value[both], df$p_value[both])][1:2]
  expect_setequal(rings[[1L]]$data$feature_id, top)
})

test_that("the per-layer axes say log2FC for group comparisons and 'effect' otherwise", {
  b <- fx_hits_bundle()
  g <- plot_integration(b, view = "effect_pair")
  expect_identical(g$labels$x, "log2FC (prot)")
  b$params$diff_params$rna$analysis_type <- "continuous"
  g <- plot_integration(b, view = "top_hits")
  expect_identical(g$labels$x, "effect")
})
