# Enrichment across several comparisons.

three_group_realistic <- function() {
  # G2 carries the G2M signal; G3 -- the control's samples again, with
  # fresh noise -- carries a different set.
  inp <- realistic_input(n_per_group = 4L, signal = "G2M")
  other <- setdiff(names(REAL_GENE_SETS), "G2M")[[1L]]
  set.seed(9)
  ctrl <- inp$expr_mat[, inp$meta_df$group == "G1", drop = FALSE]
  m <- ctrl + matrix(stats::rnorm(length(ctrl), sd = 0.6), nrow(ctrl))
  hit <- inp$feature_df$feature_symbol %in% REAL_GENE_SETS[[other]]
  m[hit, ] <- m[hit, ] + 1.5
  colnames(m) <- paste0("T", seq_len(ncol(m)))
  meta <- inp$meta_df[inp$meta_df$group == "G1", , drop = FALSE]
  rownames(meta) <- colnames(m)
  meta$group <- "G3"
  x <- omics_input(cbind(inp$expr_mat, m), rbind(inp$meta_df, meta),
                   inp$feature_df, omics_type = "proteomics",
                   assay_type = "normalized_intensity")
  list(input = x, other = other)
}

test_that("run_enrichment refuses a stacked multi-contrast table", {
  fx <- three_group_realistic()
  b <- run_diff(fx$input, method = "limma", group_col = "group",
                control_group = "G1", case_group = c("G2", "G3"))
  expect_error(run_enrichment(b), "compare_enrichment")
})

test_that("compare_enrichment enriches each comparison and finds each one's own pathway", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("msigdbr")
  fx <- three_group_realistic()
  b <- run_diff(fx$input, method = "limma", group_col = "group",
                control_group = "G1", case_group = c("G2", "G3"))
  ce <- compare_enrichment(b, type = "ora", database = "hallmark",
                           effect_cutoff = 0.5)
  expect_identical(ce$analysis_name, "compare_enrichment")
  df <- ce$results$enrich_result_df
  expect_setequal(unique(df$comparison), c("G2_vs_G1", "G3_vs_G1"))
  top <- function(cmp) {
    sub <- df[df$comparison == cmp, ]
    sub$pathway_id[which.min(sub$adj_p_value)]
  }
  expect_match(top("G2_vs_G1"), "G2M", fixed = TRUE)
  expect_false(grepl("G2M", top("G3_vs_G1"), fixed = TRUE))
  expect_s3_class(plot_enrichment_comparison(ce), "ggplot")
  expect_error(plot_enrichment_comparison(b), "compare_enrichment")
  expect_error(compare_enrichment(b, comparisons = "nope"), "Not in")
})

test_that("the script repeats the comparison enrichment with its settings", {
  fx <- three_group_realistic()
  b <- run_diff(fx$input, method = "limma", group_col = "group",
                control_group = "G1", case_group = c("G2", "G3"))
  ce <- new_analysis_bundle("compare_enrichment",
                            params = list(type = "ora", database = "hallmark",
                                          direction = "both", p_cutoff = 0.05,
                                          comparison = c("G2_vs_G1", "G3_vs_G1")))
  proj <- omics_project("p", list(proteomics = fx$input))
  proj$bundles <- list(diff = b, enrich_compare = ce)
  txt <- paste(export_script(proj, include_plots = FALSE), collapse = "\n")
  expect_match(txt, "enrich_compare <- compare_enrichment(\n  diff", fixed = TRUE)
  expect_match(txt, 'database  = "hallmark"', fixed = TRUE)
  expect_false(grepl("comparison", sub("(?s).*compare_enrichment\\(", "", txt, perl = TRUE)))
})

# ---- the comparison plot's layout and legends ----------------------------
# A stacked table as compare_enrichment() leaves it: three treatments
# against one control, ORA with up and down tested separately, and one
# pathway far stronger than the rest (as in the realistic project, where
# INFLAMMATORY RESPONSE sits at -log10 p = 248).

compare_fixture <- function(adj_p = c(1e-248, 1e-16, 1e-12, 1e-9, 1e-8, 1e-6,
                                      1e-5, 1e-4, 1e-3),
                            comparisons = c("TreatA_vs_Control", "TreatB_vs_Control",
                                            "TreatC_vs_Control")) {
  n <- length(adj_p)
  df <- data.frame(
    database = "hallmark", result_type = "ora",
    comparison = rep_len(comparisons, n),
    pathway_id = paste0("P", seq_len(n)), pathway_name = paste("pathway", seq_len(n)),
    effect = NA_real_, effect_type = NA_character_,
    direction = rep_len(c("up", "down", "up"), n),
    p_value = adj_p / 10, adj_p_value = adj_p, q_value = NA_real_,
    gene_set_size = 100, overlap_size = 20, overlap_features = "A",
    leading_features = "A", source_label = "ora", stringsAsFactors = FALSE)
  new_analysis_bundle("compare_enrichment",
                      params = list(type = "ora", direction = "separate",
                                    database = "hallmark", comparison = comparisons),
                      results = list(enrich_result_df = df))
}

test_that("the subtitle is two short lines and the title spans the whole plot", {
  # At report size the one-line subtitle ran into the legend title
  # beside the panel ("...tested separately" over "found among").
  p <- plot_enrichment_comparison(compare_fixture())
  lines <- strsplit(p$labels$subtitle, "\n", fixed = TRUE)[[1L]]
  expect_length(lines, 2L)
  expect_true(all(nchar(lines) <= 48L))
  expect_match(lines[[2L]], "tested separately", fixed = TRUE)
  expect_identical(p$theme$plot.title.position, "plot")
  # The legends start at the top of the space beside the panel, under
  # the subtitle, not centred on the panel's height.
  expect_identical(p$theme$legend.justification.right, "top")
  expect_identical(p$theme$legend.box.just, "left")
})

test_that("a control shared by every comparison is said once, in the subtitle", {
  p <- plot_enrichment_comparison(compare_fixture())
  expect_identical(levels(p$data$.col), c("TreatA", "TreatB", "TreatC"))
  expect_match(p$labels$subtitle, "each vs Control", fixed = TRUE)
  # Different controls keep the full "A vs B" under each column.
  p2 <- plot_enrichment_comparison(compare_fixture(
    comparisons = c("B_vs_A", "C_vs_D")))
  expect_identical(levels(p2$data$.col), c("B vs A", "C vs D"))
  expect_false(grepl("each vs", p2$labels$subtitle, fixed = TRUE))
  # Long treatment names are shortened, on one line.
  cols <- comparison_columns(c("Compound alpha high dose 10 uM 24 h_vs_Vehicle",
                               "Short_vs_Vehicle"))
  expect_true(all(nchar(cols$labels) <= 24L))
  expect_false(any(grepl("\n", cols$labels, fixed = TRUE)))
  expect_identical(cols$control, "Vehicle")
})

test_that("the filled/hollow key appears only when something is hollow", {
  # Every pathway significant: a one-entry key ("adjusted p < 0.05: yes")
  # said nothing.
  p <- plot_enrichment_comparison(compare_fixture())
  expect_identical(p$scales$get_scales("shape")$guide, "none")
  # Some hollow: the key says what filled and hollow mean, in words.
  p <- plot_enrichment_comparison(compare_fixture(), p_cutoff = 1e-7)
  sh <- p$scales$get_scales("shape")
  expect_s3_class(sh$guide, "Guide")
  expect_identical(unname(sh$labels[c("yes", "no")]),
                   c("adjusted p < 1e-07", "not significant"))
  expect_true(all(c(16, 1) %in% ggplot2::ggplot_build(p)$data[[1L]]$shape))
})

test_that("point size is the plot's p, capped so one outlier does not set the scale", {
  p <- plot_enrichment_comparison(compare_fixture())
  sz <- p$scales$get_scales("size")
  expect_identical(sz$name, "-log10(adjusted p)")
  # Capped below the outlier at 248, the top key marked "at least".
  lim <- sz$get_limits()
  expect_lt(lim[2L], 50)
  b <- ggplot2::ggplot_build(p)
  labs <- sz$get_labels()
  expect_match(labs[length(labs)], "^\u2265 ")
  # The outlier is squished to the largest dot, not dropped.
  sizes <- b$data[[1L]]$size
  expect_false(anyNA(sizes))
  expect_equal(max(sizes), 6)
  # Mid-range pathways now differ in size (they all came out alike when
  # 248 set the top).
  mid <- sizes[order(-p$data$.neglog)][2:6]
  expect_gt(diff(range(mid)), 1)
  # Raw p when told so.
  p_raw <- plot_enrichment_comparison(compare_fixture(), p_preference = "raw")
  expect_identical(p_raw$scales$get_scales("size")$name, "-log10(p)")
})
