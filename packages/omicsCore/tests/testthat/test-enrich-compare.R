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
