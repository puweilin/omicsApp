# The dot plot maps point size to a column that GSEA does not have.
#
# standardize_enrich_results() fills `overlap_size` from `Count` or
# `GeneRatio`; fgsea emits neither, so for a GSEA bundle the column is
# entirely NA. Mapping size to it gave every point an NA size, and
# geom_point drops those -- the panel drew its axes, its facet strips and
# its legend, and not one dot.
#
# That is the failure worth a test: it does not error, and an empty
# enrichment panel is exactly what a result with no significant pathways
# looks like. The bug hides as an answer.

drawn_points <- function(df, p_col = "p_value") {
  g <- ggplot2::layer_grob(plot_enrich_dot(df, p_col))[[1L]]
  if (inherits(g, "zeroGrob")) 0L else length(g$x)
}

enrich_df <- function(overlap_size, gene_set_size = c(120, 80, 200, 60, 45)) {
  data.frame(
    database      = "GO_BP",
    pathway_id    = paste0("GO:", 1:5),
    pathway_name  = paste("pathway", 1:5),
    p_value       = c(1e-4, 2e-3, 8e-3, 2e-2, 4e-2),
    adj_p_value   = c(2e-3, 3e-2, 9e-2, 2e-1, 4e-1),
    effect        = c(2.1, -1.8, 1.5, -1.2, 1.1),
    gene_set_size = gene_set_size,
    overlap_size  = overlap_size,
    direction     = c("up", "down", "up", "down", "up"),
    stringsAsFactors = FALSE
  )
}

test_that("GSEA results plot every point despite having no overlap size", {
  expect_equal(drawn_points(enrich_df(rep(NA_real_, 5))), 5L)
})

test_that("ORA results still size by overlap", {
  expect_equal(drawn_points(enrich_df(c(15, 9, 22, 7, 5))), 5L)
})

test_that("points are drawn even when neither size column is usable", {
  expect_equal(
    drawn_points(enrich_df(rep(NA_real_, 5), rep(NA_real_, 5))),
    5L
  )
})

size_scale_name <- function(p) {
  s <- Filter(function(x) "size" %in% x$aesthetics, p$scales$scales)
  if (length(s) == 0L) NA_character_ else s[[1L]]$name
}

test_that("the size legend is named for the column actually mapped", {
  # A GSEA panel labelled "overlap" would be claiming an overlap count
  # that the method never computed.
  expect_identical(size_scale_name(plot_enrich_dot(enrich_df(rep(NA_real_, 5)),
                                                   "p_value")),
                   "set size")
  expect_identical(size_scale_name(plot_enrich_dot(enrich_df(c(15, 9, 22, 7, 5)),
                                                   "p_value")),
                   "overlap")
})

test_that("no size scale is added when nothing can be mapped to it", {
  expect_true(is.na(size_scale_name(
    plot_enrich_dot(enrich_df(rep(NA_real_, 5), rep(NA_real_, 5)), "p_value"))))
})
