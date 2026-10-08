# What the enrichment dot plot encodes.
#
# ORA used to put the overlap count on x *and* on point size, so one
# number was read twice and a size legend repeated the axis. Colour was
# -log10(p) on a scale whose top was set by the most extreme pathway:
# with one pathway at -log10 p = 248 and the rest between 2 and 15, the
# rest were one grey and the colour told the reader nothing.

ora_df <- function(p = c(1e-248, 4e-16, 5e-12, 2e-9, 2e-9, 4e-3, 1e-2, 3e-2),
                   gene_set_size = c(152, 156, 64, 151, 143, 74, 150, 111),
                   overlap_size = c(151, 37, 21, 28, 27, 11, 16, 12),
                   direction = "up") {
  n <- length(p)
  data.frame(
    database = "hallmark", result_type = "ora", comparison = "B_vs_A",
    pathway_id = paste0("P", seq_len(n)), pathway_name = paste("pathway", seq_len(n)),
    effect = NA_real_, effect_type = NA_character_,
    direction = rep_len(direction, n), p_value = p, adj_p_value = p,
    q_value = NA_real_, gene_set_size = gene_set_size, overlap_size = overlap_size,
    overlap_features = "A", leading_features = NA_character_, source_label = "ora",
    stringsAsFactors = FALSE)
}

gsea_df <- function() {
  data.frame(
    database = "hallmark", result_type = "gsea", comparison = "B_vs_A",
    pathway_id = paste0("P", 1:4), pathway_name = paste("pathway", 1:4),
    effect = c(2.2, 1.7, -1.5, -1.9), effect_type = "nes",
    direction = c("up", "up", "down", "down"),
    p_value = c(1e-9, 1e-4, 1e-2, 1e-3), adj_p_value = c(1e-8, 1e-3, 3e-2, 5e-3),
    q_value = NA_real_, gene_set_size = c(40, 20, 10, 60), overlap_size = NA_real_,
    overlap_features = NA_character_, leading_features = "A", source_label = "gsea",
    stringsAsFactors = FALSE)
}

enrich_bundle <- function(df, type = "ora") {
  new_analysis_bundle("run_enrichment",
                      params = list(type = type, direction = "separate"),
                      results = list(enrich_result_df = df))
}

# The values each aesthetic carries, as plotted.
mapped <- function(p) {
  lapply(p$mapping, function(m) rlang::eval_tidy(m, p$data))
}

colour_scale <- function(p) p$scales$get_scales("colour")

test_that("ORA puts the share of the pathway on x and the overlap on size", {
  df <- ora_df()
  p <- plot_enrichment(enrich_bundle(df), top_n = 20L)
  m <- mapped(p)
  ord <- match(p$data$pathway_id, df$pathway_id)
  expect_equal(m$x, df$overlap_size[ord] / df$gene_set_size[ord])
  expect_equal(m$size, df$overlap_size[ord])
  expect_equal(m$colour, -log10(df$adj_p_value[ord]))
  expect_identical(p$labels$x, "% of pathway genes in the list")
  expect_identical(p$scales$get_scales("size")$name, "genes in list")
})

test_that("no number is mapped to two aesthetics", {
  check <- function(p) {
    m <- mapped(p)
    m <- m[setdiff(names(m), "y")]
    for (i in seq_along(m)) for (j in seq_along(m)) if (i < j) {
      expect_false(isTRUE(all.equal(m[[i]], m[[j]])),
                   info = paste(names(m)[i], "and", names(m)[j]))
    }
  }
  check(plot_enrichment(enrich_bundle(ora_df())))
  fb <- ora_df()
  fb$gene_set_size <- NA_real_
  check(plot_enrichment(enrich_bundle(fb)))
  check(plot_enrichment(enrich_bundle(gsea_df(), "gsea")))
})

test_that("without the pathway size, x falls back to significance and says so", {
  df <- ora_df()
  df$gene_set_size <- NA_real_
  p <- plot_enrichment(enrich_bundle(df))
  m <- mapped(p)
  ord <- match(p$data$pathway_id, df$pathway_id)
  expect_equal(m$x, -log10(df$adj_p_value[ord]))
  expect_identical(p$labels$x, "-log10(adjusted p)")
  # x already is the significance: colouring by it too would repeat it.
  expect_null(p$mapping$colour)
  expect_equal(m$size, df$overlap_size[ord])
  expect_equal(length(ggplot2::layer_grob(p)[[1L]]$x), nrow(df))
})

test_that("one pathway without its size makes the whole panel fall back, not lose a row", {
  df <- ora_df()
  df$gene_set_size[3] <- NA_real_
  p <- plot_enrichment(enrich_bundle(df))
  expect_identical(p$labels$x, "-log10(adjusted p)")
  expect_equal(length(ggplot2::layer_grob(p)[[1L]]$x), nrow(df))
})

test_that("one extreme pathway does not set the top of the colour scale", {
  df <- ora_df()
  p <- plot_enrichment(enrich_bundle(df))
  sc <- colour_scale(p)
  signif <- -log10(df$adj_p_value)
  expect_lt(sc$limits[2L], max(signif))
  # The cap is the largest of the ordinary values, rounded down.
  expect_equal(sc$limits[2L], floor(sort(signif, decreasing = TRUE)[2L]))
  # Squished, not censored: out-of-range values move to the end.
  expect_equal(sc$oob(c(-1, 5, 100), c(0, 10)), c(0, 5, 10))
  expect_equal(utils::tail(sc$breaks, 1L), sc$limits[2L])
  expect_match(utils::tail(sc$labels, 1L), "^\u2265 ")
  expect_identical(sc$name, "-log10(adjusted p)")

  # The extreme pathway is drawn in the top colour, not dropped; the
  # ordinary ones spread over the scale instead of all sitting at grey.
  built <- ggplot2::ggplot_build(p)$data[[1L]]
  ord <- match(p$data$pathway_id, df$pathway_id)
  top <- built$colour[ord == 1L]
  expect_identical(toupper(top), toupper(omics_colors$scale_high))
  expect_gt(length(unique(built$colour)), 5L)
})

test_that("nothing is capped when no pathway stands out", {
  df <- ora_df(p = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2),
               gene_set_size = c(50, 60, 70, 80, 90), overlap_size = c(10, 9, 8, 7, 6))
  sc <- colour_scale(plot_enrichment(enrich_bundle(df)))
  expect_equal(sc$limits, c(2, 6))
  expect_s3_class(sc$labels, "waiver")
})

test_that("a few pathways with one far ahead are still capped", {
  # Three values are too few for quartiles to call the top one an
  # outlier; it is still 20 times the next.
  expect_equal(signif_limits(c(2, 10, 200)), list(limits = c(2, 10), capped = TRUE))
  expect_equal(signif_limits(c(2, 3, 5))$capped, FALSE)
})

test_that("one pathway, or pathways all at the same p, still draw in colour", {
  one <- plot_enrichment(enrich_bundle(ora_df()[1, ]))
  expect_s3_class(ggplot2::ggplot_build(one), "ggplot_built")
  expect_equal(length(ggplot2::layer_grob(one)[[1L]]$x), 1L)
  expect_equal(colour_scale(one)$limits[1L], 0)
  # The size key has its one value rather than no key at all.
  expect_equal(ggplot2::ggplot_build(one)$plot$scales$get_scales("size")$get_breaks(), 151)

  same <- ora_df(p = rep(1e-3, 4), gene_set_size = rep(50, 4), overlap_size = 5:8)
  ps <- plot_enrichment(enrich_bundle(same))
  expect_s3_class(ggplot2::ggplot_build(ps), "ggplot_built")
  expect_equal(colour_scale(ps)$limits, c(0, 3))
  expect_false(anyNA(ggplot2::ggplot_build(ps)$data[[1L]]$colour))
})

test_that("GSEA keeps the NES on x with the zero line, and set size on size", {
  df <- gsea_df()
  p <- plot_enrichment(enrich_bundle(df, "gsea"))
  m <- mapped(p)
  ord <- match(p$data$pathway_id, df$pathway_id)
  expect_equal(m$x, df$effect[ord])
  expect_equal(m$size, df$gene_set_size[ord])
  expect_identical(p$scales$get_scales("size")$name, "set size")
  expect_identical(p$labels$x, "normalized enrichment score (NES)")
  expect_true(any(vapply(p$layers, function(l) inherits(l$geom, "GeomVline"), logical(1))))
})

test_that("the gsea_dot view uses the same significance colours", {
  df <- gsea_df()
  p <- plot_enrichment(enrich_bundle(df, "gsea"), view = "gsea_dot")
  sc <- colour_scale(p)
  # It coloured the raw p, small = red, under the column's name.
  expect_identical(sc$name, "-log10(adjusted p)")
  built <- ggplot2::ggplot_build(p)$data[[1L]]
  ord <- match(p$data$pathway_id, df$pathway_id)
  expect_identical(toupper(built$colour[ord == 1L]), toupper(omics_colors$scale_high))
  expect_identical(toupper(built$colour[ord == 3L]), toupper(omics_colors$scale_low))
})

test_that("ORA's up and down lists share the x axis", {
  df <- ora_df(direction = c("up", "down"))
  p <- plot_enrichment(enrich_bundle(df))
  expect_setequal(as.character(p$data$.list), c("Up-regulated genes", "Down-regulated genes"))
  expect_false(p$facet$params$free$x)
  expect_true(p$facet$params$free$y)
})

test_that("a size key keeps at least two values when it can", {
  # 10..91 asked for three pretty breaks gives 0, 50, 100: one inside.
  expect_gte(length(size_breaks(c(10, 91))), 2L)
  expect_lte(length(size_breaks(c(10, 40))), 3L)
  expect_equal(size_breaks(c(151, 151)), 151)
})
