# The enrichment bar view.
#
# Its fill was the database, in ggplot's default hue -- salmon for the
# usual single database -- with the legend hidden, so the colour meant
# nothing. And bar length is -log10(p): one pathway at 248 left every
# other bar a sliver at the left of the panel.

layer_of <- function(p, geom) {
  Filter(function(l) inherits(l$geom, geom), p$layers)
}

built_layer <- function(p, geom) {
  b <- ggplot2::ggplot_build(p)
  i <- which(vapply(p$layers, function(l) inherits(l$geom, geom), logical(1)))
  b$data[[i[1L]]]
}

bars <- function(p) {
  d <- built_layer(p, "GeomCol")
  d[order(d$PANEL, d$y), ]
}

test_that("ORA bars are coloured by gene list, without a legend", {
  df <- ora_df(direction = c("up", "down"))
  p <- plot_enrichment(enrich_bundle(df), view = "bar")
  expect_identical(sort(unique(p$data$.fill[p$data$direction == "up"])), omics_colors$up)
  expect_identical(sort(unique(p$data$.fill[p$data$direction == "down"])), omics_colors$down)
  b <- built_layer(p, "GeomCol")
  expect_setequal(toupper(b$fill), toupper(c(omics_colors$up, omics_colors$down)))
  # The facet strips name the list; a legend would repeat them.
  expect_identical(p$scales$get_scales("fill")$guide, "none")
})

test_that("GSEA bars are coloured by the sign of the NES, with a legend on top", {
  df <- gsea_df()
  p <- plot_enrichment(enrich_bundle(df, "gsea"), view = "bar")
  ord <- match(p$data$pathway_id, df$pathway_id)
  expect_identical(p$data$.fill,
                   ifelse(df$effect[ord] > 0, omics_colors$up, omics_colors$down))
  sc <- p$scales$get_scales("fill")
  expect_identical(sc$guide, "legend")
  expect_identical(sc$labels, c("NES > 0 (up)", "NES < 0 (down)"))
  # fit_to_width() leaves a legend on top where it is.
  expect_identical(p$theme$legend.position, "top")
})

test_that("results with no direction take one neutral colour", {
  df <- ora_df(direction = NA_character_)
  p <- plot_enrichment(enrich_bundle(df), view = "bar")
  expect_identical(unique(p$data$.fill), omics_colors$fg_dark)
  expect_identical(p$scales$get_scales("fill")$guide, "none")
})

test_that("an outlying bar is cut at the cap, marked, and labelled with its value", {
  df <- ora_df()
  signif <- -log10(df$adj_p_value)
  top <- sort(signif, decreasing = TRUE)[2L]
  p <- plot_enrichment(enrich_bundle(df), view = "bar")
  ord <- match(p$data$pathway_id, df$pathway_id)
  # The outlier is drawn to the longest ordinary bar; the others in full.
  expect_equal(p$data$.bar, pmin(signif[ord], top))
  expect_equal(max(built_layer(p, "GeomCol")$xmax), top)
  expect_identical(p$labels$x, "-log10(adjusted p)")

  txt <- built_layer(p, "GeomText")
  expect_identical(txt$label, "248")
  expect_equal(txt$x, top)

  # Two white slashes across that bar, near its end.
  seg <- built_layer(p, "GeomSegment")
  expect_equal(nrow(seg), 2L)
  expect_true(all(seg$colour == "white"))
  expect_true(all(seg$x > 0.75 * top & seg$xend < top))
  # The outlier is the most significant of eight: the top row.
  expect_true(all(abs((seg$y + seg$yend) / 2 - 8) < 1e-8))
  expect_equal(txt$y, 8)
})

test_that("the break mark is drawn in the panel of the bar it cuts", {
  # The outlier in the down list, the second panel: its y is its rank
  # among that panel's pathways, not among all of them.
  df <- ora_df(direction = c("down", "up", "down", "up", "down", "up", "up", "up"))
  p <- plot_enrichment(enrich_bundle(df), view = "bar")
  layout <- ggplot2::ggplot_build(p)$layout$layout
  down_panel <- layout$PANEL[layout$.list == "Down-regulated genes"]
  # In that panel, the outlier is the most significant of three, so the
  # top row: y = 3.
  seg <- built_layer(p, "GeomSegment")
  txt <- built_layer(p, "GeomText")
  expect_true(all(seg$PANEL == down_panel))
  expect_true(all(txt$PANEL == down_panel))
  expect_equal(txt$y, 3)
  expect_true(all(abs((seg$y + seg$yend) / 2 - 3) < 1e-8))
  col <- built_layer(p, "GeomCol")
  expect_equal(max(col$y[col$PANEL == down_panel]), 3)
})

test_that("no bar is cut when nothing stands out", {
  df <- ora_df(p = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2),
               gene_set_size = rep(50, 5), overlap_size = c(10, 9, 8, 7, 6))
  p <- plot_enrichment(enrich_bundle(df), view = "bar")
  expect_equal(sort(p$data$.bar), 2:6)
  expect_length(layer_of(p, "GeomSegment"), 0L)
  expect_length(layer_of(p, "GeomText"), 0L)
})

test_that("a bar only a little past the fence is drawn in full", {
  # 9 is beyond the quartile fence of these values, but not three times
  # the next: cutting it would make the plot harder to read, not easier.
  df <- ora_df(p = 10^-c(9, 4.1, 4, 1.7, 1.7), gene_set_size = rep(50, 5),
               overlap_size = c(10, 9, 8, 7, 6))
  p <- plot_enrichment(enrich_bundle(df), view = "bar")
  expect_equal(max(p$data$.bar), 9)
  expect_length(layer_of(p, "GeomSegment"), 0L)
})

test_that("a single pathway draws one full bar", {
  p <- plot_enrichment(enrich_bundle(ora_df()[1, ]), view = "bar")
  col <- built_layer(p, "GeomCol")
  expect_equal(nrow(col), 1L)
  expect_equal(col$xmax, 248, tolerance = 1e-6)
  expect_length(layer_of(p, "GeomSegment"), 0L)
})

test_that("the bar view draws at phone width, the value label unclipped", {
  phone <- ggplot2::theme(text = ggplot2::element_text(size = 9),
                          legend.position = "bottom", legend.direction = "horizontal",
                          legend.box = "vertical")
  for (b in list(enrich_bundle(ora_df(direction = c("up", "down"))),
                 enrich_bundle(gsea_df(), "gsea"))) {
    p <- plot_enrichment(b, view = "bar")
    expect_identical(p$coordinates$clip, "off")
    path <- tempfile(fileext = ".png")
    grDevices::png(path, width = 300, height = 420, res = 96)
    expect_no_error(print(p + phone))
    grDevices::dev.off()
    unlink(path)
  }
})
