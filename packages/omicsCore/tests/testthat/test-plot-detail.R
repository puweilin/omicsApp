# One feature by group (plot_feature_expression) and the top hits across
# the samples (plot_heatmap): the two figures the Differential view draws
# for a selected gene.

detail_counts <- function() {
  set.seed(11)
  groups <- rep(c("WT", "KO", "Rescue"), each = 4)
  samp <- paste0("S", seq_along(groups))
  ids <- paste0("ENSG", 1:30)
  m <- matrix(stats::rpois(30 * 12, 200), 30, dimnames = list(ids, samp))
  # A deep library: twice the reads in every gene of S1.
  m[, "S1"] <- m[, "S1"] * 2L
  m[1, groups == "KO"] <- m[1, groups == "KO"] * 4L
  omics_input(
    m, data.frame(group = groups, dose = c(0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3),
                  row.names = samp, stringsAsFactors = FALSE),
    data.frame(feature_id = ids, feature_symbol = c("IL6", paste0("G", 2:30)),
               row.names = ids, stringsAsFactors = FALSE),
    omics_type = "rnaseq", assay_type = "raw_count")
}

detail_prot <- function() {
  set.seed(12)
  groups <- rep(c("Control", "TreatA", "TreatB"), each = 4)
  samp <- paste0("P", seq_along(groups))
  ids <- paste0("P", sprintf("%03d", 1:60))
  m <- matrix(stats::rnorm(60 * 12, 22, 0.3), 60, dimnames = list(ids, samp))
  m[1:20, groups == "TreatA"] <- m[1:20, groups == "TreatA"] + 2
  m[21:30, groups == "TreatA"] <- m[21:30, groups == "TreatA"] - 2
  m[5, 2] <- NA
  omics_input(
    m, data.frame(group = groups, row.names = samp, stringsAsFactors = FALSE),
    data.frame(feature_id = ids, feature_symbol = paste0("GENE", 1:60), row.names = ids,
               stringsAsFactors = FALSE),
    omics_type = "proteomics", assay_type = "normalized_intensity")
}

# ---- plot_feature_expression -------------------------------------------

test_that("counts are drawn as log2(CPM + 1) on whole-library sizes, and the axis says so", {
  inp <- detail_counts()
  p <- plot_feature_expression(inp, "IL6", group_by = "group")
  expect_identical(p$labels$y, "log2(CPM + 1)")
  expect_identical(p$labels$title, "IL6")
  b <- ggplot2::ggplot_build(p)
  pts <- b$data[[2]]
  lib <- colSums(inp$expr_mat)
  want <- log2(inp$expr_mat["ENSG1", ] / lib * 1e6 + 1)
  expect_equal(sort(pts$y), sort(unname(want)), tolerance = 1e-9)
  # Not log2(count + 1): S1's doubled library does not lift it above the
  # other WT samples.
  wt <- want[c("S1", "S2", "S3", "S4")]
  expect_lt(abs(wt[["S1"]] - mean(wt[-1])), 0.5)
})

test_that("group_levels picks the comparison's groups, reference first", {
  inp <- detail_counts()
  p <- plot_feature_expression(inp, "ENSG1", group_by = "group",
                               group_levels = c("WT", "KO"))
  b <- ggplot2::ggplot_build(p)
  expect_identical(b$layout$panel_params[[1]]$x$get_labels(), c("WT", "KO"))
  expect_identical(nrow(b$data[[2]]), 8L)
  # Each group keeps the colour the PCA gives it over all three groups.
  pal <- group_colours(inp$meta_df$group)
  pts <- b$data[[2]]
  cols <- unique(data.frame(x = round(pts$x), colour = pts$colour))
  expect_identical(cols$colour[order(cols$x)], unname(pal[c("WT", "KO")]))
  expect_error(plot_feature_expression(inp, "ENSG1", "group", group_levels = "nope"),
               "None of `group_levels`")
})

test_that("intensities already on a log scale are drawn as they are", {
  inp <- detail_prot()
  p <- plot_feature_expression(inp, "P001", group_by = "group")
  expect_identical(p$labels$y, "log2 intensity")
  y <- ggplot2::ggplot_build(p)$data[[2]]$y
  expect_equal(sort(y), sort(unname(inp$expr_mat["P001", ])))
})

test_that("a numeric group_by is drawn as a scatter with a line", {
  inp <- detail_counts()
  p <- plot_feature_expression(inp, "IL6", group_by = "dose")
  geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_true("GeomSmooth" %in% geoms)
  expect_false("GeomBoxplot" %in% geoms)
  expect_identical(p$labels$x, "dose")
})

test_that("two features sharing a symbol get a panel each", {
  inp <- detail_prot()
  inp$feature_df$feature_symbol[2] <- "GENE1"
  p <- plot_feature_expression(inp, c("P001", "P002"), group_by = "group")
  expect_setequal(levels(p$data$feature), c("GENE1 (P001)", "GENE1 (P002)"))
})

# ---- plot_heatmap ------------------------------------------------------

test_that("the ggplot heatmap draws the comparison's samples, grouped, rows z-scored", {
  inp <- detail_prot()
  b <- run_diff(inp, method = "ttest", analysis_type = "group", group_col = "group",
                control_group = "Control", case_group = "TreatA")
  ids <- paste0("P", sprintf("%03d", 1:30))
  p <- plot_heatmap(b, input = inp, features = ids, group_by = "group",
                    group_levels = c("Control", "TreatA"), engine = "ggplot")
  expect_s3_class(p, "ggplot")
  tiles <- p$data
  # 30 features x 8 samples (TreatB left out).
  expect_identical(nrow(tiles), 30L * 8L)
  # Each row is centred on its mean across the samples shown.
  by_row <- tapply(tiles$value, tiles$y, mean, na.rm = TRUE)
  expect_true(all(abs(by_row) < 1e-9))
  # The missing cell stays missing (drawn grey), not 0.
  expect_identical(sum(is.na(tiles$value)), 1L)
  # Control samples sit left of TreatA, with a gap between.
  bar <- p$layers[[2]]$data
  expect_identical(as.character(bar$group), rep(c("Control", "TreatA"), each = 4))
  expect_true(all(diff(bar$x)[4] > 1))
  # The group bar's colours are the groups' (the PCA's), with a legend.
  pal <- group_colours(inp$meta_df$group)
  expect_identical(unname(p$layers[[2]]$aes_params$fill),
                   unname(pal[as.character(bar$group)]))
  built <- ggplot2::ggplot_build(p)
  expect_true("group" %in% unlist(lapply(built$plot$scales$scales, function(s) s$name)))
})

test_that("the highlighted feature is outlined and named in bold, also when rows go unnamed", {
  inp <- detail_prot()
  b <- run_diff(inp, method = "ttest", analysis_type = "group", group_col = "group",
                control_group = "Control", case_group = "TreatA")
  ids <- paste0("P", sprintf("%03d", 1:30))
  text_layer <- function(p) {
    l <- Filter(function(l) inherits(l$geom, "GeomText"), p$layers)
    l[[1]]$data
  }
  p <- plot_heatmap(b, input = inp, features = ids, group_by = "group",
                    group_levels = c("Control", "TreatA"), highlight = "P007",
                    engine = "ggplot", show_rownames = TRUE)
  txt <- text_layer(p)
  expect_identical(txt$face[txt$text == "GENE7"], "bold")
  expect_true(all(txt$face[txt$text != "GENE7" & txt$text != "group"] == "plain"))
  expect_identical(sum(txt$text %in% paste0("GENE", 1:30)), 30L)
  rects <- Filter(function(l) identical(class(l$geom)[1], "GeomRect"), p$layers)
  expect_length(rects, 1L)

  # Rows unnamed: only the highlighted one keeps its name.
  p2 <- plot_heatmap(b, input = inp, features = ids, group_by = "group",
                     group_levels = c("Control", "TreatA"), highlight = "GENE7",
                     engine = "ggplot", show_rownames = FALSE)
  txt2 <- text_layer(p2)
  expect_setequal(txt2$text, c("GENE7", "group"))
  # Highlighting by symbol finds the same row.
  expect_identical(txt2$face[txt2$text == "GENE7"], "bold")
})

test_that("the heatmap's counts use whole-library sizes", {
  inp <- detail_counts()
  b <- run_diff(inp, method = "ttest", analysis_type = "group", group_col = "group",
                control_group = "WT", case_group = "KO")
  p <- plot_heatmap(b, input = inp, features = c("ENSG1", "ENSG2", "ENSG3"),
                    scale = "none", engine = "ggplot", cluster_rows = FALSE,
                    cluster_cols = FALSE)
  lib <- colSums(inp$expr_mat)
  want <- log2(inp$expr_mat["ENSG1", ] / lib * 1e6 + 1)
  got <- p$data$value[p$data$y == 3]
  expect_equal(unname(got), unname(want), tolerance = 1e-9)
})

test_that("plot_heatmap checks its new arguments by name", {
  inp <- detail_prot()
  expect_error(plot_heatmap(inp, group_levels = "Control"), "`group_levels` needs `group_by`")
  expect_error(plot_heatmap(inp, group_by = "nope", engine = "ggplot"), "`group_by` not found")
  expect_error(plot_heatmap(inp, group_by = "group", group_levels = "nope", engine = "ggplot"),
               "None of `group_levels`")
  expect_error(plot_heatmap(inp, highlight = 1), "highlight")
})

test_that("the ComplexHeatmap heatmap takes the groups and the highlight too", {
  skip_if_not_installed("ComplexHeatmap")
  inp <- detail_prot()
  h <- plot_heatmap(inp, n_top = 10L, group_by = "group", highlight = "P001",
                    engine = "ComplexHeatmap")
  expect_s4_class(h, "Heatmap")
  expect_identical(ncol(h@matrix), 12L)
})
