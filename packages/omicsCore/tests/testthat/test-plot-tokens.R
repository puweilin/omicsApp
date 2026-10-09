# The sample-group colours: told apart by readers with red-green colour
# blindness, never mistaken for up or down, and the same group the same
# colour in every figure.

# OKLab, and the Machado, Oliveira & Fernandes (2009) simulations at
# severity 1 -- the model the dataviz palette rules are calibrated to.
# Distances are Euclidean in OKLab x 100.
oklab_of <- function(lin) {
  lms <- matrix(c(0.4122214708, 0.5363325363, 0.0514459929,
                  0.2119034982, 0.6806995451, 0.1073969566,
                  0.0883024619, 0.2817188376, 0.6299787005), 3, byrow = TRUE) %*% lin
  lms <- sign(lms) * abs(lms)^(1 / 3)
  drop(matrix(c(0.2104542553, 0.7936177850, -0.0040720468,
                1.9779984951, -2.4285922050, 0.4505937099,
                0.0259040371, 0.7827717662, -0.8086757660), 3, byrow = TRUE) %*% lms)
}
linear_rgb <- function(hex) {
  s <- grDevices::col2rgb(hex)[, 1] / 255
  ifelse(s <= 0.04045, s / 12.92, ((s + 0.055) / 1.055)^2.4)
}
cvd_matrix <- list(
  protan = matrix(c(0.152286, 1.052583, -0.204868,
                    0.114503, 0.786281, 0.099216,
                    -0.003882, -0.048116, 1.051998), 3, byrow = TRUE),
  deutan = matrix(c(0.367322, 0.860646, -0.227968,
                    0.280085, 0.672501, 0.047413,
                    -0.011820, 0.042940, 0.968881), 3, byrow = TRUE))
lab_as_seen <- function(hex, kind = NULL) {
  lin <- linear_rgb(hex)
  if (!is.null(kind)) lin <- pmin(pmax(drop(cvd_matrix[[kind]] %*% lin), 0), 1)
  oklab_of(lin)
}
delta_e <- function(a, b, kind = NULL) {
  100 * sqrt(sum((lab_as_seen(a, kind) - lab_as_seen(b, kind))^2))
}
# Normal vision, and the worse of protanopia and deuteranopia.
separation <- function(a, b) {
  c(normal = delta_e(a, b),
    cvd = min(delta_e(a, b, "protan"), delta_e(a, b, "deutan")))
}

test_that("group_palette gives eight distinct colours, then repeats them", {
  expect_identical(group_palette(0), character(0))
  p8 <- group_palette(8)
  expect_length(p8, 8L)
  expect_true(all(grepl("^#[0-9A-F]{6}$", p8)))
  expect_false(anyDuplicated(p8) > 0L)
  expect_identical(group_palette(3), p8[1:3])
  expect_identical(group_palette(11), c(p8, p8[1:3]))
  expect_error(group_palette(-1), "whole number")
})

test_that("the first six group colours are told apart, every pair, with and without CVD", {
  pal <- group_palette(6)
  for (i in 1:5) for (j in (i + 1):6) {
    s <- separation(pal[i], pal[j])
    expect_gte(s[["cvd"]], 8, label = paste("CVD", pal[i], pal[j]))
    expect_gte(s[["normal"]], 15, label = paste("normal", pal[i], pal[j]))
  }
})

test_that("neighbouring group colours are told apart, also where they repeat", {
  pal <- group_palette(9)  # the ninth is the first again, beside the eighth
  for (i in 1:8) {
    s <- separation(pal[i], pal[i + 1L])
    expect_gte(s[["cvd"]], 8, label = paste("CVD", pal[i], pal[i + 1L]))
    expect_gte(s[["normal"]], 15, label = paste("normal", pal[i], pal[i + 1L]))
  }
})

test_that("no group colour looks like up or down", {
  for (g in group_palette(8)) {
    for (dir in c(omics_colors$up, omics_colors$down)) {
      s <- separation(g, dir)
      expect_gte(s[["normal"]], 20, label = paste(g, "vs", dir))
      expect_gte(s[["cvd"]], 10, label = paste(g, "vs", dir, "under CVD"))
    }
  }
})

test_that("group colours are coloured enough to be colours, and not too light", {
  for (g in group_palette(8)) {
    lab <- lab_as_seen(g)
    expect_gte(sqrt(lab[2]^2 + lab[3]^2), 0.10, label = paste("chroma of", g))
    expect_true(lab[1] >= 0.43 && lab[1] <= 0.77, label = paste("lightness of", g))
  }
})

test_that("colours are dealt reference first, then by the column's order", {
  v <- c("TreatB", "Control", "TreatA", NA, "TreatA")
  expect_identical(group_order(v), c("Control", "TreatA", "TreatB"))
  expect_identical(group_order(v, reference = "TreatA"), c("TreatA", "Control", "TreatB"))
  # A reference that is not a group changes nothing.
  expect_identical(group_order(v, reference = "Vehicle"), group_order(v))
  # A factor's own order, without groups that have no samples.
  f <- factor(c("hi", "lo", "mid"), levels = c("lo", "mid", "hi", "none"))
  expect_identical(group_order(f), c("lo", "mid", "hi"))
  pal <- group_colours(v, reference = "Control")
  expect_identical(unname(pal), group_palette(3))
  expect_identical(names(pal), c("Control", "TreatA", "TreatB"))
})

# A four-group layer whose recorded reference, Vehicle, sorts last: a
# figure that sorted the groups would deal it the last colour.
tok_input <- function(groups = c("Vehicle", "A", "B", "C"), n_per = 3L) {
  set.seed(11)
  ids <- paste0("S", seq_len(length(groups) * n_per))
  meta <- data.frame(group = rep(groups, each = n_per),
                     batch = rep(c("b1", "b2", "b3"), length.out = length(ids)),
                     row.names = ids)
  m <- matrix(stats::rnorm(40 * length(ids), 20, 0.5), 40,
              dimnames = list(paste0("G", 1:40), ids))
  for (k in seq_along(groups)[-1L]) m[1:10, meta$group == groups[k]] <- m[1:10, meta$group == groups[k]] + k
  x <- omics_input(m, meta, data.frame(feature_id = rownames(m), feature_symbol = rownames(m)),
                   omics_type = "proteomics", assay_type = "normalized_intensity")
  set_study_design(x, "group", groups[[1L]])
}

# The colour each group's points or tiles are drawn in.
drawn_colours <- function(df, group, colour = df$colour) {
  u <- unique(data.frame(g = as.character(group), c = colour))
  stats::setNames(u$c, u$g)
}

test_that("a group is the same colour in the PCA, the boxplot and the heatmap", {
  inp <- tok_input()
  want <- stats::setNames(group_palette(4), c("Vehicle", "A", "B", "C"))
  q <- run_qc(inp, outlier_method = "none", impute_method = "none")

  pca <- plot_qc(q, view = "pca", color_by = "group", reference = "Vehicle")
  b <- ggplot2::ggplot_build(pca)
  got <- drawn_colours(b$data[[1]], pca$data$group)
  expect_identical(got[names(want)], want)
  # The legend lists the reference first, with its own shape.
  guide <- ggplot2::get_guide_data(pca, "colour")
  expect_identical(as.character(guide$.value), names(want))
  expect_identical(ggplot2::get_guide_data(pca, "shape")$shape, c(16L, 17L, 15L, 18L))
  # Not one of ggplot's default hues, the first a salmon red.
  expect_false(scales::hue_pal()(4)[1] %in% b$data[[1]]$colour)

  pp <- plot_pca(inp, color_by = "group")  # reference from the layer's design
  got <- drawn_colours(ggplot2::ggplot_build(pp)$data[[1]], pp$data$group)
  expect_identical(got[names(want)], want)

  # The comparison's two groups, its reference first: each keeps its colour.
  fe <- plot_feature_expression(inp, "G1", group_by = "group", group_levels = c("B", "C"))
  pts <- ggplot2::ggplot_build(fe)$data[[2]]
  got <- drawn_colours(pts, c("B", "C")[round(pts$x)])
  expect_identical(got[c("B", "C")], want[c("B", "C")])
  # Every group, the reference first on the axis too.
  fe <- plot_feature_expression(inp, "G1", group_by = "group")
  expect_identical(levels(fe$data$.group), names(want))

  d <- run_diff(inp, method = "ttest", analysis_type = "group", group_col = "group",
                control_group = "Vehicle", case_group = "A")
  hm <- plot_heatmap(d, input = inp, features = paste0("G", 1:12), group_by = "group",
                     group_levels = c("Vehicle", "A"), engine = "ggplot")
  bar <- hm$layers[[2]]
  got <- drawn_colours(bar$data, bar$data$group, unname(bar$aes_params$fill))
  expect_identical(got[c("Vehicle", "A")], want[c("Vehicle", "A")])
  hm_all <- plot_heatmap(inp, n_top = 12, group_by = "group", engine = "ggplot")
  expect_identical(levels(hm_all$layers[[2]]$data$group), names(want))
})

test_that("a factor group column keeps its own order in every figure", {
  inp <- tok_input(c("high", "mid", "low"))
  inp$design <- NULL
  inp$meta_df$group <- factor(inp$meta_df$group, levels = c("low", "mid", "high"))
  want <- stats::setNames(group_palette(3), c("low", "mid", "high"))
  pp <- plot_pca(inp, color_by = "group")
  got <- drawn_colours(ggplot2::ggplot_build(pp)$data[[1]], pp$data$group)
  expect_identical(got[names(want)], want)
  fe <- plot_feature_expression(inp, "G1", group_by = "group")
  expect_identical(levels(fe$data$.group), names(want))
  pts <- ggplot2::ggplot_build(fe)$data[[2]]
  got <- drawn_colours(pts, names(want)[round(pts$x)])
  expect_identical(got[names(want)], want)
  hm <- plot_heatmap(inp, n_top = 10, group_by = "group", engine = "ggplot")
  bar <- hm$layers[[2]]
  expect_identical(levels(bar$data$group), names(want))
  got <- drawn_colours(bar$data, bar$data$group, unname(bar$aes_params$fill))
  expect_identical(got[names(want)], want)
})

test_that("a boxplot coloured by another column uses the group colours for it", {
  inp <- tok_input()
  fe <- plot_feature_expression(inp, "G1", group_by = "group", color_by = "batch")
  pts <- ggplot2::ggplot_build(fe)$data[[2]]
  expect_setequal(unique(pts$colour), group_palette(3))
})

test_that("nine or more groups: colours repeat, and the PCA's shapes tell them apart", {
  groups <- c("Ctrl", paste0("T", sprintf("%02d", 1:9)))
  inp <- tok_input(groups, n_per = 2L)
  q <- run_qc(inp, outlier_method = "none", impute_method = "none")
  p <- plot_qc(q, view = "pca", color_by = "group", reference = "Ctrl")
  b <- ggplot2::ggplot_build(p)
  key <- unique(data.frame(g = p$data$group, colour = b$data[[1]]$colour,
                           shape = b$data[[1]]$shape))
  expect_identical(nrow(key), 10L)
  # No two groups share both colour and shape.
  expect_false(anyDuplicated(paste(key$colour, key$shape)) > 0L)
  expect_identical(key$colour[key$g == "T08"], key$colour[key$g == "Ctrl"])
  expect_false(identical(key$shape[key$g == "T08"], key$shape[key$g == "Ctrl"]))
  expect_identical(group_shapes(10), c(rep(16L, 8L), 17L, 17L))
  # The heatmap and the boxplot draw too.
  expect_s3_class(plot_heatmap(inp, n_top = 10, group_by = "group", engine = "ggplot"), "ggplot")
  expect_s3_class(plot_feature_expression(inp, "G1", group_by = "group"), "ggplot")
})
