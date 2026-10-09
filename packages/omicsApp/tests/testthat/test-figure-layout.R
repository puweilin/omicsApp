# How the app lays its figures out (R/ui_helpers.R fit_to_width() and the
# cards' heights): no title over a figure its card already names, a
# phone PCA whose key leaves the points room, and the overlap plot's
# compact layout on a phone.

layout_session <- function(id, width) {
  cd <- list()
  cd[[paste0("output_", id, "_width")]] <- width
  list(clientData = cd, ns = function(x) x)
}

long_group_qc <- function() {
  groups <- c("Vehicle control", "Compound alpha high dose 10 uM 24 h",
              "Compound beta high dose 10 uM 24 h", "Compound gamma low dose 1 uM 24 h",
              "Knockout + rescue construct", "Short")
  set.seed(5)
  ids <- sprintf("S%02d", seq_len(length(groups) * 4L))
  meta <- data.frame(group = rep(groups, each = 4L), row.names = ids)
  m <- matrix(stats::rnorm(200L * length(ids), 20, 1), 200L,
              dimnames = list(paste0("G", 1:200), ids))
  inp <- omicsCore::omics_input(m, meta, data.frame(feature_id = rownames(m),
                                                    feature_symbol = rownames(m)),
                                omics_type = "proteomics",
                                assay_type = "normalized_intensity")
  omicsCore::run_qc(inp, outlier_method = "none", impute_method = "none")
}

# The panel's height, in CSS px, when `p` is drawn at width x height.
panel_height_px <- function(p, width, height) {
  f <- tempfile(fileext = ".png")
  grDevices::png(f, width = width, height = height, res = PLOT_RES)
  on.exit({ grDevices::dev.off(); unlink(f) })
  gt <- ggplot2::ggplotGrob(p)
  row <- unique(gt$layout$t[grepl("^panel", gt$layout$name)])
  others <- grid::convertHeight(sum(gt$heights[-row]), "in", valueOnly = TRUE)
  (height / PLOT_RES - others) * PLOT_RES
}

test_that("a figure's title is dropped in its card, and its subtitle kept", {
  p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) +
    ggplot2::geom_point() +
    ggplot2::labs(title = "PCA of cleaned input", subtitle = "120 features excluded")
  wide <- fit_to_width("pca", p, session = layout_session("pca", 640))
  expect_null(wide$labels$title)
  expect_identical(wide$labels$subtitle, "120 features excluded")
  narrow <- fit_to_width("pca", p, session = layout_session("pca", 293))
  expect_null(narrow$labels$title)
  # A title that is the figure's content stays.
  kept <- fit_to_width("pca", p, session = layout_session("pca", 640), keep_title = TRUE)
  expect_identical(kept$labels$title, "PCA of cleaned input")
  # The figure itself -- what a download saves -- keeps its title.
  expect_identical(p$labels$title, "PCA of cleaned input")
  # A patchwork loses its overall title; its panels keep theirs, which
  # say which panel is which.
  pw <- patchwork::wrap_plots(p, p) +
    patchwork::plot_annotation(title = "Overall", subtitle = "thresholds")
  out <- fit_to_width("x", pw, session = layout_session("x", 640))
  expect_null(out$patches$annotation$title)
  expect_identical(out$patches$annotation$subtitle, "thresholds")
  expect_identical(out[[1]]$labels$title, "PCA of cleaned input")
})

test_that("on a phone, a PCA's group key is two to a row and leaves the points room", {
  qc <- long_group_qc()
  p <- omicsCore::plot_qc(qc, view = "pca", color_by = "group")
  n <- pca_key_entries(p)
  expect_identical(n, 6L)
  h <- pca_plot_px(n, narrow = TRUE)
  expect_gt(h, 360L)
  expect_identical(pca_plot_px(n, narrow = FALSE), 360L)
  phone <- fit_to_width("pca", p, session = layout_session("pca", 293))
  colour <- phone$scales$get_scales("colour")
  expect_equal(colour$guide$params$ncol, 2L)
  labs <- ggplot2::ggplot_build(phone)$plot$scales$get_scales("colour")$get_labels()
  expect_true(all(lengths(strsplit(labs, "\n", fixed = TRUE)) <= 2L))
  # Colour and shape keep the same names and guide, so they stay one key.
  expect_equal(phone$scales$get_scales("shape")$guide$params$ncol, 2L)
  # The figure the card was given -- and the download -- is untouched.
  expect_equal(p$scales$get_scales("colour")$guide$params$ncol, 1)
  # At least 200 px of points on a 390 px phone (a 293 px plot); the key
  # under them used to leave about 70.
  expect_gte(panel_height_px(phone, 293, h), 200)
  # On a desktop the key sits beside the points and they keep the card.
  desk <- fit_to_width("pca", p, session = layout_session("pca", 640))
  expect_equal(desk$scales$get_scales("colour")$guide$params$ncol, 1)
  expect_gte(panel_height_px(desk, 640, 360), 250)
})

test_that("a phone's key keeps one column for one or two entries and colour bars", {
  df <- data.frame(x = 1:4, y = 1:4, g = c("a", "a", "b", "b"), v = 1:4)
  two <- ggplot2::ggplot(df, ggplot2::aes(x, y, colour = g)) + ggplot2::geom_point() +
    ggplot2::scale_colour_discrete(guide = ggplot2::guide_legend(ncol = 1))
  out <- fit_to_width("p", two, session = layout_session("p", 293))
  expect_equal(out$scales$get_scales("colour")$guide$params$ncol, 1)
  cont <- ggplot2::ggplot(df, ggplot2::aes(x, y, colour = v)) + ggplot2::geom_point()
  expect_s3_class(fit_to_width("p", cont, session = layout_session("p", 293)), "ggplot")
  expect_true(is.na(pca_key_entries(cont)))
  expect_identical(pca_plot_px(NA_integer_, TRUE), 360L)
  expect_identical(pca_plot_px(0L, TRUE), 360L)
})

test_that("two-column key names are re-wrapped to half the width", {
  expect_identical(narrow_key_chars(293), 17L)
  expect_gte(narrow_key_chars(100), 8L)
  # A plot that writes its row names into a wide left margin (the
  # heatmap) leaves its key less room.
  p <- ggplot2::ggplot() + omicsCore::theme_omicsCore()
  expect_equal(wide_margin_px(p), 0)
  wide <- p + ggplot2::theme(plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 60))
  expect_equal(wide_margin_px(wide), (60 - 5.5) / 72 * PLOT_RES, tolerance = 0.01)
})

test_that("the overlap card is sized for the compact labels on a phone", {
  cmps <- c("Compound alpha high dose 10 uM 24 h_vs_Vehicle control",
            "Compound beta high dose 10 uM 24 h_vs_Vehicle control",
            "Knockout + rescue construct_vs_Vehicle control",
            "Short_vs_Vehicle control", "Other treatment_vs_Vehicle control")
  desk <- overlap_plot_px(cmps, narrow = FALSE)
  phone <- overlap_plot_px(cmps, narrow = TRUE)
  expect_identical(desk, 270L + label_rows_px(cmps, 10L))
  # Treatment names alone, at most two lines: shorter than the desktop's
  # three-line "<treatment> vs <control>" rows.
  expect_lt(phone, desk)
  expect_identical(phone, as.integer(290L + sum(10L + 15L *
    omicsCore::comparison_label_lines(cmps, compact = TRUE))))
})
