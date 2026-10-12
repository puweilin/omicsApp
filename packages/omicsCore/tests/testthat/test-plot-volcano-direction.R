# The volcano colours its hits by direction, as every other figure in the
# app does (up red, down blue). It used to colour every hit red, so a
# down-regulated gene looked like an up-regulated one.

vd_input <- function(seed = 1) {
  set.seed(seed)
  n <- 6L
  m <- matrix(stats::rnorm(200 * 2 * n, 10, 0.3), 200,
              dimnames = list(paste0("P", 1:200), paste0("s", 1:(2 * n))))
  case <- (n + 1L):(2L * n)
  m[1:20, case] <- m[1:20, case] + 2     # clearly up in B
  m[21:30, case] <- m[21:30, case] - 1   # clearly down in B
  meta <- data.frame(group = rep(c("A", "B"), each = n), row.names = colnames(m))
  omics_input(m, meta, data.frame(feature_id = rownames(m), feature_symbol = rownames(m)),
              omics_type = "proteomics", assay_type = "normalized_intensity")
}

vd_diff <- function() {
  suppressMessages(run_diff(vd_input(), method = "ttest", analysis_type = "group",
                            group_col = "group", control_group = "A", case_group = "B"))
}

colour_of <- function(p, label) {
  sc <- ggplot2::ggplot_build(p)$plot$scales$get_scales("colour")
  unname(sc$map(label))
}

test_that("up and down hits get their own colour and a counted legend entry", {
  b <- vd_diff()
  df <- b$results$diff_result_df
  sig <- !is.na(df$adj_p_value) & df$adj_p_value < 0.05
  n_up <- sum(sig & df$effect > 0)
  n_down <- sum(sig & df$effect < 0)
  expect_gt(n_up, 0L)
  expect_gt(n_down, 0L)

  p <- plot_volcano(b)
  labels <- c(sprintf("up (%d)", n_up), sprintf("down (%d)", n_down), "not significant")
  # Hits first in the legend, in plain words, each with its count.
  expect_identical(levels(p$data$.class), labels)
  expect_identical(sum(p$data$.class == labels[[1]]), n_up)
  expect_identical(sum(p$data$.class == labels[[2]]), n_down)
  expect_identical(colour_of(p, labels[[1]]), omics_colors$up)
  expect_identical(colour_of(p, labels[[2]]), omics_colors$down)
  expect_identical(colour_of(p, labels[[3]]), omics_colors$ns)
  # Every down point really is below zero, every up point above.
  expect_true(all(p$data$effect[p$data$.class == labels[[2]]] < 0))
  expect_true(all(p$data$effect[p$data$.class == labels[[1]]] > 0))
})

test_that("the grey cloud is drawn first, so the hits sit on top of it", {
  p <- plot_volcano(vd_diff(), top_n = 0L)
  points <- Filter(function(l) inherits(l$geom, "GeomPoint"), p$layers)
  expect_length(points, 2L)
  expect_true(all(points[[1]]$data$.class == "not significant"))
  expect_false(any(points[[2]]$data$.class == "not significant"))
  # Smaller and fainter, so the hits stand out.
  expect_lt(points[[1]]$aes_params$size, points[[2]]$aes_params$size)
  expect_lt(points[[1]]$aes_params$alpha, points[[2]]$aes_params$alpha)
})

test_that("a class with no members keeps its legend entry", {
  b <- vd_diff()
  # Nothing reaches |log2FC| >= 1.5 downwards (down is ~ -1).
  p <- plot_volcano(b, effect_threshold = 1.5)
  expect_true("down (0)" %in% levels(p$data$.class))
  expect_true("down (0)" %in% ggplot2::get_guide_data(p, "colour")$.label)
})

test_that("the MA plot colours by direction the same way", {
  b <- vd_diff()
  v <- plot_volcano(b, effect_threshold = 0.5)
  m <- plot_ma(b, effect_threshold = 0.5)
  expect_identical(as.character(m$data$.class), as.character(v$data$.class))
  expect_identical(levels(m$data$.class), levels(v$data$.class))
})

test_that("an F statistic has no direction: its hits are just significant", {
  skip_if_not_installed("limma")
  set.seed(4)
  m <- matrix(stats::rnorm(100 * 9, 10, 0.3), 100,
              dimnames = list(paste0("P", 1:100), paste0("s", 1:9)))
  m[1:10, 4:6] <- m[1:10, 4:6] + 2
  meta <- data.frame(group = rep(c("A", "B", "C"), each = 3), row.names = colnames(m))
  inp <- omics_input(m, meta, data.frame(feature_id = rownames(m), feature_symbol = rownames(m)),
                     omics_type = "proteomics", assay_type = "normalized_intensity")
  b <- suppressMessages(run_diff(inp, method = "limma", analysis_type = "anova",
                                 group_col = "group"))
  expect_identical(unique(b$results$diff_result_df$effect_type), "F_statistic")
  p <- plot_volcano(b, effect_threshold = 1)
  lv <- levels(p$data$.class)
  expect_length(lv, 2L)
  expect_match(lv[[1]], "^significant \\(\\d+\\)$")
  expect_identical(lv[[2]], "not significant")
  expect_false(any(grepl("up|down", lv)))
  # F is never negative: no mirrored axis, no line at -cut.
  vl <- Filter(function(l) inherits(l$geom, "GeomVline"), p$layers)
  expect_identical(vl[[1]]$data$xintercept, 1)
  expect_gte(ggplot2::layer_scales(p)$x$range$range[[1]], 0)
})

test_that("a continuous variable's hits are positive or negative, not up or down", {
  b <- suppressMessages(run_diff_continuous(realistic_input(), method = "lm",
                                            continuous_col = "age"))
  lv <- levels(plot_volcano(b, p_threshold = 0.5)$data$.class)
  expect_match(lv[[1]], "^positive \\(")
  expect_match(lv[[2]], "^negative \\(")
})

test_that("the x axis is symmetric about zero, and drops no point", {
  b <- vd_diff()
  x_range <- function(p) ggplot2::ggplot_build(p)$layout$panel_params[[1]]$x.range
  p <- plot_volcano(b)
  rng <- x_range(p)
  expect_equal(rng[[1]], -rng[[2]])
  expect_gte(rng[[2]], max(abs(b$results$diff_result_df$effect), na.rm = TRUE))
  # A threshold beyond every point still shows both its lines.
  p <- plot_volcano(b, effect_threshold = 5)
  expect_equal(x_range(p)[[1]], -x_range(p)[[2]])
  expect_gte(x_range(p)[[2]], 5)
  built <- ggplot2::ggplot_build(p)
  is_pt <- vapply(p$layers, function(l) inherits(l$geom, "GeomPoint"), logical(1))
  n_drawn <- sum(vapply(built$data[is_pt], nrow, integer(1)))
  expect_identical(n_drawn, sum(!is.na(b$results$diff_result_df$effect) &
                                  !is.na(b$results$diff_result_df$adj_p_value)))
})

test_that("an effect threshold draws dashed lines at both signs", {
  b <- vd_diff()
  vl <- function(p) Filter(function(l) inherits(l$geom, "GeomVline"), p$layers)
  expect_length(vl(plot_volcano(b)), 0L)
  lines <- vl(plot_volcano(b, effect_threshold = 0.263))
  expect_length(lines, 1L)
  expect_setequal(lines[[1]]$data$xintercept, c(-0.263, 0.263))
  expect_identical(lines[[1]]$aes_params$linetype, "dashed")
  # And the horizontal p line stays.
  expect_length(Filter(function(l) inherits(l$geom, "GeomHline"),
                       plot_volcano(b, effect_threshold = 0.263)$layers), 1L)
})

test_that("the caption states the effect cut the way it is applied", {
  p <- plot_volcano(vd_diff(), effect_threshold = 0.263)
  expect_match(p$labels$caption, "adjusted p < 0.05, |mean difference| >= 0.263", fixed = TRUE)
})
