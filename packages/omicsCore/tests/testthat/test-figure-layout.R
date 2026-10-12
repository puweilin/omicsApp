# Figures laid out for what they hold: the overlap (UpSet) plot on a
# phone, the titles a report reader sees, and report figures as tall as
# their rows.

fl_input <- function(groups, n_per = 4L, n_feat = 160L) {
  set.seed(11)
  ids <- paste0("S", seq_len(length(groups) * n_per))
  meta <- data.frame(group = rep(groups, each = n_per), row.names = ids)
  m <- matrix(stats::rnorm(n_feat * length(ids), 20, 0.5), n_feat,
              dimnames = list(paste0("G", seq_len(n_feat)), ids))
  # Each group shifts its own block of features and a shared one, so the
  # comparisons have hits apart and together.
  for (k in seq_along(groups)[-1L]) {
    own <- (k - 1L) * 12L + seq_len(12L)
    m[c(own, 1:6), meta$group == groups[k]] <- m[c(own, 1:6), meta$group == groups[k]] + 2
  }
  omics_input(m, meta, data.frame(feature_id = rownames(m), feature_symbol = rownames(m)),
              omics_type = "proteomics", assay_type = "normalized_intensity")
}

long_groups <- c("Vehicle control", "Compound alpha high dose 10 uM 24 h",
                 "Compound beta high dose 10 uM 24 h", "Knockout + rescue construct",
                 "Compound gamma low dose 1 uM 24 h", "Short")

overlap_row_labels <- function(ov) {
  dots <- ggplot2::ggplot_build(ov[[2]])
  dots$layout$panel_params[[1]]$y$get_labels()
}

test_that("the compact overlap plot names treatments alone when they share a control", {
  skip_if_not_installed("limma")
  d <- run_diff(fl_input(long_groups), method = "limma", group_col = "group",
                control_group = "Vehicle control", case_group = long_groups[-1])
  full <- plot_diff_overlap(d)
  comp <- plot_diff_overlap(d, compact = TRUE)
  # The full figure keeps "<treatment> vs <control>"; the compact one
  # drops the control from every row and says it once, under the title.
  expect_true(all(grepl("vs Vehicle control \\(\\d+\\)$", overlap_row_labels(full))))
  labs <- overlap_row_labels(comp)
  expect_false(any(grepl("Vehicle", labs)))
  expect_true(all(lengths(strsplit(labs, "\n", fixed = TRUE)) <= 2L))
  expect_true("Short (" %in% substr(labs, 1, 7))
  expect_match(comp$patches$annotation$subtitle, "each vs Vehicle control", fixed = TRUE)
  expect_false(grepl("each vs", full$patches$annotation$subtitle, fixed = TRUE))
  # The title is the figure's, over both panels, not the bar panel's.
  expect_identical(comp$patches$annotation$title, "Shared hits between comparisons")
  expect_null(comp[[1]]$labels$title)
})

test_that("the compact overlap plot draws the six largest combinations and says so", {
  skip_if_not_installed("limma")
  groups <- long_groups[1:4]
  d <- run_diff(fl_input(groups), method = "limma", group_col = "group",
                control_group = groups[1], contrasts = pairwise_contrasts(groups))
  expect_length(diff_comparisons(d), 6L)
  sets <- diff_hit_sets(d)
  ids <- unique(unlist(sets))
  n_combos <- length(unique(vapply(ids, function(i)
    paste(as.integer(vapply(sets, function(s) i %in% s, logical(1))), collapse = ""), "")))
  skip_if(n_combos <= 6L, "the fixture has too few combinations")
  comp <- plot_diff_overlap(d, compact = TRUE)
  bars <- ggplot2::ggplot_build(comp[[1]])$data[[1]]
  expect_equal(nrow(bars), 6L)
  expect_match(comp$patches$annotation$subtitle,
               sprintf("top 6 of %d combinations", n_combos), fixed = TRUE)
  # No shared control: both sides stay, each over at most two lines.
  labs <- overlap_row_labels(comp)
  expect_true(all(grepl("vs ", labs, fixed = TRUE)))
  expect_true(all(lengths(strsplit(labs, "\n", fixed = TRUE)) <= 4L))
  # The full figure draws every combination up to 20 and adds no note
  # while it draws them all.
  full <- plot_diff_overlap(d)
  expect_equal(nrow(ggplot2::ggplot_build(full[[1]])$data[[1]]), min(20L, n_combos))
  if (n_combos <= 20L) expect_false(grepl("combinations", full$patches$annotation$subtitle))
  # comparison_label_lines() sizes a card for these labels.
  # (The rows run bottom to top, the comparisons top to bottom.)
  expect_identical(comparison_label_lines(diff_comparisons(d), compact = TRUE),
                   rev(lengths(strsplit(labs, "\n", fixed = TRUE))))
})

test_that("bar counts are shortened only where columns are narrow", {
  expect_identical(overlap_count_label(c(7L, 120L, 1234L, 15400L)),
                   c("7", "120", "1234", "15400"))
  expect_identical(overlap_count_label(c(7L, 120L, 1234L, 15400L), short = TRUE),
                   c("7", "120", "1.2k", "15k"))
})

test_that("compact row labels are at most two lines a side", {
  cmps <- c("Compound alpha high dose 10 uM 24 h_vs_Vehicle control",
            "Short_vs_Vehicle control")
  expect_identical(comparison_label_lines(cmps, compact = TRUE), c(2L, 1L))
  pw <- c("Compound beta high dose 10 uM 24 h_vs_Compound alpha high dose 10 uM 24 h",
          "Short_vs_Vehicle control")
  expect_identical(comparison_label_lines(pw, compact = TRUE), c(4L, 3L))
  # Without compact, as before: the lines wrap_comparison() uses.
  expect_identical(comparison_label_lines(cmps), c(3L, 2L))
})

test_that("figure titles are in plain words", {
  skip_if_not_installed("limma")
  inp <- fl_input(c("Control", "A", "B"))
  q <- run_qc(inp, outlier_method = "none", impute_method = "none")
  expect_identical(plot_qc(q, view = "pca", color_by = "group")$labels$title,
                   "Samples on the first two principal components")
  d <- run_diff(inp, method = "limma", group_col = "group", control_group = "Control",
                case_group = c("A", "B"))
  d1 <- select_comparison(d, diff_comparisons(d)[1])
  expect_identical(plot_volcano(d1)$labels$title, "Volcano plot")
  expect_identical(plot_pca(inp, color_by = "group")$labels$title,
                   "Samples on the first two principal components")
})

# The PNGs a report embeds, as width and height in pixels, in order.
report_png_sizes <- function(project) {
  path <- tempfile(fileext = ".html")
  suppressMessages(export_report(project, path, format = "html"))
  html <- paste(readLines(path, warn = FALSE), collapse = "\n")
  imgs <- regmatches(html, gregexpr("data:image/png;base64,[A-Za-z0-9+/=]+", html))[[1]]
  # The IHDR chunk: width and height as big-endian 32-bit integers.
  be <- function(b) sum(as.integer(b) * 256^(3:0))
  raws <- lapply(imgs, function(s) jsonlite::base64_dec(sub("data:image/png;base64,", "", s)))
  data.frame(w = vapply(raws, function(r) be(r[17:20]), numeric(1)),
             h = vapply(raws, function(r) be(r[21:24]), numeric(1)))
}

test_that("report figures are as tall as their rows", {
  skip_if_not_installed("rmarkdown")
  skip_if_not(rmarkdown::pandoc_available(), "pandoc unavailable")
  skip_if_not_installed("limma")
  skip_if_not_installed("jsonlite")
  sizes <- lapply(list(three = long_groups[c(1, 6, 4)], six = long_groups), function(g) {
    inp <- fl_input(g)
    d <- run_diff(inp, method = "limma", group_col = "group", control_group = g[1],
                  case_group = g[-1])
    p <- omics_project("Heights", experiments = list(proteomics = inp))
    p$bundles <- list(diff = d)
    report_png_sizes(p)
  })
  # In order: the bars per comparison, their overlap, then one volcano a
  # comparison.
  expect_equal(nrow(sizes$three), 2L + 2L)
  expect_equal(nrow(sizes$six), 2L + 5L)
  # Five long comparison names take more height than two short ones, in
  # both figures that list them; the volcanoes keep one size.
  expect_gt(sizes$six[1, "h"], sizes$three[1, "h"])
  expect_gt(sizes$six[2, "h"], sizes$three[2, "h"])
  expect_equal(unique(sizes$three[-(1:2), "h"]), unique(sizes$six[-(1:2), "h"]))
  # The width is the report's, whatever the height.
  expect_length(unique(c(sizes$three[, "w"], sizes$six[, "w"])), 1L)
})
