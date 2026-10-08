# "Label top hits" on the volcano: plotly has no label repulsion, and
# twenty labels at one fixed offset piled up wherever the top hits sat
# together -- which on real data is always. The labels are laid out in
# R instead (volcano_annotations() in mod_diff_results.R); these check
# the layout at the card's desktop size and at phone width.

# A figure-shaped list carrying only what volcano_geometry() reads.
fake_volcano_fig <- function(x_half, y_max) {
  list(x = list(layout = list(
    margin = list(t = 38.6, r = 7.3, b = 37.3, l = 43.1),
    xaxis = list(range = c(-x_half, x_half)),
    yaxis = list(range = c(-0.03 * y_max, 1.03 * y_max)))))
}

# The box each label occupies, in px of the plot area: its text at
# ~6.5 px a character, one 13 px line high, from where it is anchored.
label_boxes <- function(ann, geom, x_range) {
  y0 <- geom$y_range[[1]]; y1 <- geom$y_range[[2]]
  do.call(rbind, lapply(ann, function(a) {
    tx <- (a$ax - x_range[[1]]) / diff(x_range) * geom$w
    ty <- (1 - (a$ay - y0) / (y1 - y0)) * geom$h
    w <- nchar(a$text) * 6.5
    data.frame(text = a$text,
               left = if (identical(a$xanchor, "left")) tx else tx - w,
               right = if (identical(a$xanchor, "left")) tx + w else tx,
               top = ty - 6.5, bottom = ty + 6.5)
  }))
}

expect_no_box_overlap <- function(boxes) {
  if (NROW(boxes) < 2L) return(invisible(succeed()))
  pairs <- utils::combn(nrow(boxes), 2L)
  hit <- apply(pairs, 2L, function(p) {
    a <- boxes[p[[1]], ]; b <- boxes[p[[2]], ]
    a$left < b$right && b$left < a$right && a$top < b$bottom && b$top < a$bottom
  })
  expect_false(any(hit), info = paste(
    apply(pairs[, hit, drop = FALSE], 2L, function(p) paste(boxes$text[p], collapse = " / ")),
    collapse = "; "))
}

check_layout <- function(bundle, p_col, width) {
  df <- bundle$results$diff_result_df
  y <- -log10(pmax(df[[p_col]], .Machine$double.xmin))
  fig <- fake_volcano_fig(1.05 * max(abs(df$effect), na.rm = TRUE), max(y, na.rm = TRUE))
  geom <- volcano_geometry(fig, width)
  lab <- volcano_annotations(bundle, 20L, p_col, geom)
  ann <- lab$annotations
  expect_gt(length(ann), 0L)
  expect_lte(length(ann), 20L)

  # Nothing overlaps.
  boxes <- label_boxes(ann, geom, lab$x_range)
  expect_no_box_overlap(boxes)
  # And every label is inside the figure (a few px of margin allowed).
  expect_true(all(boxes$left >= -geom$w * 0.02 & boxes$right <= geom$w * 1.05))
  expect_true(all(boxes$top >= -6 & boxes$bottom <= geom$h + 6))

  # Each leader line ends on its own point: the head is a real feature.
  for (a in ann) {
    row <- df[df$feature_symbol == a$text, , drop = FALSE]
    expect_equal(nrow(row), 1L, info = a$text)
    expect_equal(a$x, row$effect, tolerance = 1e-12, info = a$text)
    expect_equal(a$y, -log10(row[[p_col]]), tolerance = 1e-9, info = a$text)
    # Up to the right of its point, down to the left.
    if (a$x >= 0) {
      expect_identical(a$xanchor, "left")
      expect_gt(a$ax, a$x - 1e-9)
    } else {
      expect_identical(a$xanchor, "right")
      expect_lt(a$ax, a$x + 1e-9)
    }
  }
  # The labels chosen are the most significant ones that fit, in order.
  ranks <- match(vapply(ann, `[[`, "", "text"),
                 df$feature_symbol[order(df[[p_col]])])
  expect_false(is.unsorted(ranks))
  invisible(ann)
}

# Twenty hits crowded at the top of the cloud on both sides, the case
# the old fixed offsets drew as one block of text.
crowded_bundle <- function() {
  set.seed(9)
  n <- 400
  eff <- c(stats::rnorm(n, 0, 0.3), runif(10, 1.2, 1.4), -runif(10, 1.2, 1.4))
  p <- c(stats::runif(n, 0.05, 1), 10^-runif(20, 9, 10))
  sym <- c(paste0("GENE", seq_len(n)), paste0("UPHIT", 1:10), paste0("DOWNHIT", 1:10))
  list(results = list(diff_result_df = data.frame(
    feature_id = sym, feature_symbol = sym, effect = eff, p_value = p,
    adj_p_value = pmin(p * 5, 1), stringsAsFactors = FALSE)))
}

test_that("crowded top hits are labelled without overlap, at desktop and phone width", {
  b <- crowded_bundle()
  for (w in c(620, 300)) {
    ann <- check_layout(b, "adj_p_value", w)
    sides <- vapply(ann, function(a) a$x >= 0, logical(1))
    expect_true(any(sides) && any(!sides), info = w)
  }
})

test_that("the demo's volcano labels do not overlap, by adjusted or raw p", {
  proj <- tutorial_project()
  shiny::testServer(diff_view_server, args = list(current_project = shiny::reactiveVal(proj)), {
    session$setInputs(layer = "proteomics", group_col = "group", control = "Control",
                      case = "TreatA", method = "limma", rerun = 1)
    b <- shown_bundle()
    for (w in c(620, 300)) {
      check_layout(b, "adj_p_value", w)
      check_layout(b, "p_value", w)
    }
  })
})

test_that("the app's figure carries the laid-out labels and the wider x range", {
  proj <- tutorial_project()
  shiny::testServer(diff_view_server, args = list(current_project = shiny::reactiveVal(proj)), {
    session$setInputs(layer = "proteomics", group_col = "group", control = "Control",
                      case = "TreatA", method = "limma", rerun = 1)
    plain <- jsonlite::fromJSON(output$volcano, simplifyVector = FALSE)$x
    session$setInputs(label_top = TRUE)
    fig <- jsonlite::fromJSON(output$volcano, simplifyVector = FALSE)$x
    ann <- fig$layout$annotations
    expect_gt(length(ann), 0L)
    expect_true(all(vapply(ann, function(a) identical(a$axref, "x") &&
                             identical(a$ayref, "y"), logical(1))))
    # Symmetric still, and no narrower than without labels.
    xr <- unlist(fig$layout$xaxis$range)
    expect_equal(xr[[1]], -xr[[2]])
    expect_gte(xr[[2]], unlist(plain$layout$xaxis$range)[[2]] - 1e-9)
  })
})

test_that("labels that cannot fit are left off, least significant first", {
  # A column holds only so many lines: the stacker drops from the end.
  st <- stack_volcano_labels(want = rep(10, 30), rank = 1:30, lo = 8, hi = 100, line = 16)
  expect_identical(sort(st$keep), 1:6)
  expect_true(all(diff(st$y) >= 16 - 1e-9))
  # A label pushed further than allowed from its point is left off too.
  st <- stack_volcano_labels(want = rep(10, 5), rank = 1:5, lo = 8, hi = 300, line = 16,
                             max_shift = 40)
  expect_identical(sort(st$keep), 1:3)
})
