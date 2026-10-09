# Hovering (or tapping) a point in the PCA and the integration plots names
# it (R/plot_hover.R).
#
# The events are built the way the browser builds them: the plot is drawn
# at a card's size, Shiny works out where its panel sits in the image
# (the coordmap the browser receives), and a point's data coordinates are
# turned into the pixel the pointer would be on.

hover_event_at <- function(p, x, y, width = 600, height = 360, dx = 0, dy = 0) {
  f <- tempfile(fileext = ".png")
  grDevices::png(f, width = width, height = height, res = PLOT_RES)
  # Measured on the same device it is drawn on, as renderPlot() does.
  coordmap <- tryCatch({
    built <- shiny:::custom_print.ggplot(p)
    shiny:::getGgplotCoordmap(built, width, height, PLOT_RES)
  }, finally = {
    grDevices::dev.off()
    unlink(f)
  })
  panel <- coordmap$panels[[1L]]
  lim <- panel$domain$discrete_limits
  xn <- if (length(lim$x)) match(as.character(x), unlist(lim$x)) else x
  yn <- if (length(lim$y)) match(as.character(y), unlist(lim$y)) else y
  px <- shiny:::scaleCoords(xn, yn, panel)
  px <- list(x = px$x + dx, y = px$y + dy)
  c(panel[c("domain", "range", "log", "mapping")],
    list(x = xn, y = yn, coords_css = px, coords_img = px,
         img_css_ratio = list(x = 1, y = 1)))
}

tip_html <- function(x) paste(as.character(x$html %||% ""), collapse = "")

test_that("a plot's x and y columns are read from its mapping", {
  df <- data.frame(PC1 = 1:3, PC2 = 3:1, grp = c("a", "b", "c"), eff = 1:3)
  col <- "eff"
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$PC1, y = .data[[col]], colour = grp)) +
    ggplot2::geom_point()
  expect_identical(plot_aes_column(p, "x"), "PC1")
  expect_identical(plot_aes_column(p, "y"), "eff")
  expect_identical(plot_aes_column(p, "colour"), "grp")
  p2 <- ggplot2::ggplot(df, ggplot2::aes(x = PC1 * 2, y = PC2)) + ggplot2::geom_point()
  expect_null(plot_aes_column(p2, "x"))
  expect_null(plot_hover_row(p2, hover_event_at(p2, 2, 3)))
})

test_that("hovering a PCA point names the sample, its group and its coordinates", {
  p <- omicsCore::plot_qc(example_qc_bundle(), view = "pca", color_by = "group")
  d <- p$data
  i <- which.max(d$PC1)
  row <- plot_hover_row(p, hover_event_at(p, d$PC1[i], d$PC2[i], dx = 3, dy = -2))
  expect_identical(row$sample_id, d$sample_id[i])
  txt <- pca_hover_text(row, p)
  expect_identical(txt$title, d$sample_id[i])
  expect_identical(txt$rows[["group"]], as.character(d$group[i]))
  expect_identical(txt$rows[["PC1"]], formatC(round(d$PC1[i], 1), format = "f", digits = 1))
  # Twenty pixels off, in the empty top-left corner of the panel: nothing.
  far <- hover_event_at(p, min(d$PC1), max(d$PC2), dx = -20, dy = -20)
  expect_null(plot_hover_row(p, far))
})

test_that("the QC view shows the hovered sample's card and hides it away from the points", {
  shiny::testServer(qc_view_server, args = list(current_project = shiny::reactiveVal(NULL)), {
    session$setInputs(missing_threshold = 0.5, outlier_method = "iqr")
    p <- pca_plot()
    d <- p$data
    i <- which.min(d$PC2)
    session$setInputs(pca_hover = hover_event_at(p, d$PC1[i], d$PC2[i]))
    html <- tip_html(output$pca_tip)
    expect_match(html, "plot-hover-tip", fixed = TRUE)
    expect_match(html, sprintf(">%s<", d$sample_id[i]), fixed = TRUE)
    group_col <- plot_aes_column(p, "colour")
    expect_match(html, sprintf("<td>%s</td>", d[[group_col]][i]), fixed = TRUE)

    session$setInputs(pca_hover = hover_event_at(p, max(d$PC1), max(d$PC2), dx = 25, dy = -25))
    expect_no_match(tip_html(output$pca_tip), "plot-hover-tip", fixed = TRUE)

    # A tap (a click) shows it too, and the pointer leaving clears it.
    session$setInputs(pca_click = hover_event_at(p, d$PC1[i], d$PC2[i], dx = 12))
    expect_match(tip_html(output$pca_tip), d$sample_id[i], fixed = TRUE)
    session$setInputs(pca_hover = NULL)
    expect_no_match(tip_html(output$pca_tip), "plot-hover-tip", fixed = TRUE)

    # A new image (another colouring) resets the tap to NULL: the card goes.
    session$setInputs(pca_click = hover_event_at(p, d$PC1[i], d$PC2[i]))
    expect_match(tip_html(output$pca_tip), d$sample_id[i], fixed = TRUE)
    session$setInputs(pca_click = NULL)
    expect_no_match(tip_html(output$pca_tip), "plot-hover-tip", fixed = TRUE)
  })
})

test_that("the card stays inside the plot, opening towards the larger side", {
  txt <- list(title = "S1", rows = c(group = "A"))
  ev <- function(x, y) list(coords_css = list(x = x, y = y))
  left_top <- as.character(plot_hover_card(txt, ev(40, 30), width = 300, height = 320))
  expect_match(left_top, "left:52px;max-width:244px;top:42px", fixed = TRUE)
  # Right half and lower half: anchored by its right and bottom edges, and
  # no wider than the room left of the pointer.
  right_low <- as.character(plot_hover_card(txt, ev(200, 300), width = 300, height = 320))
  expect_match(right_low, "right:112px;max-width:184px;bottom:32px", fixed = TRUE)
  phone_mid <- as.character(plot_hover_card(txt, ev(160, 100), width = 300, height = 320))
  expect_match(phone_mid, "right:152px;max-width:144px", fixed = TRUE)
  # Hidden from screen readers; the plot's alt text describes the figure.
  ui <- as.character(hover_plot_output("x-pca", height = "360px"))
  expect_match(ui, 'aria-hidden="true"', fixed = TRUE)
  expect_match(ui, 'data-hover-id="x-pca_hover"', fixed = TRUE)
  expect_match(ui, 'data-click-id="x-pca_click"', fixed = TRUE)
})

test_that("hovering the effect-pair scatter names the gene, both log2FCs and its class", {
  b <- example_integration_bundle()
  p <- omicsCore::plot_integration(b, view = "effect_pair", top_n = 6L)
  d <- p$data
  hit <- which(d$.class != "background")[1L]
  bg <- which(d$.class == "background")[1L]
  row <- plot_hover_row(p, hover_event_at(p, d$effect_a[hit], d$effect_b[hit]))
  txt <- effect_pair_hover_text(row, p)
  # The demo pairs some genes with two proteins: those carry the pair's id.
  name <- if (sum(d$feature_symbol == d$feature_symbol[hit]) > 1L) {
    sprintf("%s (%s)", d$feature_symbol[hit], d$feature_id[hit])
  } else d$feature_symbol[hit]
  expect_identical(txt$title, name)
  exps <- b$params$experiments
  expect_identical(names(txt$rows), paste0("log2FC (", exps, ")"))
  expect_identical(unname(txt$rows[[1L]]), sprintf("%.2f", d$effect_a[hit]))
  expect_identical(unname(txt$rows[[2L]]), sprintf("%.2f", d$effect_b[hit]))
  word <- c(up_up = "up in both", down_down = "down in both",
            up_down = sprintf("%s up, %s down", exps[1], exps[2]),
            down_up = sprintf("%s down, %s up", exps[1], exps[2]))[[as.character(d$.class[hit])]]
  expect_identical(txt$note, paste("hit in both layers ·", word))

  row <- plot_hover_row(p, hover_event_at(p, d$effect_a[bg], d$effect_b[bg]))
  expect_identical(effect_pair_hover_text(row, p)$note, "not a hit in both layers")
})

test_that("the integration view serves the scatter's card", {
  shiny::testServer(integration_view_server, args = list(
    current_project = shiny::reactiveVal(NULL), diff_bundle = shiny::reactiveVal(NULL)), {
    session$setInputs(rerun = 0)
    p <- omicsCore::plot_integration(example_integration_bundle(), view = "effect_pair",
                                     top_n = 6L)
    d <- p$data
    i <- nrow(d)
    session$setInputs(scatter_hover = hover_event_at(p, d$effect_a[i], d$effect_b[i]))
    expect_match(tip_html(output$scatter_tip), sprintf(">%s", d$feature_symbol[i]), fixed = TRUE)
    lim <- max(abs(c(d$effect_a, d$effect_b)), na.rm = TRUE)
    session$setInputs(scatter_hover = hover_event_at(p, -lim, lim, dx = 6, dy = 6))
    expect_no_match(tip_html(output$scatter_tip), "plot-hover-tip", fixed = TRUE)
  })
})

test_that("hovering either dot of a top hit gives the gene's effect in both layers", {
  b <- example_integration_bundle()
  p <- omicsCore::plot_integration(b, view = "top_hits", top_n = 12L)
  d <- p$data
  i <- nrow(d)
  row <- plot_hover_row(p, hover_event_at(p, d$effect[i], d$.label[i]))
  txt <- top_hits_hover_text(row, p)
  expect_identical(txt$title, as.character(d$.label[i]))
  mine <- d[d$.label == d$.label[i], ]
  expect_identical(unname(txt$rows), sprintf("%.2f", mine$effect))
  expect_identical(names(txt$rows), sprintf("log2FC (%s)", mine$layer))
  # Between two rows, off both dots: nothing.
  expect_null(plot_hover_row(p, hover_event_at(p, max(d$effect), d$.label[i], dx = 40)))
})

hover_integration_bundle <- function(method, df, ...) {
  omicsCore::new_analysis_bundle(
    "run_integration",
    params = list(method = method, experiments = c("proteomics", "rnaseq"), ...),
    results = list(integration_df = df))
}

hover_integration_df <- function(n, ...) {
  data.frame(
    feature_id = paste0("f", seq_len(n)), feature_symbol = paste0("G", seq_len(n)),
    result_type = "integration", experiments = "proteomics|rnaseq", comparison = "B_vs_A",
    effect = seq(-0.9, 0.9, length.out = n), effect_type = "spearman_r",
    statistic = 0, statistic_type = "none",
    p_value = 10^-seq(1, 6, length.out = n), adj_p_value = 10^-seq(0.5, 5, length.out = n),
    direction = "up", quadrant = NA_character_,
    is_significant = rep(c(FALSE, TRUE), length.out = n), source_label = "x",
    ..., stringsAsFactors = FALSE)
}

test_that("hovering the correlation scatter gives the gene, its r and adjusted p", {
  b <- hover_integration_bundle("correlation", hover_integration_df(10))
  p <- omicsCore::plot_integration(b, view = "scatter")
  d <- p$data
  i <- which(d$feature_symbol == "G7")
  row <- plot_hover_row(p, hover_event_at(p, d$effect[i], d$.neglog10p[i]))
  txt <- correlation_hover_text(row, p)
  expect_identical(txt$title, "G7")
  expect_identical(txt$rows[["Spearman r"]], sprintf("%.2f", d$effect[i]))
  expect_identical(txt$rows[["adjusted p"]], formatC(d$adj_p_value[i], format = "e", digits = 1))
})

test_that("hovering an ActivePathways dot gives the pathway in full and its adjusted p", {
  df <- hover_integration_df(5, evidence = c("shared", "unique", "combined", "shared", "unique"))
  df$feature_symbol <- paste0("HALLMARK_A_VERY_LONG_PATHWAY_NAME_NUMBER_", 1:5)
  df$effect <- -log10(df$adj_p_value)
  df$direction <- c("up", "down", "mixed", "up", "down")
  p <- omicsCore::plot_integration(hover_integration_bundle("active_pathways", df),
                                   view = "dotplot")
  d <- p$data
  i <- which(d$feature_symbol == df$feature_symbol[3])
  row <- plot_hover_row(p, hover_event_at(p, d$effect[i], d$.label[i]))
  txt <- active_pathways_hover_text(row, p)
  # In the axis's words -- no collection prefix, no underscores -- but
  # whole where the axis shortens it.
  expect_identical(txt$title, "A VERY LONG PATHWAY NAME NUMBER 3")
  # The axis wraps long names over lines rather than cutting them.
  axis <- gsub("\n", " ", sub("…$|\\.\\.\\.$", "", as.character(d$.label[i])), fixed = TRUE)
  expect_true(startsWith(txt$title, axis))
  expect_identical(hover_pathway_name("REACTOME_CELL_CYCLE"), "CELL CYCLE")
  expect_identical(hover_pathway_name("my own set"), "my own set")
  expect_identical(txt$rows[["adjusted p"]], hover_p(df$adj_p_value[3]))
  expect_identical(txt$rows[["found by"]], "only combined")
  expect_identical(txt$rows[["direction"]], "mixed / layers disagree")
})

test_that("p values and effects read as in the tables", {
  expect_identical(hover_p(0.0123), "0.012")
  expect_identical(hover_p(3.21e-5), "3.2e-05")
  expect_identical(hover_p(NA), "—")
  expect_identical(hover_num(-0.026, 1L), "0.0")
  expect_identical(hover_num(1.234), "1.23")
})
