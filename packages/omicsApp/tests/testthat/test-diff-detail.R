# A gene selected in the Differential view -- by a row of Top hits or a
# point on the volcano -- is drawn by group in the Selected feature card
# and outlined in the Heatmap card; and the volcano of a large result is
# thinned for the browser without losing a hit.

detail_project <- function() {
  set.seed(7)
  groups <- rep(c("DrugA", "DMSO", "DrugB"), each = 4)
  samp <- paste0("S", seq_along(groups))
  ids <- paste0("P", 1:60)
  m <- matrix(stats::rnorm(60 * length(groups), 12, 0.3), 60,
              dimnames = list(ids, samp))
  m[1:8, groups == "DrugA"] <- m[1:8, groups == "DrugA"] + 2.5
  m[9:14, groups == "DrugA"] <- m[9:14, groups == "DrugA"] - 2.5
  inp <- omicsCore::omics_input(
    m, data.frame(treatment = groups, row.names = samp, stringsAsFactors = FALSE),
    data.frame(feature_id = ids, feature_symbol = paste0("G", 1:60), row.names = ids,
               stringsAsFactors = FALSE),
    omics_type = "proteomics", assay_type = "normalized_intensity")
  omicsCore::omics_project("detail", list(proteomics = inp))
}

run_detail <- function(session) {
  session$setInputs(group_col = "treatment", control = "DMSO", case = "DrugA",
                    method = "limma", fdr_cut = 0.05, fc_cut = 0.263, p_kind = "adj",
                    rerun = 1)
  session$elapse(300)
}

test_that("a row of Top hits selects its feature for the card, the heatmap and the volcano", {
  shiny::testServer(diff_view_server,
                    args = list(current_project = shiny::reactiveVal(detail_project())), {
    run_detail(session)
    # Nothing selected yet: the card asks for a click, and draws nothing.
    expect_null(selected_feature())
    expect_match(output$feature_info$html, "Click a row in Top hits")
    expect_false(detail$has_feature())

    invisible(output$hits)
    hits <- results$hits_df()
    expect_gte(nrow(hits), 10L)
    session$setInputs(hits_rows_selected = 2L)
    id <- hits$feature_id[[2]]
    expect_identical(selected_feature(), id)
    expect_true(detail$has_feature())

    # The card names it, with its log2FC and adjusted p for this comparison.
    html <- output$feature_info$html
    expect_match(html, hits$feature_symbol[[2]], fixed = TRUE)
    expect_match(html, sprintf("log2FC %+.2f", hits$effect[[2]]), fixed = TRUE)
    expect_match(html, "adjusted p", fixed = TRUE)
    expect_match(html, "passes the current thresholds", fixed = TRUE)

    # Its values in the comparison's two groups, the reference first.
    p <- detail$feature_plot()
    built <- ggplot2::ggplot_build(p)
    expect_identical(built$layout$panel_params[[1]]$x$get_labels(), c("DMSO", "DrugA"))
    expect_identical(p$labels$y, "log2 intensity")
    expect_null(p$labels$title)
    expect_match(output$feature_plot$src, "^data:image/png")

    # The heatmap: the hits, the comparison's samples, the feature outlined.
    h <- detail$heatmap_plot()
    expect_identical(length(unique(h$data$y)), min(nrow(hits), HEATMAP_MAX_ROWS))
    expect_identical(length(unique(h$data$x)), 8L)
    txt <- Filter(function(l) inherits(l$geom, "GeomText"), h$layers)[[1]]$data
    expect_identical(txt$face[txt$text == hits$feature_symbol[[2]]], "bold")
    # Drawn once the card is in view (test-lazy-heatmap.R).
    session$setInputs(heatmap_visible = TRUE)
    expect_match(output$heatmap$src, "^data:image/png")
    expect_match(output$heatmap_note$html, "scaled to its own mean")

    # Clear empties the card.
    session$setInputs(clear_feature = 1)
    expect_null(selected_feature())
    expect_match(output$feature_info$html, "Click a row in Top hits")
  })
})

test_that("a click on the volcano selects the feature drawn there", {
  shiny::testServer(diff_view_server,
                    args = list(current_project = shiny::reactiveVal(detail_project())), {
    run_detail(session)
    df <- marked()
    # A feature nowhere near the thresholds.
    i <- which(!df$is_significant & df$adj_p_value > 0.5)[[1]]
    id <- df$feature_id[[i]]
    pt <- list(list(curveNumber = 0, pointNumber = 3, x = df$effect[[i]],
                    y = -log10(df$adj_p_value[[i]])))
    click_id <- paste0("plotly_click-", session$ns("volcano"))
    do.call(session$rootScope()$setInputs,
            stats::setNames(list(as.character(jsonlite::toJSON(pt, auto_unbox = TRUE))),
                            click_id))
    expect_identical(selected_feature(), id)
    # A feature that is not a hit is shown too, and says so.
    expect_match(output$feature_info$html, "does not pass the current thresholds")
  })
})

test_that("nothing significant: the heatmap says so instead of drawing", {
  shiny::testServer(diff_view_server,
                    args = list(current_project = shiny::reactiveVal(detail_project())), {
    run_detail(session)
    session$setInputs(fdr_cut = 1e-300)
    session$elapse(300)
    expect_identical(length(detail$heatmap_hits()$ids), 0L)
    expect_false(detail$has_heatmap())
    expect_match(output$heatmap_note$html, "No feature passes the current thresholds")
  })
})

test_that("the selection follows the feature across comparisons, and is dropped with the result", {
  proj <- detail_project()
  shiny::testServer(diff_view_server,
                    args = list(current_project = shiny::reactiveVal(proj)), {
    session$setInputs(group_col = "treatment", control = "DMSO", case = c("DrugA", "DrugB"),
                      method = "limma", rerun = 1)
    session$elapse(300)
    selected_feature("P3")
    session$setInputs(comparison = "DrugB_vs_DMSO")
    expect_identical(selected_feature(), "P3")
    expect_identical(diff_detail_design(shown_bundle(), active()$input)$levels,
                     c("DMSO", "DrugB"))
    diff_bundle(NULL)
    session$flushReact()
    expect_null(selected_feature())
  })
})

test_that("the comparison's groups: a pair, a weighted contrast, a continuous design", {
  inp <- detail_project()$experiments$proteomics
  pair <- list(params = list(group_col = "treatment", control_group = "DMSO",
                             case_group = "DrugB"))
  expect_identical(diff_detail_design(pair, inp), list(column = "treatment",
                                                      levels = c("DMSO", "DrugB")))
  weighted <- list(params = list(group_col = "treatment", control_group = "DMSO",
                                 contrasts = "(DrugA + DrugB)/2 - DMSO"))
  expect_identical(diff_detail_design(weighted, inp)$levels, c("DMSO", "DrugA", "DrugB"))
  inp$meta_df$dose <- seq_len(nrow(inp$meta_df))
  cont <- list(params = list(analysis_type = "continuous", continuous_col = "dose"))
  expect_identical(diff_detail_design(cont, inp), list(column = "dose", levels = NULL))
  expect_null(diff_detail_design(list(params = list(group_col = "nope")), inp)$column)
})

# ---- thinning the volcano for the browser --------------------------------

big_volcano_bundle <- function(n = 60000L) {
  b <- omicsCore::run_diff(detail_project()$experiments$proteomics, method = "ttest",
                           analysis_type = "group", group_col = "treatment",
                           control_group = "DMSO", case_group = "DrugA")
  df <- b$results$diff_result_df
  set.seed(3)
  df <- df[rep(seq_len(nrow(df)), length.out = n), ]
  df$feature_id <- paste0("F", seq_len(n))
  df$feature_symbol <- paste0("GENE", seq_len(n))
  # A typical RNA-seq screen: most genes unchanged, a few percent moved.
  de <- stats::runif(n) < 0.05
  df$effect <- ifelse(de, stats::rnorm(n, 0, 1.6), stats::rnorm(n, 0, 0.25))
  df$p_value <- 2 * stats::pnorm(-ifelse(de, abs(df$effect) * 3, abs(stats::rnorm(n))))
  df$adj_p_value <- stats::p.adjust(df$p_value, "BH")
  rownames(df) <- NULL
  b$results$diff_result_df <- df
  b
}

test_that("the thinned volcano keeps every hit, label and selection, at a fraction of the size", {
  b <- big_volcano_bundle()
  df <- b$results$diff_result_df
  sig <- df$adj_p_value < 0.05 & abs(df$effect) >= 0.263
  labelled <- volcano_label_ids(df, "adj_p_value", 20L)
  # A grey point in the dense middle of the cloud, selected.
  sel <- df$feature_id[which(!sig & abs(df$effect) < 0.05 & df$adj_p_value > 0.9)[1]]
  shown <- volcano_thin_mask(df, "adj_p_value", sig, c(labelled, sel))

  expect_true(all(shown[sig]))
  expect_true(all(shown[df$feature_id %in% c(labelled, sel)]))
  expect_lt(sum(shown & !sig), sum(!sig) / 10)
  # The axes span what they did.
  expect_true(shown[which.max(df$effect)] && shown[which.min(df$effect)])
  # A result under the limit is sent whole.
  expect_true(all(volcano_thin_mask(df[1:4000, ], "adj_p_value", sig[1:4000])))

  payload <- function(shown) {
    fig <- volcano_figure(b, "adj_p_value", 0.05, 0.263, label_top = TRUE, width = 640,
                          webgl = TRUE, shown = shown, selected = sel)
    w <- suppressWarnings(plotly::plotly_build(fig))
    list(bytes = nchar(htmlwidgets:::toJSON(w$x), type = "bytes"), x = w$x)
  }
  full <- payload(NULL)
  thin <- payload(shown)
  expect_lt(thin$bytes, 1024^2)
  expect_lt(thin$bytes, full$bytes / 4)

  # Every hit is in the figure, by name in its hover text.
  hover <- unlist(lapply(thin$x$data, function(tr) tr$text))
  names_drawn <- sub("<br>.*", "", hover)
  expect_true(all(df$feature_symbol[sig] %in% names_drawn))
  # The selected feature carries its ring, the last trace.
  last <- thin$x$data[[length(thin$x$data)]]
  expect_identical(last$name, "selected")
  # The labels are still laid out.
  expect_gt(length(thin$x$layoutAttrs[[length(thin$x$layoutAttrs)]]$annotations %||%
                     thin$x$layout$annotations), 0L)
})

test_that("the legend says the grey cloud was thinned, and by how much", {
  html <- as.character(volcano_legend(10, 5, ns_shown = 4210, ns_total = 58000))
  expect_match(html, "not significant (thinned for display: 4,210 of 58,000 shown)", fixed = TRUE)
  expect_match(as.character(volcano_legend(10, 5)), ">\\s*not significant\\s*<")
})

test_that("a click finds the nearest feature; the ring is drawn only for a real point", {
  df <- data.frame(feature_id = c("a", "b", "c"), feature_symbol = c("A", "", NA),
                   effect = c(-1, 0, 2), adj_p_value = c(0.01, 0.5, 1e-6),
                   stringsAsFactors = FALSE)
  expect_identical(volcano_click_feature(df, "adj_p_value", 1.98, 6.01), "c")
  expect_identical(volcano_click_feature(df, "adj_p_value", -0.9, 2.1), "a")
  expect_null(volcano_click_feature(df, "adj_p_value", NULL, 1))
  expect_null(volcano_mark_trace(df, "adj_p_value", "zz"))
  expect_null(volcano_mark_trace(df, "adj_p_value", NULL))
  tr <- volcano_mark_trace(df, "adj_p_value", "b")
  expect_identical(tr$text, list("b"))
  expect_equal(tr$y[[1]], -log10(0.5))
})

test_that("the heatmap card grows with its rows, and makes room for a phone's legends", {
  expect_identical(heatmap_height(2), 260)
  expect_gt(heatmap_height(50), heatmap_height(20))
  expect_identical(heatmap_height(20, narrow = TRUE) - heatmap_height(20), 150)
})
