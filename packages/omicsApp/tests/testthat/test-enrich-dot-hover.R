# Hovering (or tapping) a dot of the enrichment dot plot names the
# pathway in full, with its gene list, overlap and p (R/plot_hover.R,
# R/mod_enrich_results.R). The plot is faceted by database and gene
# list, so the event must be read in the panel it came from.

# A pathway found in both lists of one database, drawn at the same place
# in the two panels: second of two rows, a quarter of its genes found.
# Only the panel tells the two dots apart.
shared_name <- "REGULATION OF CYSTEINE-TYPE ENDOPEPTIDASE ACTIVITY INVOLVED IN APOPTOTIC PROCESS"
hover_ora_bundle <- function() {
  df <- data.frame(
    database = "go_bp", result_type = "ora", comparison = "Treat_vs_Ctrl",
    pathway_id = c("GO_SHARED", "GO_UP", "GO_SHARED", "GO_DOWN"),
    pathway_name = c(shared_name, "AEROBIC RESPIRATION", shared_name, "NEUTROPHIL DEGRANULATION"),
    effect = NA_real_, effect_type = NA_character_,
    direction = c("up", "up", "down", "down"),
    p_value = c(1e-7, 1e-4, 1e-5, 1e-3), adj_p_value = c(1e-6, 1.234e-3, 1.8e-248, 0.0213),
    q_value = NA_real_,
    gene_set_size = c(40, 30, 20, 50), overlap_size = c(10, 6, 5, 9),
    overlap_features = "A/B", leading_features = NA_character_, source_label = "ora",
    stringsAsFactors = FALSE)
  omicsCore::new_analysis_bundle(
    "run_enrichment", input_info = list(omics_type = "proteomics"),
    params = list(type = "ora", database = "go_bp", organism = "Homo sapiens",
                  direction = "separate", comparison = "Treat_vs_Ctrl"),
    results = list(enrich_result_df = df))
}

hover_gsea_bundle <- function() {
  b <- hover_ora_bundle()
  df <- b$results$enrich_result_df[c(1, 2, 4), ]
  df$result_type <- "gsea"
  df$effect <- c(2.345, -1.5, 1.2)
  df$effect_type <- "nes"
  df$direction <- ifelse(df$effect > 0, "up", "down")
  df$overlap_size <- NA_real_
  b$params$type <- "gsea"
  b$results$enrich_result_df <- df
  b
}

# The event the browser sends for a point of `row` (a row of the plot's
# data), from the coordmap renderPlot() sent with the image: the panel
# whose facet values are the row's.
dot_hover_event <- function(coordmap, row) {
  for (panel in coordmap$panels) {
    pv <- panel$panel_vars
    same <- vapply(names(pv), function(v) {
      col <- panel$mapping[[v]]
      identical(as.character(pv[[v]]), as.character(row[[col]]))
    }, logical(1))
    if (!all(same)) next
    xcol <- plot_aes_column(list(mapping = list(x = rlang::parse_expr(panel$mapping$x))), "x")
    yn <- match(as.character(row$.label), unlist(panel$domain$discrete_limits$y))
    px <- shiny:::scaleCoords(row[[xcol]], yn, panel)
    return(c(panel[c("domain", "range", "log", "mapping")], pv,
             list(x = row[[xcol]], y = yn, coords_css = px, coords_img = px,
                  img_css_ratio = list(x = 1, y = 1))))
  }
  stop("no panel for that row")
}

dot_tip <- function(x) paste(as.character(x$html %||% ""), collapse = "")

test_that("hovering a dot names its pathway, list, overlap and p, in the right panel", {
  b <- hover_ora_bundle()
  shiny::testServer(enrich_view_server, args = list(diff_bundle = shiny::reactiveVal(NULL)), {
    enrich_bundle(b); is_demo(FALSE)
    session$setInputs(show_p = "adjusted", show_cutoff = 0.05)
    session$flushReact()
    p <- results$dot_plot()
    d <- p$data
    cm <- output$dot$coordmap
    # One panel per gene list.
    expect_length(cm$panels, 2L)
    # The shared pathway's dot in the down panel sits where its dot in
    # the up panel does, within each panel: read in either panel's
    # coordinates, the two rows land on the same pixel, and only the
    # panel tells them apart.
    up <- d[d$pathway_id == "GO_SHARED" & d$direction == "up", ]
    down <- d[d$pathway_id == "GO_SHARED" & d$direction == "down", ]
    ev_up <- dot_hover_event(cm, up)
    ev_down <- dot_hover_event(cm, down)
    expect_identical(c(ev_up$x, ev_up$y), c(ev_down$x, ev_down$y))

    session$setInputs(dot_hover = ev_down)
    html <- dot_tip(output$dot_tip)
    # The name in full, not wrapped as on the axis.
    expect_match(html, sprintf(">%s<", shared_name), fixed = TRUE)
    expect_match(html, "<td>Down</td>", fixed = TRUE)
    expect_match(html, "5 of 20 genes", fixed = TRUE)
    expect_match(html, "1.80e-248", fixed = TRUE)

    session$setInputs(dot_hover = ev_up)
    html <- dot_tip(output$dot_tip)
    expect_match(html, "<td>Up</td>", fixed = TRUE)
    expect_match(html, "10 of 40 genes", fixed = TRUE)
    expect_match(html, "1.00e-06", fixed = TRUE)

    # A tap shows it too; away from every dot, nothing.
    session$setInputs(dot_hover = NULL)
    session$setInputs(dot_click = ev_down)
    expect_match(dot_tip(output$dot_tip), "5 of 20 genes", fixed = TRUE)
    far <- ev_down
    far$coords_css$x <- far$coords_css$x - 60
    far$coords_img <- far$coords_css
    session$setInputs(dot_click = far)
    expect_no_match(dot_tip(output$dot_tip), "plot-hover-tip", fixed = TRUE)
  })
})

test_that("the card names the p-value the plot is coloured by", {
  b <- hover_ora_bundle()
  shiny::testServer(enrich_view_server, args = list(diff_bundle = shiny::reactiveVal(NULL)), {
    enrich_bundle(b); is_demo(FALSE)
    session$setInputs(show_p = "raw", show_cutoff = 0.05)
    session$flushReact()
    d <- results$dot_plot()$data
    row <- d[d$pathway_id == "GO_UP", ]
    session$setInputs(dot_hover = dot_hover_event(output$dot$coordmap, row))
    html <- dot_tip(output$dot_tip)
    expect_match(html, "<th>p</th>", fixed = TRUE)
    expect_match(html, "1.00e-04", fixed = TRUE)
    expect_no_match(html, "adjusted p", fixed = TRUE)
  })
})

test_that("a GSEA dot gives its NES and pathway size", {
  b <- hover_gsea_bundle()
  shiny::testServer(enrich_view_server, args = list(diff_bundle = shiny::reactiveVal(NULL)), {
    enrich_bundle(b); is_demo(FALSE)
    session$setInputs(show_p = "adjusted", show_cutoff = 0.05)
    session$flushReact()
    d <- results$dot_plot()$data
    row <- d[d$pathway_id == "GO_SHARED", ]
    session$setInputs(dot_hover = dot_hover_event(output$dot$coordmap, row))
    html <- dot_tip(output$dot_tip)
    expect_match(html, sprintf(">%s<", shared_name), fixed = TRUE)
    expect_match(html, "+2.35", fixed = TRUE)
    expect_match(html, "<th>NES</th>", fixed = TRUE)
    expect_match(html, "40 genes", fixed = TRUE)
    expect_no_match(html, "gene list", fixed = TRUE)
  })
})

test_that("the dot plot's card grows with long names and several panels", {
  short <- omicsCore::plot_enrichment(example_enrich_bundle(), view = "dot", top_n = 12L)
  # Never below the card's old height; the keys sit under the panel, so
  # even a short result may need a little more.
  expect_gte(enrich_dot_height(short), ENRICH_DOT_MIN_PX)
  expect_identical(enrich_dot_height(NULL), ENRICH_DOT_MIN_PX)
  # Three panels of long names need more room than one of short ones,
  # and a phone adds room for the legends under the plot.
  b <- hover_ora_bundle()
  df <- b$results$enrich_result_df
  more <- df[rep(seq_len(nrow(df)), 3), ]
  more$pathway_id <- paste0(more$pathway_id, rep(1:3, each = nrow(df)))
  more$pathway_name <- paste(more$pathway_name, rep(c("ONE", "TWO", "THREE"), each = nrow(df)))
  more$database <- rep(c("go_bp", "hallmark", "reactome"), each = nrow(df))
  b$results$enrich_result_df <- more
  long <- omicsCore::plot_enrichment(b, view = "dot", top_n = 12L)
  expect_gt(enrich_dot_height(long), enrich_dot_height(short))
  expect_gt(enrich_dot_height(long, narrow = TRUE), enrich_dot_height(long) - 100)
  # The card's plot takes the height renderPlot() gives it.
  ui <- as.character(enrich_dot_card(shiny::NS("enrich")))
  expect_match(ui, 'data-hover-id="enrich-dot_hover"', fixed = TRUE)
  expect_match(ui, "height:auto", fixed = TRUE)
})
