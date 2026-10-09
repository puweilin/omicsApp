# The Selected feature card and the Heatmap card of the differential view.
#
# A feature is selected by clicking its row in Top hits or its point on
# the volcano (mod_diff_results.R); `selected` holds its feature id. The
# Selected feature card draws that feature's values in each group of the
# comparison on screen; the Heatmap card draws the comparison's top hits
# across its samples, the selected feature outlined.
#
# Called from inside diff_view_server()'s moduleServer(), like the other
# mod_diff_* files, so the outputs keep their ids ("diff-feature_plot",
# "diff-heatmap").

# The most rows the heatmap draws; past HEATMAP_NAMED_MAX they go unnamed
# (all but the selected one), where the names would overlap.
HEATMAP_MAX_ROWS <- 50L
HEATMAP_NAMED_MAX <- 40L

diff_detail_server <- function(input, output, session, active, shown_bundle, marked,
                               selected) {
  # The samples the comparison compares, and in which order to draw them.
  design <- shiny::reactive({
    b <- shown_bundle()
    shiny::req(omicsCore::is_analysis_bundle(b))
    diff_detail_design(b, active()$input)
  })

  # A feature the result on screen does not have (another layer's, after
  # the layer changed) is not selected any more. Another comparison of
  # the same layer has the same features, so the selection carries over
  # and the user can follow one gene from comparison to comparison.
  shiny::observeEvent(shown_bundle(), {
    id <- selected()
    b <- shown_bundle()
    if (is.null(id)) return()
    if (!omicsCore::is_analysis_bundle(b) || !id %in% b$results$diff_result_df$feature_id) {
      selected(NULL)
    }
  }, ignoreNULL = FALSE)
  shiny::observeEvent(input$clear_feature, selected(NULL))

  # The row of the result for the selected feature, with its
  # significance at the current thresholds.
  selected_row <- shiny::reactive({
    id <- selected()
    if (is.null(id)) return(NULL)
    df <- marked()
    i <- match(id, df$feature_id)
    if (is.na(i)) NULL else df[i, , drop = FALSE]
  })

  # The card shows its figure only while a feature is selected.
  has_feature <- shiny::reactive(!is.null(selected_row()))
  shiny::observe({
    show <- isTRUE(tryCatch(has_feature(), error = function(e) FALSE))
    tryCatch(shinyjs::toggle("feature_panel", condition = show), error = function(e) NULL)
  })

  output$feature_info <- shiny::renderUI({
    row <- selected_row()
    if (is.null(row)) {
      return(htmltools::tags$p(
        class = "muted", style = "font-size:13px;margin:4px 0",
        "Click a row in Top hits, or a point on the volcano, to see that feature's",
        "values in each group."))
    }
    b <- shown_bundle()
    name <- feature_display_name(row)
    verdict <- if (isTRUE(row$is_significant)) {
      if (row$effect > 0) "up, passes the current thresholds"
      else "down, passes the current thresholds"
    } else "does not pass the current thresholds"
    htmltools::tags$div(
      style = "display:flex;flex-wrap:wrap;align-items:baseline;gap:4px 14px;margin-bottom:4px",
      htmltools::tags$strong(style = "font-size:16px", name),
      if (!identical(name, row$feature_id))
        htmltools::tags$span(class = "muted", style = "font-size:12px;font-family:var(--font-mono)", row$feature_id),
      htmltools::tags$span(style = "font-size:13px;font-family:var(--font-mono)",
                           sprintf("%s %+.2f", omicsCore::effect_label(b), row$effect)),
      htmltools::tags$span(style = "font-size:13px;font-family:var(--font-mono)",
                           sprintf("adjusted p %s", format_p(row$adj_p_value))),
      htmltools::tags$span(
        style = sprintf("font-size:12px;color:%s",
                        if (!isTRUE(row$is_significant)) "var(--fg-3)"
                        else if (row$effect > 0) omics_colors$up else omics_colors$down),
        verdict),
      shiny::actionLink(session$ns("clear_feature"), "Clear", style = "font-size:12px")
    )
  })

  output$feature_note <- shiny::renderText({
    shiny::req(selected_row())
    d <- design()
    paste(c(
      if (length(d$levels)) sprintf("Each point is a sample; %s first.", d$levels[[1L]])
      else "Each point is a sample.",
      expression_scale_note(active()$input$assay_type)), collapse = " ")
  })

  # The figure, built here and drawn below: whatever saves the card's
  # figure takes this same ggplot.
  feature_plot <- shiny::reactive({
    row <- selected_row()
    shiny::req(row)
    inp <- active()$input
    d <- design()
    shiny::validate(
      shiny::need(inherits(inp, "omics_input") && !is.null(d$column),
                  "This layer's sample values are not available."),
      shiny::need(row$feature_id %in% rownames(inp$expr_mat),
                  "This feature is not in the layer's data table."))
    omicsCore::plot_feature_expression(inp, row$feature_id, group_by = d$column,
                                       group_levels = d$levels) +
      # The card says which feature it is, with its numbers.
      ggplot2::labs(title = NULL)
  })
  # Not through fit_to_width(): the figure has no legend to move, and the
  # card is narrow even on a desktop (beside the table), where the smaller
  # text it sets for a phone made the group names hard to read.
  output$feature_plot <- shiny::renderPlot(
    res = PLOT_RES, alt = "The selected feature's value in each sample, by group",
    feature_plot())

  # ---- heatmap -------------------------------------------------------

  # The hits at the current thresholds, most significant first, at most
  # HEATMAP_MAX_ROWS of them.
  heatmap_hits <- shiny::reactive({
    df <- marked()
    sig <- df[df$is_significant, , drop = FALSE]
    sig <- sig[order(sig$adj_p_value, sig$p_value), , drop = FALSE]
    inp <- active()$input
    ids <- sig$feature_id
    if (inherits(inp, "omics_input")) ids <- ids[ids %in% rownames(inp$expr_mat)]
    list(ids = utils::head(ids, HEATMAP_MAX_ROWS), n_sig = nrow(sig))
  })

  # Drawn only with two hits or more; the note says why otherwise.
  has_heatmap <- shiny::reactive({
    !is.null(shown_bundle()) && length(heatmap_hits()$ids) >= 2L
  })
  shiny::observe({
    show <- isTRUE(tryCatch(has_heatmap(), error = function(e) FALSE))
    tryCatch(shinyjs::toggle("heatmap_panel", condition = show), error = function(e) NULL)
  })

  output$heatmap_note <- shiny::renderUI({
    if (is.null(shown_bundle())) {
      return(htmltools::tags$p(class = "muted", style = "font-size:13px;margin:4px 0",
                               "Run the analysis to see the top hits across the samples."))
    }
    h <- heatmap_hits()
    if (length(h$ids) < 2L) {
      return(htmltools::tags$p(
        class = "muted", style = "font-size:13px;margin:4px 0",
        if (h$n_sig == 0L) "No feature passes the current thresholds, so there is no heatmap to draw."
        else "Only one feature passes the current thresholds; the heatmap needs at least two.",
        "Loosen the thresholds in the Parameters card to see more."))
    }
    n <- length(h$ids)
    htmltools::tags$p(
      class = "muted", style = "font-size:12px;margin:0 0 4px",
      if (h$n_sig > n) sprintf("The %d most significant of %s hits (by adjusted p; at most %d are drawn).",
                               n, format(h$n_sig, big.mark = ","), HEATMAP_MAX_ROWS)
      else sprintf("All %d hits, most significant first, clustered.", n),
      "Each row is scaled to its own mean and spread across the samples shown:",
      "red above the feature's mean, blue below.",
      if (n > HEATMAP_NAMED_MAX)
        sprintf("Rows are not named past %d; the selected feature still is.", HEATMAP_NAMED_MAX))
  })

  heatmap_plot <- shiny::reactive({
    b <- shown_bundle()
    h <- heatmap_hits()
    shiny::req(omicsCore::is_analysis_bundle(b), length(h$ids) >= 2L)
    inp <- active()$input
    d <- design()
    shiny::validate(shiny::need(inherits(inp, "omics_input"),
                                "This layer's sample values are not available."))
    omicsCore::plot_heatmap(
      b, input = inp, features = h$ids,
      group_by = d$column, group_levels = d$levels,
      highlight = selected(), engine = "ggplot", title = "",
      show_rownames = length(h$ids) <= HEATMAP_NAMED_MAX) +
      ggplot2::labs(title = NULL)
  })
  output$heatmap <- shiny::renderPlot(
    res = PLOT_RES, alt = "Heatmap of the top hits across the samples, grouped",
    height = function() {
      w <- session$clientData[[paste0("output_", session$ns("heatmap"), "_width")]]
      heatmap_height(length(heatmap_hits()$ids),
                     narrow = is.numeric(w) && length(w) && w < NARROW_PLOT_PX)
    },
    fit_to_width("heatmap", heatmap_plot()))

  # The two figures' ggplots, for whatever saves a card's figure, and
  # what the cards' panels are shown by.
  list(feature_plot = feature_plot, heatmap_plot = heatmap_plot,
       heatmap_hits = heatmap_hits, has_feature = has_feature,
       has_heatmap = has_heatmap)
}

# The card height for a heatmap of `n` rows: a row of about 13 px, the
# group bar and sample names around them, and room under the tiles for
# the legends where a phone puts them.
heatmap_height <- function(n, narrow = FALSE) {
  max(260, 110 + 13 * n) + if (isTRUE(narrow)) 150 else 0
}

# The meta column and groups the comparison on screen compares, the
# reference (control) group first; for a continuous design the variable
# and no groups. `column` is NULL when the layer has no such column.
diff_detail_design <- function(b, inp) {
  params <- b$params
  meta <- if (inherits(inp, "omics_input")) inp$meta_df else NULL
  if (identical(params$analysis_type, "continuous")) {
    col <- params$continuous_col
    return(list(column = if (length(col) == 1L && col %in% colnames(meta)) col,
                levels = NULL))
  }
  col <- params$group_col
  if (length(col) != 1L || is.null(meta) || !col %in% colnames(meta)) {
    return(list(column = NULL, levels = NULL))
  }
  present <- sort(unique(stats::na.omit(as.character(meta[[col]]))))
  control <- params$control_group
  case <- params$case_group
  lv <- if (length(control) == 1L && length(case) == 1L) {
    c(control, case)
  } else {
    # A weighted contrast ("(TreatA + TreatB)/2 - Control") has no single
    # case group: the groups its expression names, the reference first.
    spec <- paste(params$contrasts %||% "", collapse = " ")
    named <- present[vapply(present, function(g) grepl(g, spec, fixed = TRUE), logical(1))]
    c(control, setdiff(if (length(named)) named else present, control))
  }
  lv <- lv[lv %in% present]
  list(column = col, levels = if (length(lv)) lv else present)
}

# A feature's name as the table shows it: its symbol, or its id.
feature_display_name <- function(row) {
  sym <- row$feature_symbol %||% NA_character_
  if (length(sym) != 1L || is.na(sym) || !nzchar(sym)) row$feature_id else sym
}

format_p <- function(p) {
  if (length(p) != 1L || is.na(p)) return("—")
  formatC(p, digits = 2, format = if (p < 1e-3) "e" else "g")
}

# What the figure's values are, where the axis title alone does not say
# how they were made.
expression_scale_note <- function(assay_type) {
  switch(as.character(assay_type %||% ""),
         raw_count = "Counts are shown as log2(CPM + 1): scaled to each sample's library size, then logged.",
         raw_intensity = "Intensities are shown as log2(intensity + 1).",
         tpm = "TPM values are shown as log2(TPM + 1).",
         fpkm = "FPKM values are shown as log2(FPKM + 1).",
         NULL)
}

diff_feature_card <- function(ns) {
  bslib::card(
    bslib::card_header(
      htmltools::tags$h3(class = "card-title", "Selected feature"),
      htmltools::tags$span(class = "card-sub", "its value in each group of this comparison"),
      plot_download_ui(ns("feature_plot_download"))
    ),
    bslib::card_body(
      shiny::uiOutput(ns("feature_info")),
      shinyjs::hidden(htmltools::tags$div(
        id = ns("feature_panel"),
        shiny::plotOutput(ns("feature_plot"), height = "320px"),
        shiny::textOutput(ns("feature_note"), container = function(...)
          htmltools::tags$p(class = "muted", style = "font-size:12px;margin:4px 0 0", ...))
      ))
    )
  )
}

diff_heatmap_card <- function(ns) {
  bslib::card(
    bslib::card_header(
      htmltools::tags$h3(class = "card-title", "Heatmap"),
      htmltools::tags$span(class = "card-sub", "top hits across the samples, by group"),
      plot_download_ui(ns("heatmap_download"))
    ),
    bslib::card_body(
      shiny::uiOutput(ns("heatmap_note")),
      shinyjs::hidden(htmltools::tags$div(
        id = ns("heatmap_panel"),
        shiny::plotOutput(ns("heatmap"), height = "auto")
      ))
    )
  )
}
