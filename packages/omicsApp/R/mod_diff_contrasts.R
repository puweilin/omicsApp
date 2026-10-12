# The Comparisons card of the differential view: hits per comparison and
# their overlap, shown when one run fitted several contrasts.
#
# Split out of mod_diff_view.R and called from inside its module server,
# so the outputs keep their ids ("diff-contrast_plot"). Returns the
# per-comparison summary the table is drawn from, and the two figures.
diff_contrasts_server <- function(input, output, session, diff_bundle, comparisons,
                                  fdr_cut_d, fc_cut_d) {
  output$contrast_summary <- shiny::renderUI({
    if (length(comparisons()) < 2L) return(NULL)
    bslib::card(
      bslib::card_header(
        htmltools::tags$h3(class = "card-title", "Comparisons"),
        htmltools::tags$span(class = "card-sub",
                             "hits per comparison at the current thresholds"),
        plot_download_ui(session$ns("contrast_plot_download"))
      ),
      bslib::card_body(
        htmltools::tags$div(
          class = "row-grid r-6-6",
          shiny::plotOutput(session$ns("contrast_plot"),
                            height = paste0(160 + label_rows_px(comparisons(), 20L), "px")),
          DT::DTOutput(session$ns("contrast_table"))
        ),
        # Which comparisons share their hits. The table's "also in
        # another" says how many; this says with which.
        # The overlap figure's own download sits on its control row: the
        # card header's saves the bars above.
        htmltools::tags$div(
          class = "inline-control plot-download-row",
          shiny::radioButtons(session$ns("overlap_dir"), label = "Overlap of",
                              choices = c("all hits" = "any", "up" = "up",
                                          "down" = "down"),
                              selected = "any", inline = TRUE),
          plot_download_ui(session$ns("overlap_plot_download"))
        ),
        # Its height is the server's (overlap_plot_px()): a phone's
        # labels are shorter.
        shiny::plotOutput(session$ns("overlap_plot"), height = "auto")
      )
    )
  })

  # The two figures, built here and drawn below: the card's downloads
  # save these same ggplots.
  overlap_figure <- function(compact) {
    b <- diff_bundle()
    shiny::req(omicsCore::is_analysis_bundle(b), length(comparisons()) > 1L)
    omicsCore::plot_diff_overlap(
      b, p_cutoff = fdr_cut_d(),
      p_preference = if (identical(input$p_kind %||% "adj", "raw")) "raw" else "adjusted",
      effect_cutoff = fc_cut_d(),
      direction = input$overlap_dir %||% "any",
      compact = compact)
  }
  # The full figure is the one a download saves, whatever the screen.
  overlap_plot <- shiny::reactive(overlap_figure(compact = FALSE))
  # On a phone, five comparisons' names took most of the width and the
  # counts over twenty narrow bars ran into each other: there the card
  # draws the compact layout -- the six largest combinations, treatment
  # names alone when the comparisons share a control.
  overlap_narrow <- shiny::reactive(is_narrow_width(plot_output_width("overlap_plot", session)))
  overlap_shown <- shiny::reactive({
    if (overlap_narrow()) overlap_figure(compact = TRUE) else overlap_plot()
  })
  output$overlap_plot <- shiny::renderPlot(
    res = PLOT_RES, alt = "Overlap of the significant features between comparisons",
    height = function() overlap_plot_px(comparisons(), overlap_narrow()),
    fit_to_width("overlap_plot", overlap_shown()))

  contrast_summary_df <- shiny::reactive({
    b <- diff_bundle()
    shiny::req(omicsCore::is_analysis_bundle(b), length(comparisons()) > 1L)
    omicsCore::summarize_diff_contrasts(
      b, p_cutoff = fdr_cut_d(),
      p_preference = if (identical(input$p_kind %||% "adj", "raw")) "raw" else "adjusted",
      effect_cutoff = fc_cut_d())
  })

  contrast_plot <- shiny::reactive({
    b <- diff_bundle()
    shiny::req(omicsCore::is_analysis_bundle(b), length(comparisons()) > 1L)
    omicsCore::plot_diff_contrasts(
      b, p_cutoff = fdr_cut_d(),
      p_preference = if (identical(input$p_kind %||% "adj", "raw")) "raw" else "adjusted",
      effect_cutoff = fc_cut_d())
  })
  output$contrast_plot <- shiny::renderPlot(
    res = PLOT_RES, alt = "Number of up- and down-regulated features per comparison",
    fit_to_width("contrast_plot", contrast_plot()))

  output$contrast_table <- DT::renderDT({
    s <- contrast_summary_df()
    out <- data.frame(
      Comparison = gsub("_vs_", " vs ", s$comparison),
      Up = s$n_up, Down = s$n_down,
      `Also in another` = s$n_shared,
      check.names = FALSE, stringsAsFactors = FALSE)
    DT::datatable(out, rownames = FALSE, selection = "none",
                  options = list(dom = "t", pageLength = 50))
  }, server = TRUE)

  list(summary = contrast_summary_df, contrast_plot = contrast_plot,
       overlap_plot = overlap_plot)
}

# The overlap plot's height in CSS px: the bars, then a row per
# comparison as tall as its label. A phone's labels are the compact
# figure's (comparison_label_lines(compact = TRUE)), in a phone's smaller
# text (a line about 15 px), and its subtitle has a line or two more
# ("top 6 of 14 combinations").
overlap_plot_px <- function(comparisons, narrow = FALSE) {
  if (!isTRUE(narrow)) return(270L + label_rows_px(comparisons, 10L))
  lines <- omicsCore::comparison_label_lines(comparisons, compact = TRUE)
  as.integer(290L + sum(10L + 15L * lines))
}
