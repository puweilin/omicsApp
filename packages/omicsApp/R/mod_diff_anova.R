# The global-test card of the differential view.
#
# "Does this feature differ between any of the groups?" -- one test per
# feature over every group at once, before (or instead of) reading the
# comparisons one by one. Kept beside the comparisons rather than
# replacing them: it has no direction and no fold change, so nothing
# downstream (enrichment, integration) can use it.
#
# Split out of mod_diff_view.R and called from inside its module server,
# so the ids stay "diff-run_anova", "diff-anova_table". The result, its
# error and its busy flag are the module's (a change of layer resets
# them); this returns the sorted hit table.
diff_anova_server <- function(input, output, session, active, default_contrast, levels_,
                              anova_bundle, anova_error, anova_running, anova_epoch,
                              p_col, p_label, fdr_cut_d) {

  output$anova_card <- shiny::renderUI({
    if (length(levels_()) < 3L) return(NULL)
    b <- anova_bundle()
    bslib::card(
      bslib::card_header(
        htmltools::tags$h3(class = "card-title", "Any difference between groups"),
        htmltools::tags$span(class = "card-sub",
                             "one global test per feature (ANOVA / LRT)"),
        info_tip(paste("Tests all groups of the column at once: a small p says the",
                       "feature differs somewhere among them, not where. Use it to",
                       "screen, then read the comparisons for the direction."))
      ),
      bslib::card_body(
        htmltools::tags$div(
          style = "display:flex;gap:12px;align-items:center;flex-wrap:wrap",
          disabled_if(shiny::actionButton(session$ns("run_anova"),
                                          if (is.null(b)) "Run global test" else "Re-run global test",
                                          class = "btn btn-sm btn-outline-primary"),
                      shiny::isolate(anova_running())),
          shiny::uiOutput(session$ns("anova_summary"), inline = TRUE)
        ),
        if (!is.null(anova_error())) {
          notice("The global test could not run", kind = "error",
                 technical = anova_error())
        },
        if (!is.null(b)) DT::DTOutput(session$ns("anova_table"), fill = FALSE),
        if (!is.null(b)) shiny::downloadButton(session$ns("download_anova"),
                                               "Download all (CSV)",
                                               class = "btn btn-sm btn-ghost")
      )
    )
  })

  shiny::observeEvent(input$run_anova, {
    a <- active()
    shiny::req(a$input)
    d <- default_contrast()
    group_col <- input$group_col %||% d$group_col
    method <- input$method %||% "auto"
    if (!method %in% c("limma", "edger", "deseq2")) method <- "auto"
    covariates <- input$covariates
    my_run <- anova_epoch$start()
    set_button_busy("run_anova", TRUE, anova_running)
    run_async(
      detached_call(
        function() {
          omicsCore::run_diff(input = inp, method = method,
                              analysis_type = "anova", group_col = group_col,
                              covariates = covariates)
        },
        inp = a$input, method = method, group_col = group_col,
        covariates = if (length(covariates)) covariates else NULL
      ),
      on_success = function(bundle) {
        if (anova_epoch$is_last_started(my_run)) set_button_busy("run_anova", FALSE, anova_running)
        if (!anova_epoch$is_current(my_run)) return(invisible())
        anova_error(NULL)
        anova_bundle(with_layer(bundle, a))
      },
      on_error = function(msg) {
        if (anova_epoch$is_last_started(my_run)) set_button_busy("run_anova", FALSE, anova_running)
        if (anova_epoch$is_current(my_run)) anova_error(msg)
      },
      message = "Running the global test..."
    )
  })

  anova_hits <- shiny::reactive({
    b <- anova_bundle()
    shiny::req(b)
    df <- b$results$diff_result_df
    df[order(df[[p_col()]], na.last = TRUE), , drop = FALSE]
  })

  output$anova_summary <- shiny::renderUI({
    df <- anova_hits()
    n <- sum(df[[p_col()]] < fdr_cut_d(), na.rm = TRUE)
    htmltools::tags$span(
      htmltools::tags$strong(format(n, big.mark = ",")),
      sprintf(" of %s features differ between the groups of '%s' (%s < %.3f)",
              format(nrow(df), big.mark = ","), anova_bundle()$params$group_col,
              p_label(), fdr_cut_d()))
  })

  output$anova_table <- DT::renderDT({
    df <- anova_hits()
    out <- data.frame(
      Feature = df$feature_symbol %||% df$feature_id,
      Statistic = signif(df$statistic, 3),
      p = signif(df[[p_col()]], 3),
      check.names = FALSE, stringsAsFactors = FALSE)
    names(out)[2] <- df$statistic_type[1] %||% "Statistic"
    names(out)[3] <- p_label()
    DT::datatable(out, rownames = FALSE, selection = "none",
                  options = list(pageLength = 10, dom = "ftip"))
  }, server = TRUE)

  output$download_anova <- shiny::downloadHandler(
    filename = function() sprintf("global_test_%s.csv", active()$tag %||% "layer"),
    content = function(file) {
      b <- anova_bundle()
      shiny::req(b)
      utils::write.csv(b$results$diff_result_df, file, row.names = FALSE, fileEncoding = "UTF-8")
    }
  )

  anova_hits
}
