# The results half of the differential view: header, notices, stat
# cards, volcano, hit table and the full-table downloads.
#
# Split out of mod_diff_view.R. These are plain functions called from
# inside diff_view_server()'s moduleServer(), not modules, so every
# output keeps its id ("diff-volcano", "diff-download_csv").

diff_results_server <- function(input, output, session, navigate, active, shown_bundle,
                                diff_bundle, diff_error, comparisons, settings_changed,
                                marked, p_col, p_label, fdr_cut_d, fc_cut_d) {
  output$header <- shiny::renderUI({
    a <- active()
    b <- shown_bundle()
    source_note <- htmltools::tags$span(
      class = "muted",
      if (a$is_demo) "demo project (built-in)"
      else sprintf("layer = %s", a$tag)
    )

    # Before the first result there is nothing to summarise. Three em
    # dashes separated by middots said that in a way that read as
    # damage rather than as an empty state, which is how it was
    # reported.
    if (is.null(b)) {
      return(view_header(
        title    = "Differential",
        subtitle = htmltools::tagList(
          htmltools::tags$span(class = "muted", "no result yet"),
          htmltools::HTML(" &middot; "),
          source_note
        )
      ))
    }

    stale <- !is.null(diff_error()) && !identical(diff_error(), CANCELLED_MESSAGE)
    comparison <- if (length(b$params$comparison)) gsub("_vs_", " vs ", b$params$comparison) else
      sprintf("%s vs %s",
              b$params$case_group %||% input$case %||% "case",
              b$params$control_group %||% input$control %||% "control")
    n_cmp <- length(comparisons())
    view_header(
      title    = "Differential",
      actions  = if (is.function(navigate)) {
        shiny::actionButton(session$ns("go_enrich"), "Next: Enrichment \u2192",
                            class = "btn btn-ghost")
      },
      subtitle = htmltools::tagList(
        omics_display(b$input_info$omics_type),
        if (n_cmp > 1L) htmltools::tagList(
          htmltools::HTML(" &middot; "),
          sprintf("%d comparisons", n_cmp)),
        htmltools::HTML(" &middot; "),
        comparison,
        htmltools::HTML(" &middot; "),
        b$params$method,
        htmltools::HTML(" &middot; "),
        source_note,
        if (stale) htmltools::tagList(
          htmltools::HTML(" &middot; "),
          htmltools::tags$span(style = "color:var(--warn);font-weight:600",
                               "previous result \u2014 the latest run failed"))
      )
    )
  })

  output$notices <- shiny::renderUI({
    err <- diff_error()
    missing_engines <- diff_missing_engines()
    tagged <- htmltools::tagList()
    if (identical(err, CANCELLED_MESSAGE)) {
      tagged <- htmltools::tagAppendChild(tagged, notice(
        "Cancelled", "Press Run analysis to start again.", kind = "info"))
    } else if (!is.null(err)) {
      # The engine's own sentence when it is one: hiding "'control'
      # is not a group. Groups: Control, TreatA" behind "See the
      # technical details below" made users open a fold to read it.
      plain <- !looks_internal_error(err)
      tagged <- htmltools::tagAppendChild(
        tagged,
        notice(title  = if (is.null(diff_bundle())) "The differential analysis could not run"
                        else "The latest run failed; the result below is from the previous run",
               detail = if (plain) err else diff_error_hint(err),
               kind   = "error",
               technical = if (!plain) err)
      )
    }
    if (isTRUE(settings_changed()) && is.null(err)) {
      tagged <- htmltools::tagAppendChild(tagged, notice(
        "The settings have changed since this result was computed",
        "The plots and tables below still show the previous result. Press Re-run to update them.",
        kind = "warn"))
    }
    # What the engine said about the result: an ignored covariate, a
    # scale conversion, genes too low to test.
    warns <- diff_bundle()$warnings
    if (length(warns)) {
      tagged <- htmltools::tagAppendChild(tagged, notice(
        title = "Notes on this result",
        detail = htmltools::tags$ul(lapply(unique(warns), htmltools::tags$li)),
        kind = "info"))
    }
    if (is.null(diff_bundle()) && is.null(err)) {
      tagged <- htmltools::tagAppendChild(
        tagged,
        notice(title  = "No result yet",
               detail = paste("Check the layer, the control group and the groups to",
                              "compare in the Parameters card, then press Run analysis."),
               kind   = "info")
      )
    }
    if (length(missing_engines)) {
      tagged <- htmltools::tagAppendChild(
        tagged,
        notice(
          title  = "Some engines are unavailable",
          detail = sprintf(
            "Not installed: %s. Install with `omicsCore::install_optional()`.",
            paste(missing_engines, collapse = ", ")
          ),
          kind = "info"
        )
      )
    }
    tagged
  })

  output$stats <- shiny::renderUI({
    b <- shown_bundle()
    if (is.null(b)) return(NULL)
    # A weighted contrast has no single case group; its "up" is the
    # direction of the contrast as written.
    case_lbl <- b$params$case_group %||% "contrast"
    if (length(case_lbl) > 1L) case_lbl <- "case"
    df <- marked()
    sig <- df[df$is_significant, , drop = FALSE]
    up_n   <- sum(sig$effect > 0, na.rm = TRUE)
    down_n <- sum(sig$effect < 0, na.rm = TRUE)
    top    <- if (nrow(sig) > 0L) sig[which.max(abs(sig$effect)), ] else NULL
    top_value <- if (is.null(top)) "\u2014" else as.character(top$feature_symbol[1L])
    top_trend <- if (is.null(top)) "no features pass thresholds"
                 else sprintf("%s %+.2f \u00B7 %s %.2g", omicsCore::effect_label(b),
                              top$effect[1L], p_label(), top[[p_col()]][1L])
    htmltools::tags$div(
      class = "stat-grid",
      {
        # Features with no p-value were not tested (too many missing
        # values, or set aside as too low to test); counting them as
        # tested overstated the screen.
        n_tested <- sum(!is.na(df$p_value))
        stat_card(
          label = "Tested features",
          value = if (n_tested < nrow(df))
                    sprintf("%s / %s", format(n_tested, big.mark = ","), format(nrow(df), big.mark = ","))
                  else format(nrow(df), big.mark = ","),
          trend = if (n_tested < nrow(df))
                    sprintf("%s not testable (missing values or too few counts)",
                            format(nrow(df) - n_tested, big.mark = ","))
                  else sprintf("%s, %s", b$params$method, b$params$comparison %||% "\u2014"),
          mono  = TRUE
        )
      },
      stat_card(
        label  = sprintf("Up in %s", case_lbl),
        value  = up_n,
        trend  = sprintf("%s > %.2f \u00B7 %s < %.3g", omicsCore::effect_label(b),
                         fc_cut_d(), p_label(), fdr_cut_d()),
        accent = "up"
      ),
      stat_card(
        label  = sprintf("Down in %s", case_lbl),
        value  = down_n,
        trend  = sprintf("%s < -%.2f \u00B7 %s < %.3g", omicsCore::effect_label(b),
                         fc_cut_d(), p_label(), fdr_cut_d()),
        accent = "down"
      ),
      stat_card(
        label = "Top hit",
        value = top_value,
        trend = top_trend,
        mono  = TRUE
      )
    )
  })

  output$volcano <- plotly::renderPlotly({
    b <- shown_bundle()
    shiny::validate(shiny::need(b, "Press Run analysis to draw the volcano."))
    # Deliberately not given the slider values. The volcano is drawn
    # at plot_volcano()'s own defaults, which is what an exported
    # report and an exported script also produce -- so the figure a
    # reader is shown is the figure they can reproduce, and a
    # screenshot does not depend on where a control happened to be.
    #
    # The sliders still drive the hit table and the stat cards, where
    # sweeping a threshold is the useful thing to do; the figure is
    # the stable reference next to them.
    # The labels are added as plotly annotations, not drawn by
    # plot_volcano(): ggplotly() cannot convert ggrepel's text layer and
    # dropped it with a warning, so "Label top 20" labelled nothing.
    p <- omicsCore::plot_volcano(b, top_n = 0L)
    fig <- plotly::ggplotly(p, tooltip = "text") |> drop_hoveron()
    if (isTRUE(input$label_top)) {
      fig <- plotly::layout(fig, annotations = volcano_annotations(b, 20L))
    }
    # WebGL rather than one SVG node per point: 60,000 genes painted in
    # 0.5 s instead of 4 s, with the same points and hover text. Only
    # where the browser has WebGL: without it (remote desktops, some
    # locked-down or GPU-less machines) plotly drew a grey box saying
    # "WebGL is not supported" and no volcano at all.
    if (!identical(session$rootScope()$input$omics_webgl, FALSE)) {
      fig <- plotly::toWebGL(fig)
    }
    plotly::config(fig, displaylogo = FALSE,
                   modeBarButtonsToRemove = c("lasso2d", "select2d"))
  })

  # An empty card before the first run read as a broken one.
  output$hits_empty <- shiny::renderUI({
    if (!is.null(diff_bundle())) return(NULL)
    htmltools::tags$p(class = "muted", style = "font-size:13px;margin:4px 0",
                      "Run the analysis to see the features that pass the thresholds.")
  })

  output$hits <- DT::renderDT({
    df <- marked()
    sig <- df[df$is_significant, , drop = FALSE]
    sig <- sig[order(-abs(sig$effect)), , drop = FALSE]
    out <- data.frame(
      Feature   = sig$feature_symbol,
      Effect    = round(sig$effect, 3),
      p         = signif(sig[[p_col()]], 3),
      Direction = ifelse(sig$effect > 0, "up", "down"),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    # The columns are named for what they hold: the p-value the mask
    # was read from, and the effect in the words the cards use.
    names(out)[2] <- omicsCore::effect_label(shown_bundle())
    names(out)[3] <- p_label()
    DT::datatable(
      out,
      rownames  = FALSE,
      selection = "single",
      options   = list(
        pageLength = 10,
        dom        = "ftip",
        language   = list(emptyTable = "No feature passes the current thresholds."),
        scrollX    = TRUE,
        columnDefs = list(list(className = "dt-right", targets = 1:2))
      )
    )
  }, server = TRUE)
  invisible()
}

diff_downloads_server <- function(input, output, session, active, diff_bundle, p_col,
                                  fdr_cut_d, fc_cut_d) {
  # The whole result, every feature and every comparison, with the
  # significance at the current thresholds: the table a user takes to
  # a paper or a colleague. The hit list on screen is one page of one
  # comparison.
  full_table <- function() {
    b <- diff_bundle()
    shiny::req(b)
    df <- b$results$diff_result_df
    pv <- df[[p_col()]]
    df$significant <- !is.na(pv) & !is.na(df$effect) &
      pv < fdr_cut_d() & abs(df$effect) >= fc_cut_d()
    keep <- intersect(c("comparison", "feature_id", "feature_symbol", "effect",
                        "effect_type", "statistic", "p_value", "adj_p_value",
                        "base_mean", "significant", "method"), names(df))
    df[, keep, drop = FALSE]
  }
  table_name <- function(ext) {
    sprintf("differential_%s_%s.%s", active()$tag %||% "layer",
            format(Sys.Date(), "%Y%m%d"), ext)
  }
  output$download_csv <- shiny::downloadHandler(
    filename = function() table_name("csv"),
    content = function(file) utils::write.csv(full_table(), file, row.names = FALSE,
                                              fileEncoding = "UTF-8")
  )
  output$download_xlsx <- shiny::downloadHandler(
    filename = function() table_name("xlsx"),
    content = function(file) {
      df <- full_table()
      sheets <- split(df, df$comparison %||% "result")
      names(sheets) <- substr(gsub("[\\[\\]*?/:]", "_", names(sheets)), 1L, 31L)
      openxlsx::write.xlsx(sheets, file)
    }
  )
  invisible()
}

# Which Bioconductor diff backends aren't installed in this R
# session. Returned as a character vector for the notices strip.
# DESeq2 / edgeR / limma are the ones run_diff() can dispatch to.
diff_missing_engines <- function() {
  engines <- c(limma = "limma", DESeq2 = "DESeq2", edgeR = "edgeR")
  missing <- vapply(engines, function(pkg) !has_pkg(pkg), logical(1))
  names(engines)[missing]
}

diff_volcano_card <- function(ns) {
  bslib::card(
    bslib::card_header(
      htmltools::tags$h3(class = "card-title", "Volcano"),
      # Says plainly that the thresholds do not reach this figure. The
      # cut it was drawn at is in the plot's own caption, so a
      # screenshot carries it too.
      htmltools::tags$span(
        class = "card-sub",
        "fixed thresholds \u00B7 sliders filter the table below")
    ),
    bslib::card_body(
      plotly::plotlyOutput(ns("volcano"), height = "360px"),
      htmltools::tags$div(
        class = "legend",
        legend_swatch("significant", omics_colors$up),
        legend_swatch("ns", omics_colors$ns)
      )
    )
  )
}

diff_hits_card <- function(ns) {
  bslib::card(
    bslib::card_header(
      htmltools::tags$h3(class = "card-title", "Top hits"),
      htmltools::tags$span(class = "card-sub",
                           "largest changes first, within current thresholds")
    ),
    bslib::card_body(
      shiny::uiOutput(ns("hits_empty")),
      DT::DTOutput(ns("hits"), fill = FALSE),
      htmltools::tags$div(
        style = "display:flex;gap:8px;margin-top:8px;flex-wrap:wrap",
        shiny::downloadButton(ns("download_csv"), "Download all results (CSV)",
                              class = "btn btn-sm btn-ghost"),
        shiny::downloadButton(ns("download_xlsx"), "Excel, one sheet per comparison",
                              class = "btn btn-sm btn-ghost")
      )
    )
  )
}

# An R-internal error (a subscript, a missing object) rather than a
# sentence written for the user.
looks_internal_error <- function(msg) {
  grepl(paste0("Error in|subscript|object '.*' not found|non-numeric argument|",
               "argument is of length zero|missing value where|unused argument|",
               "could not find function|values must be length|replacement has|",
               "non-conformable|NA/NaN/Inf|invalid 'type'|undefined columns"),
        msg %||% "")
}

# A plain-language reading of the errors a run most often ends in; the
# message itself stays available under "Technical details".
diff_error_hint <- function(msg) {
  msg <- msg %||% ""
  if (grepl("not a level of", msg, fixed = TRUE)) {
    return("One of the chosen groups is not in this layer's group column.")
  }
  if (grepl("confounded", msg, fixed = TRUE)) {
    return("A covariate cannot be separated from the groups being compared; remove it.")
  }
  if (grepl("missing values", msg, fixed = TRUE)) {
    return("A covariate or pairing column has empty cells for some samples.")
  }
  if (grepl("paired design", msg, fixed = TRUE)) {
    return("The pairing column does not pair the samples one-to-one across the groups.")
  }
  if (grepl("required|not installed", msg)) {
    return("The chosen method needs a package that is not installed; pick another method.")
  }
  "See the technical details below."
}

# ggplotly() sets `hoveron` on its traces; scattergl has no such
# attribute, so once toWebGL() converts them plotly warns about it on
# every build. Dropped first, it is never there to warn about.
# The most significant features of a result as plotly annotations, at
# the coordinates plot_volcano() draws them (adjusted p when the result
# has it).
volcano_annotations <- function(bundle, n) {
  df <- bundle$results$diff_result_df
  p_col <- if ("adj_p_value" %in% names(df) && any(!is.na(df$adj_p_value))) "adj_p_value" else "p_value"
  df <- df[!is.na(df[[p_col]]) & !is.na(df$effect), , drop = FALSE]
  top <- utils::head(df[order(df[[p_col]]), , drop = FALSE], n)
  if (!nrow(top)) return(list())
  lab <- ifelse(is.na(top$feature_symbol) | !nzchar(top$feature_symbol),
                top$feature_id, top$feature_symbol)
  lapply(seq_len(nrow(top)), function(i) list(
    x = top$effect[[i]], y = -log10(max(top[[p_col]][[i]], .Machine$double.xmin)),
    text = lab[[i]], showarrow = TRUE, arrowhead = 0, arrowwidth = 0.8,
    arrowcolor = "#9AA3AE", ax = if (top$effect[[i]] >= 0) 24 else -24, ay = -14,
    font = list(size = 11, color = "#1A2541")))
}

drop_hoveron <- function(fig) {
  fig$x$data <- lapply(fig$x$data, function(tr) {
    tr$hoveron <- NULL
    tr
  })
  fig
}
