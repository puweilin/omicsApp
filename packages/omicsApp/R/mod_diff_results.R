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

  # plotly drops a ggplot caption, so the card says what "significant"
  # means on the figure: the same cut as the table beside it.
  output$volcano_cut <- shiny::renderText({
    b <- shown_bundle()
    if (is.null(b)) return("")
    volcano_cut_text(p_label(), fdr_cut_d(), omicsCore::effect_label(b), fc_cut_d())
  })

  # The figure's legend, with the counts plot_volcano() puts in its own
  # (the same mask as the stat cards, so the numbers match them too).
  output$volcano_legend <- shiny::renderUI({
    if (is.null(shown_bundle())) return(volcano_legend())
    df <- marked()
    volcano_legend(
      up_n   = sum(df$is_significant & df$effect > 0, na.rm = TRUE),
      down_n = sum(df$is_significant & df$effect < 0, na.rm = TRUE),
      continuous = startsWith(as.character(df$analysis_type[1L] %||% ""), "continuous"))
  })

  output$volcano <- plotly::renderPlotly({
    b <- shown_bundle()
    shiny::validate(shiny::need(b, "Press Run analysis to draw the volcano."))
    # Drawn at the thresholds the controls set, so its colours agree
    # with the hit table and the stat cards beside it. The same values
    # are saved with the project (params$display_thresholds), so the
    # report and the exported script draw this same figure.
    # The labels are added as plotly annotations, not drawn by
    # plot_volcano(): ggplotly() cannot convert ggrepel's text layer and
    # dropped it with a warning, so "Label top 20" labelled nothing.
    p <- omicsCore::plot_volcano(b, top_n = 0L, p_basis = volcano_p_basis(p_col()),
                                 p_threshold = fdr_cut_d(),
                                 effect_threshold = volcano_effect_cut(fc_cut_d()))
    # plotly's own legend is hidden: the card's legend below the plot
    # says the same, with the counts, and does not take width from the
    # plot on a phone.
    fig <- plotly::ggplotly(p, tooltip = "text") |> drop_hoveron() |>
      plotly::layout(showlegend = FALSE)
    if (isTRUE(input$label_top)) {
      fig <- plotly::layout(fig, annotations = volcano_annotations(b, 20L, p_col()))
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
      # The cut the figure is drawn at, in words: plotly does not show
      # the caption plot_volcano() writes it into.
      shiny::textOutput(ns("volcano_cut"), container = function(...)
        htmltools::tags$span(class = "card-sub", ...))
    ),
    bslib::card_body(
      plotly::plotlyOutput(ns("volcano"), height = "360px"),
      shiny::uiOutput(ns("volcano_legend"))
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

# The volcano's significance cut, from the view's controls, in the form
# plot_volcano() takes it. A |log2FC| cutoff of 0 is no cutoff: drawn,
# it was a pair of dashed lines on top of the zero line.
volcano_p_basis <- function(p_col) {
  if (identical(p_col, "p_value")) "raw" else "adjusted"
}
volcano_effect_cut <- function(fc_cut) {
  if (is.null(fc_cut) || !is.finite(fc_cut) || fc_cut <= 0) NULL else fc_cut
}

# "significant = adjusted p < 0.05 and |log2FC| >= 0.263", in the words
# the controls use.
volcano_cut_text <- function(p_label, p_cut, effect_name, fc_cut) {
  txt <- sprintf("significant = %s < %s", p_label, format(p_cut, digits = 3))
  if (!is.null(volcano_effect_cut(fc_cut))) {
    txt <- sprintf("%s and |%s| \u2265 %s", txt, effect_name, format(fc_cut, digits = 3))
  }
  txt
}

# The volcano's legend, in plot_volcano()'s classes and colours: up in
# red, down in blue, the rest grey. Counts when there is a result.
volcano_legend <- function(up_n = NULL, down_n = NULL, continuous = FALSE) {
  words <- if (isTRUE(continuous)) c("positive", "negative") else c("up", "down")
  count <- function(word, n) {
    if (is.null(n)) word else sprintf("%s (%s)", word, format(n, big.mark = ","))
  }
  htmltools::tags$div(
    class = "legend",
    legend_swatch(count(words[[1]], up_n), omics_colors$up),
    legend_swatch(count(words[[2]], down_n), omics_colors$down),
    legend_swatch("not significant", omics_colors$ns)
  )
}

# The most significant features of a result as plotly annotations, at
# the coordinates plot_volcano() draws them: the p-value column the
# figure was drawn from.
volcano_annotations <- function(bundle, n, p_col = "adj_p_value") {
  df <- bundle$results$diff_result_df
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

# ggplotly() sets `hoveron` on its traces; scattergl has no such
# attribute, so once toWebGL() converts them plotly warns about it on
# every build. Dropped first, it is never there to warn about.
drop_hoveron <- function(fig) {
  fig$x$data <- lapply(fig$x$data, function(tr) {
    tr$hoveron <- NULL
    tr
  })
  fig
}
