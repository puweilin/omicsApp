#' Differential analysis view module
#'
#' First **interactive** view in the app. Drives a volcano plot, a
#' top-hits table, and four stat cards from a cached
#' `example_diff_bundle()`; the FDR and |log2FC| sliders re-derive
#' a significance mask on every change. Method / Contrast /
#' Covariate controls render so the visual matches the mockup but
#' are deliberately inert in this slice — wiring real
#' `omicsCore::run_diff()` invocations is a later phase.
#'
#' Reference markup: `omicsApp/mockup/index.html:753-915`.
#'
#' @param id Module namespace id.
#' @keywords internal
#' @noRd
diff_view_ui <- function(id) {
  ns <- shiny::NS(id)

  htmltools::tagList(
    view_header(
      title    = "Differential",
      subtitle = "Proteomics \u00B7 G2 vs G1 \u00B7 limma \u00B7 age-adjusted"
    ),
    shiny::uiOutput(ns("stats")),
    htmltools::tags$div(
      class = "row-grid r-3-9",
      diff_params_card(ns),
      htmltools::tags$div(
        diff_volcano_card(ns),
        diff_hits_card(ns)
      )
    )
  )
}

#' @rdname diff_view_ui
#' @keywords internal
#' @noRd
diff_view_server <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {
    bundle <- shiny::reactive(example_diff_bundle())

    # Slider-derived significance mask. Recomputed on every change to
    # input$fdr_cut / input$fc_cut and shared by the stat cards,
    # volcano, and top-hits table.
    marked <- shiny::reactive({
      df <- bundle()$results$diff_result_df
      df$is_significant <- !is.na(df$adj_p_value) &
                           !is.na(df$effect) &
                           df$adj_p_value < input$fdr_cut &
                           abs(df$effect)  > input$fc_cut
      df
    })

    output$stats <- shiny::renderUI({
      df <- marked()
      sig <- df[df$is_significant, , drop = FALSE]
      up_n   <- sum(sig$effect > 0, na.rm = TRUE)
      down_n <- sum(sig$effect < 0, na.rm = TRUE)
      top    <- if (nrow(sig) > 0L) {
        sig[which.max(abs(sig$effect)), ]
      } else NULL
      top_value <- if (is.null(top)) "\u2014" else as.character(top$feature_symbol[1L])
      top_trend <- if (is.null(top)) "no features pass thresholds"
                   else sprintf("effect %+.2f \u00B7 adj.P %.2g",
                                top$effect[1L], top$adj_p_value[1L])

      htmltools::tags$div(
        class = "stat-grid",
        stat_card(
          label = "Tested features",
          value = format(nrow(df), big.mark = ","),
          trend = sprintf("%s, %s",
                          bundle()$params$method,
                          bundle()$params$comparison %||% "G2 vs G1"),
          mono  = TRUE
        ),
        stat_card(
          label  = "Up in G2",
          value  = up_n,
          trend  = sprintf("effect > %.2f \u00B7 adj.P < %.3f",
                           input$fc_cut, input$fdr_cut),
          accent = "up"
        ),
        stat_card(
          label  = "Down in G2",
          value  = down_n,
          trend  = sprintf("effect < -%.2f \u00B7 adj.P < %.3f",
                           input$fc_cut, input$fdr_cut),
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
      df <- marked()
      # Hand-roll the volcano so we can drive point colour from the
      # user mask and round-trip cleanly to plotly. plot_volcano()
      # would re-derive significance from the bundle, which is the
      # opposite of what we want here.
      df$.neglog10p <- -log10(pmax(df$adj_p_value, .Machine$double.xmin))
      df$.sig <- ifelse(df$is_significant, "significant", "ns")
      ranked <- df[order(df$adj_p_value, na.last = NA), , drop = FALSE]
      label_n <- if (isTRUE(input$label_top)) 20L else 0L
      top_ids <- if (label_n > 0L) utils::head(ranked$feature_id, label_n) else character(0)
      df$.label <- ifelse(df$feature_id %in% top_ids,
                          df$feature_symbol, NA_character_)

      p <- ggplot2::ggplot(
        df,
        ggplot2::aes(x = .data$effect, y = .data$.neglog10p,
                     color = .data$.sig,
                     text = sprintf(
                       "%s\neffect: %+.2f\nadj.P: %.2g",
                       .data$feature_symbol, .data$effect, .data$adj_p_value
                     ))
      ) +
        ggplot2::geom_point(alpha = 0.85, size = 1.8) +
        ggplot2::scale_color_manual(
          values = c(ns = "#9AA3AE", significant = "#C0392B"),
          name = NULL
        ) +
        ggplot2::geom_hline(yintercept = -log10(input$fdr_cut),
                            linetype = "dashed", color = "#9AA3AE") +
        ggplot2::geom_vline(xintercept = c(-input$fc_cut, input$fc_cut),
                            linetype = "dashed", color = "#9AA3AE") +
        ggplot2::labs(x = "log2 fold change", y = "-log10(adj.P)")

      if (any(!is.na(df$.label)) && requireNamespace("ggrepel", quietly = TRUE)) {
        p <- p + ggrepel::geom_text_repel(
          data = df,
          mapping = ggplot2::aes(x = .data$effect, y = .data$.neglog10p,
                                 label = .data$.label),
          size = 3, color = "#1A2541",
          max.overlaps = Inf, na.rm = TRUE, inherit.aes = FALSE
        )
      }

      plotly::ggplotly(p, tooltip = "text") |>
        plotly::config(displaylogo = FALSE,
                       modeBarButtonsToRemove = c("lasso2d", "select2d"))
    })

    output$hits <- DT::renderDT({
      df <- marked()
      sig <- df[df$is_significant, , drop = FALSE]
      sig <- sig[order(-abs(sig$effect)), , drop = FALSE]
      out <- data.frame(
        Feature   = sig$feature_symbol,
        Effect    = round(sig$effect, 3),
        `adj.P`   = signif(sig$adj_p_value, 3),
        Direction = ifelse(sig$effect > 0, "up", "down"),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      DT::datatable(
        out,
        rownames  = FALSE,
        selection = "single",
        options   = list(
          pageLength = 10,
          dom        = "ftip",
          scrollX    = TRUE,
          columnDefs = list(list(className = "dt-right", targets = 1:2))
        )
      )
    }, server = TRUE)
  })
}

# ---- internal helpers ------------------------------------------------

# Quiet R CMD check: the volcano ggplot uses `.data$col` selectors so
# the column names are resolved against the data frame, not the global
# environment. Listing `.data` here satisfies the static-binding check
# without pulling rlang into the call sites.
utils::globalVariables(".data")

# Convenience: %||% from rlang isn't worth a full import here.
`%||%` <- function(a, b) if (is.null(a)) b else a

diff_params_card <- function(ns) {
  bslib::card(
    bslib::card_header(
      htmltools::tags$h3(class = "card-title", "Parameters"),
      htmltools::tags$span(class = "card-sub",
                           "controls below the sliders are inert in this slice")
    ),
    bslib::card_body(
      htmltools::tags$div(
        class = "param-stack",
        param_group(
          "Method",
          shiny::selectInput(
            ns("method"), label = NULL,
            choices  = c("limma", "DESeq2", "edgeR", "ttest", "lm"),
            selected = "limma"
          )
        ),
        param_group(
          "Contrast",
          shiny::selectInput(
            ns("group_col"), label = "Group column",
            choices = c("group"), selected = "group"
          ),
          shiny::selectInput(
            ns("control"), label = "Control",
            choices = c("G1", "G2"), selected = "G1"
          ),
          shiny::selectInput(
            ns("case"), label = "Case",
            choices = c("G1", "G2"), selected = "G2"
          )
        ),
        param_group(
          "Covariates",
          htmltools::tags$div(
            class = "legend",
            pill("age", kind = "info"),
            htmltools::tags$span(class = "muted",
                                 style = "font-size:12px",
                                 "+ add\u2026")
          )
        ),
        param_group(
          "Thresholds",
          shiny::sliderInput(
            ns("fdr_cut"), label = "adj.P cutoff",
            min = 0, max = 0.2, value = 0.05, step = 0.005
          ),
          shiny::sliderInput(
            ns("fc_cut"), label = "|log2FC| cutoff",
            min = 0, max = 4, value = 1, step = 0.05
          ),
          shinyWidgets::materialSwitch(
            ns("label_top"), label = "Label top 20",
            value = FALSE, status = "primary", right = TRUE
          )
        ),
        htmltools::tags$div(
          style = "margin-top:8px",
          shiny::actionButton(
            ns("rerun"), "Re-run",
            class = "btn btn-primary",
            style = "width:100%"
          )
        )
      )
    )
  )
}

diff_volcano_card <- function(ns) {
  bslib::card(
    bslib::card_header(
      htmltools::tags$h3(class = "card-title", "Volcano"),
      htmltools::tags$span(class = "card-sub",
                           "drag to zoom \u00B7 double-click to reset")
    ),
    bslib::card_body(
      plotly::plotlyOutput(ns("volcano"), height = "360px"),
      htmltools::tags$div(
        class = "legend",
        legend_swatch("significant", "#C0392B"),
        legend_swatch("ns", "#9AA3AE"),
        htmltools::tags$span(class = "muted",
                             style = "font-size:12px",
                             "dashed lines mirror the slider thresholds")
      )
    )
  )
}

diff_hits_card <- function(ns) {
  bslib::card(
    bslib::card_header(
      htmltools::tags$h3(class = "card-title", "Top hits"),
      htmltools::tags$span(class = "card-sub",
                           "ranked by |effect| within current thresholds")
    ),
    bslib::card_body(
      DT::DTOutput(ns("hits"))
    )
  )
}
