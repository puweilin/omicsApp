#' Pathway enrichment view module
#'
#' Slice 2E: pathway dotplot (left) and enriched-sets DT (right) from
#' a synthetic GSEA-style fixture (`example_enrich_table()`). The
#' inline ORA / GSEA / GSVA tab strip mirrors the mockup but is inert
#' in this slice — clicking does nothing. Live
#' `omicsCore::run_enrichment()` invocations are deferred until the
#' Bioconductor-suggest gating story is settled.
#'
#' Reference markup: `omicsApp/mockup/index.html:918-964`.
#'
#' @param id Module namespace id.
#' @keywords internal
#' @noRd
enrich_view_ui <- function(id) {
  ns <- shiny::NS(id)

  htmltools::tagList(
    view_header(
      title    = "Pathway enrichment",
      subtitle = "Proteomics \u00B7 G2 vs G1 \u00B7 MSigDB Hallmark \u00B7 15 sets shown"
    ),
    htmltools::tags$div(
      class = "inline-tabs",
      htmltools::tags$button(type = "button",
                             class = "inline-tab active", "ORA"),
      htmltools::tags$button(type = "button",
                             class = "inline-tab", "GSEA"),
      htmltools::tags$button(type = "button",
                             class = "inline-tab", "GSVA")
    ),
    htmltools::tags$div(
      class = "row-grid r-7-5",
      enrich_dot_card(ns),
      enrich_hits_card(ns)
    )
  )
}

#' @rdname enrich_view_ui
#' @keywords internal
#' @noRd
enrich_view_server <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {
    table_data <- shiny::reactive(example_enrich_table())

    output$dot <- shiny::renderPlot({
      df <- table_data()
      top <- df[order(df$adj_p_value), , drop = FALSE]
      top <- top[seq_len(min(12L, nrow(top))), , drop = FALSE]
      top$pathway_name <- factor(top$pathway_name,
                                 levels = rev(top$pathway_name))
      ggplot2::ggplot(
        top,
        ggplot2::aes(x = .data$effect,
                     y = .data$pathway_name,
                     size = .data$overlap_size,
                     color = -log10(.data$adj_p_value))
      ) +
        ggplot2::geom_point() +
        ggplot2::scale_color_gradient(low = "#9AA3AE", high = "#C0392B") +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                            color = "#9AA3AE") +
        ggplot2::labs(x = "NES", y = NULL,
                      size = "Overlap",
                      color = "-log10(adj.P)") +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
          panel.grid.minor = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_text(color = "#1A2541")
        )
    })

    output$hits <- DT::renderDT({
      df <- table_data()
      df <- df[order(df$adj_p_value), , drop = FALSE]
      out <- data.frame(
        Pathway   = df$pathway_name,
        NES       = sprintf("%+.2f", df$effect),
        `adj.P`   = signif(df$adj_p_value, 3),
        Direction = df$direction,
        Overlap   = sprintf("%d/%d", df$overlap_size, df$gene_set_size),
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
          columnDefs = list(list(className = "dt-right",
                                 targets = c(1, 2)))
        )
      )
    }, server = TRUE)
  })
}

# ---- internal helpers ------------------------------------------------

# Same .data trick as mod_diff_view: silence R CMD check for the
# tidy-eval column selectors in the dotplot ggplot.
utils::globalVariables(".data")

enrich_dot_card <- function(ns) {
  bslib::card(
    bslib::card_header(
      htmltools::tags$h3(class = "card-title", "Pathway dotplot"),
      htmltools::tags$span(class = "card-sub",
                           "top 12 by adj.P \u00B7 size = overlap")
    ),
    bslib::card_body(
      shiny::plotOutput(ns("dot"), height = "420px")
    )
  )
}

enrich_hits_card <- function(ns) {
  bslib::card(
    bslib::card_header(
      htmltools::tags$h3(class = "card-title", "Enriched sets"),
      htmltools::tags$span(class = "card-sub",
                           "ranked by adj.P")
    ),
    bslib::card_body(
      DT::DTOutput(ns("hits"))
    )
  )
}
