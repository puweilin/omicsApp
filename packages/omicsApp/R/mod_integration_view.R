#' Multi-omics integration view module
#'
#' Slice 2E: four stat cards across the top, a mirrored-volcano /
#' log2FC-scatter row, and a shared-pathway dotplot / ActivePathways
#' DT row. Driven by a synthetic fixture
#' (`example_integration_tables()`) rather than live
#' `omicsCore::run_integration()` so the demo doesn't depend on the
#' `ActivePathways` Bioconductor Suggests. Mapping / Export buttons
#' render but are inert in this slice.
#'
#' Reference markup: `omicsApp/mockup/index.html:967-1055`.
#'
#' @param id Module namespace id.
#' @keywords internal
#' @noRd
integration_view_ui <- function(id) {
  ns <- shiny::NS(id)

  fixt <- example_integration_tables()
  conc <- fixt$concordance_df

  paired_n <- nrow(conc)
  up_n     <- sum(conc$quadrant == "up_up",   na.rm = TRUE)
  down_n   <- sum(conc$quadrant == "down_down", na.rm = TRUE)
  # Spearman over the seeded signal subset so the headline number
  # reflects "concordance among significant features", not noise.
  sig <- conc[conc$quadrant %in% c("up_up", "down_down",
                                   "up_down", "down_up"), , drop = FALSE]
  rho <- if (nrow(sig) >= 3L) {
    suppressWarnings(stats::cor(sig$effect_a, sig$effect_b,
                                method = "spearman",
                                use = "pairwise.complete.obs"))
  } else NA_real_

  htmltools::tagList(
    view_header(
      title    = "Multi-omics integration",
      subtitle = "Proteomics \u00D7 RNA-seq \u00B7 G2 vs G1 \u00B7 60 paired features",
      actions  = htmltools::tagList(
        htmltools::tags$button(type = "button",
                               class = "btn btn-ghost",
                               "Mapping"),
        htmltools::tags$button(type = "button",
                               class = "btn btn-ghost",
                               "Export integration bundle")
      )
    ),
    htmltools::tags$div(
      class = "stat-grid",
      stat_card(
        label  = "Paired features",
        value  = format(paired_n, big.mark = ","),
        trend  = "RNA \u2194 protein mapping",
        accent = "brand", mono = TRUE
      ),
      stat_card(
        label  = "Concordant \u2191",
        value  = up_n,
        trend  = "both layers up",
        accent = "up", mono = TRUE
      ),
      stat_card(
        label  = "Concordant \u2193",
        value  = down_n,
        trend  = "both layers down",
        accent = "down", mono = TRUE
      ),
      stat_card(
        label  = "Spearman \u03C1 (sig)",
        value  = if (is.na(rho)) "\u2014" else sprintf("%.2f", rho),
        trend  = "log2FC RNA vs Protein",
        mono   = TRUE
      )
    ),
    htmltools::tags$div(
      class = "row-grid r-6-6",
      integration_dual_card(ns),
      integration_scatter_card(ns)
    ),
    htmltools::tags$div(
      class = "row-grid r-6-6",
      integration_shared_card(ns),
      integration_ap_card(ns)
    )
  )
}

#' @rdname integration_view_ui
#' @keywords internal
#' @noRd
integration_view_server <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {
    fixt <- shiny::reactive(example_integration_tables())

    output$dual <- shiny::renderPlot({
      df <- fixt()$concordance_df
      df$.neglog10p <- -log10(pmax(df$adj_p_value, .Machine$double.xmin))
      ggplot2::ggplot(
        df,
        ggplot2::aes(x = .data$effect_diff,
                     y = .data$.neglog10p,
                     color = .data$quadrant)
      ) +
        ggplot2::geom_point(alpha = 0.85, size = 2) +
        ggplot2::scale_color_manual(values = integration_quadrant_colors()) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                            color = "#9AA3AE") +
        ggplot2::labs(x = "log2FC(protein) - log2FC(RNA)",
                      y = "-log10(adj.P)",
                      color = NULL) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
    })

    output$scatter <- shiny::renderPlot({
      df <- fixt()$concordance_df
      ggplot2::ggplot(
        df,
        ggplot2::aes(x = .data$effect_a,
                     y = .data$effect_b,
                     color = .data$quadrant)
      ) +
        ggplot2::geom_abline(slope = 1, intercept = 0,
                             linetype = "dashed", color = "#9AA3AE") +
        ggplot2::geom_hline(yintercept = 0, color = "#E5E7EB") +
        ggplot2::geom_vline(xintercept = 0, color = "#E5E7EB") +
        ggplot2::geom_point(alpha = 0.85, size = 2) +
        ggplot2::scale_color_manual(values = integration_quadrant_colors()) +
        ggplot2::labs(x = "log2FC RNA", y = "log2FC protein",
                      color = NULL) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
    })

    output$shared <- shiny::renderPlot({
      ap <- fixt()$active_pathways_df
      long <- data.frame(
        pathway_name = rep(ap$pathway_name, times = 2L),
        omics        = rep(c("protein", "RNA"), each = nrow(ap)),
        neglog10p    = c(-log10(ap$p_a), -log10(ap$p_b)),
        direction    = rep(ap$direction, times = 2L),
        stringsAsFactors = FALSE
      )
      long$pathway_name <- factor(long$pathway_name,
                                  levels = rev(ap$pathway_name))
      long$omics <- factor(long$omics, levels = c("protein", "RNA"))
      ggplot2::ggplot(
        long,
        ggplot2::aes(x = .data$omics,
                     y = .data$pathway_name,
                     size = .data$neglog10p,
                     color = .data$direction)
      ) +
        ggplot2::geom_point() +
        ggplot2::scale_color_manual(values = c(shared = "#1F4E96",
                                               unique = "#9AA3AE")) +
        ggplot2::scale_size_continuous(range = c(3, 10)) +
        ggplot2::labs(x = NULL, y = NULL,
                      size = "-log10(p)", color = NULL) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
          panel.grid.minor = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_text(color = "#1A2541")
        )
    })

    output$ap_table <- DT::renderDT({
      ap <- fixt()$active_pathways_df
      out <- data.frame(
        Pathway          = ap$pathway_name,
        `p (prot)`       = signif(ap$p_a, 3),
        `p (RNA)`        = signif(ap$p_b, 3),
        `p (combined)`   = signif(ap$p_combined, 3),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      DT::datatable(
        out,
        rownames  = FALSE,
        selection = "single",
        options   = list(
          pageLength = 10,
          dom        = "tip",
          scrollX    = TRUE,
          columnDefs = list(list(className = "dt-right",
                                 targets = c(1, 2, 3)))
        )
      )
    }, server = TRUE)
  })
}

# ---- internal helpers ------------------------------------------------

utils::globalVariables(".data")

integration_quadrant_colors <- function() {
  c(up_up     = "#C0392B",
    down_down = "#1F4E96",
    up_down   = "#E0A030",
    down_up   = "#7A4FA0",
    ns        = "#9AA3AE")
}

integration_dual_card <- function(ns) {
  bslib::card(
    bslib::card_header(
      htmltools::tags$h3(class = "card-title", "Mirrored volcano"),
      htmltools::tags$span(class = "card-sub",
                           "x = log2FC(prot) - log2FC(RNA)")
    ),
    bslib::card_body(
      shiny::plotOutput(ns("dual"), height = "320px")
    )
  )
}

integration_scatter_card <- function(ns) {
  bslib::card(
    bslib::card_header(
      htmltools::tags$h3(class = "card-title", "log2FC concordance"),
      htmltools::tags$span(class = "card-sub",
                           "protein vs RNA \u00B7 dashed = y=x")
    ),
    bslib::card_body(
      shiny::plotOutput(ns("scatter"), height = "320px")
    )
  )
}

integration_shared_card <- function(ns) {
  bslib::card(
    bslib::card_header(
      htmltools::tags$h3(class = "card-title", "Shared pathway dotplot"),
      htmltools::tags$span(class = "card-sub",
                           "rows = pathway \u00B7 cols = omics")
    ),
    bslib::card_body(
      shiny::plotOutput(ns("shared"), height = "340px")
    )
  )
}

integration_ap_card <- function(ns) {
  bslib::card(
    bslib::card_header(
      htmltools::tags$h3(class = "card-title",
                         "ActivePathways \u00B7 combined p"),
      htmltools::tags$span(class = "card-sub",
                           "Brown's method \u00B7 Reactome")
    ),
    bslib::card_body(
      DT::DTOutput(ns("ap_table"))
    )
  )
}
