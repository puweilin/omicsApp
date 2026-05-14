#' QC view module
#'
#' Renders the slice-2D `example_qc_bundle()` against the QC mockup:
#' a 4-card stat strip on top, a PCA scatter on the left, and a
#' missing-rate distribution on the right. Read-only — the mockup
#' has no controls here, and the bundle is built once per session.
#'
#' Reference markup: `omicsApp/mockup/index.html:686-750`.
#'
#' @param id Module namespace id.
#' @keywords internal
#' @noRd
qc_view_ui <- function(id) {
  ns <- shiny::NS(id)
  proto <- example_input("proteomics")
  n_samples  <- ncol(proto$expr_mat)
  n_features <- nrow(proto$expr_mat)

  htmltools::tagList(
    view_header(
      title    = "Quality control",
      subtitle = sprintf(
        "Proteomics \u00B7 %d features \u00B7 %d samples \u00B7 2 groups",
        n_features, n_samples
      )
    ),
    htmltools::tags$div(
      class = "stat-grid",
      stat_card(
        label  = "Samples passing",
        value  = sprintf("%d / %d", n_samples, n_samples),
        trend  = "no IQR outliers",
        accent = "ok",
        mono   = TRUE
      ),
      stat_card(
        label = "Features kept",
        value = format(n_features, big.mark = ","),
        trend = "missing-rate filter at 50%",
        mono  = TRUE
      ),
      stat_card(
        label = "Imputation",
        value = "none",
        trend = "5% NAs left visible"
      ),
      stat_card(
        label  = "Batch effect",
        value  = "low",
        trend  = "PVCA: batch < group",
        accent = "ok"
      )
    ),
    htmltools::tags$div(
      class = "row-grid r-7-5",
      bslib::card(
        bslib::card_header(
          htmltools::tags$h3(class = "card-title", "PCA"),
          htmltools::tags$span(class = "card-sub",
                               "samples projected on PC1 \u00D7 PC2")
        ),
        bslib::card_body(
          shiny::plotOutput(ns("pca"), height = "360px"),
          htmltools::tags$div(
            class = "legend",
            legend_swatch("G1", "var(--brand-600)"),
            legend_swatch("G2", "var(--omics-down)")
          )
        )
      ),
      bslib::card(
        bslib::card_header(
          htmltools::tags$h3(class = "card-title", "Missingness"),
          htmltools::tags$span(class = "card-sub",
                               "per-sample and per-feature missing rate")
        ),
        bslib::card_body(
          shiny::plotOutput(ns("missing"), height = "360px"),
          htmltools::tags$div(
            class = "muted",
            style = "font-size:12px;margin-top:6px",
            "Demo fixture: ~5% of cells set to NA at random."
          )
        )
      )
    )
  )
}

#' @rdname qc_view_ui
#' @keywords internal
#' @noRd
qc_view_server <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {
    qc_bundle <- shiny::reactive(example_qc_bundle())

    output$pca <- shiny::renderPlot(
      omicsCore::plot_qc(qc_bundle(), view = "pca", color_by = "group")
    )
    output$missing <- shiny::renderPlot(
      omicsCore::plot_qc(qc_bundle(), view = "missing")
    )
  })
}
