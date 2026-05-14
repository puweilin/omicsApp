#' Report view module
#'
#' Slice 2F: previews the cached analysis bundles that will feed the
#' final Rmd template. Real `omicsCore::export_report()` wiring is
#' deferred to Phase 4 — the "Generate HTML" / "Generate PDF" buttons
#' render but are inert in this slice. Closes Phase 2 (all 7 views
#' now populated by example data).
#'
#' Reference markup: `omicsApp/mockup/index.html:1057-1071`. The
#' mockup itself shows this view as a Phase-4 placeholder; the
#' fixture-backed bundle list here is slice-2F polish so the page
#' isn't a bare card.
#'
#' @param id Module namespace id.
#' @keywords internal
#' @noRd
report_view_ui <- function(id) {
  ns <- shiny::NS(id)

  qc  <- example_qc_bundle()
  dif <- example_diff_bundle()
  enr <- example_enrich_table()
  itg <- example_integration_tables()

  qc_n_samples  <- nrow(qc$results$qc_summary$missingness$sample_metrics)
  qc_n_features <- nrow(qc$results$qc_summary$missingness$feature_metrics)
  dif_n         <- nrow(dif$results$diff_result_df)
  enr_n         <- nrow(enr)
  itg_n         <- nrow(itg$concordance_df)

  htmltools::tagList(
    view_header(
      title    = "Report",
      subtitle = "Preview cached analysis bundles \u00B7 generation lands in Phase 4",
      actions  = htmltools::tagList(
        htmltools::tags$button(type = "button",
                               class = "btn btn-ghost",
                               "Generate HTML"),
        htmltools::tags$button(type = "button",
                               class = "btn btn-ghost",
                               "Generate PDF")
      )
    ),
    notice(
      title  = "Report generation is wired in Phase 4.",
      detail = paste(
        "This page previews the cached analysis bundles that will feed",
        "the final Rmd template. Clicking the header buttons is a no-op",
        "in this slice."
      ),
      kind   = "info"
    ),
    bslib::card(
      bslib::card_header(
        htmltools::tags$h3(class = "card-title", "Bundles available"),
        htmltools::tags$span(class = "card-sub",
                             "ready for the Phase-4 template")
      ),
      bslib::card_body(
        schema_row(
          ix         = "1",
          title      = "QC",
          desc       = sprintf("%d samples \u00B7 %d features",
                               qc_n_samples, qc_n_features),
          role       = "ok",
          confidence = 1.0,
          accent     = "ok",
          role_kind  = "ok"
        ),
        schema_row(
          ix         = "2",
          title      = "Differential",
          desc       = sprintf("%d features tested \u00B7 G2 vs G1", dif_n),
          role       = "limma",
          confidence = 1.0,
          accent     = "ok",
          role_kind  = "ok"
        ),
        schema_row(
          ix         = "3",
          title      = "Enrichment",
          desc       = sprintf("%d hallmark pathways", enr_n),
          role       = "GSEA",
          confidence = 1.0,
          accent     = "ok",
          role_kind  = "ok"
        ),
        schema_row(
          ix         = "4",
          title      = "Integration",
          desc       = sprintf("%d paired features \u00B7 combined p", itg_n),
          role       = "AP",
          confidence = 1.0,
          accent     = "ok",
          role_kind  = "ok"
        )
      )
    )
  )
}

#' @rdname report_view_ui
#' @keywords internal
#' @noRd
report_view_server <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {
    # Server body intentionally empty. Phase 4 wires
    # `omicsCore::export_report()` here behind the Generate buttons.
  })
}
