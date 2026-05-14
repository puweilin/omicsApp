#' Report view module
#'
#' Slice 3F: wires `omicsCore::export_report()` behind "Generate HTML"
#' and "Generate PDF" download buttons. Gates on `rmarkdown` availability
#' (buttons disabled + install hint when missing). Shows live bundle
#' status cards derived from `project$bundles`. Closes Phase 3 (full
#' happy path: upload \u2192 QC \u2192 diff \u2192 enrich \u2192 integration \u2192 report).
#'
#' Reference markup: `omicsApp/mockup/index.html:1057-1071`.
#'
#' @param id Module namespace id.
#' @param current_project Reactive yielding the live `omics_project`.
#' @keywords internal
#' @noRd
report_view_ui <- function(id) {
  ns <- shiny::NS(id)

  htmltools::tagList(
    shiny::uiOutput(ns("header")),
    shiny::uiOutput(ns("notices")),
    shiny::uiOutput(ns("bundle_cards"))
  )
}

#' @rdname report_view_ui
#' @keywords internal
#' @noRd
report_view_server <- function(id, current_project = shiny::reactiveVal(NULL)) {
  shiny::moduleServer(id, function(input, output, session) {

    have_rmd <- requireNamespace("rmarkdown", quietly = TRUE)
    ns <- session$ns

    # ---- header -------------------------------------------------------
    output$header <- shiny::renderUI({
      proj <- current_project()
      subtitle <- if (is.null(proj)) {
        "Demo fixture \u00B7 rmarkdown"
      } else {
        sprintf("%s \u00B7 %d experiment(s)",
                proj$name %||% "Project",
                length(proj$experiments))
      }
      view_header(
        title    = "Report",
        subtitle = subtitle,
        actions  = htmltools::tagList(
          shiny::downloadButton(
            ns("download_html"), "Generate HTML",
            class = if (have_rmd) "btn btn-primary" else "btn btn-ghost"
          ),
          shiny::downloadButton(
            ns("download_pdf"), "Generate PDF",
            class = if (have_rmd) "btn btn-ghost" else "btn btn-ghost"
          )
        )
      )
    })

    # ---- notices ------------------------------------------------------
    output$notices <- shiny::renderUI({
      tagged <- htmltools::tagList()
      if (!have_rmd) {
        tagged <- htmltools::tagAppendChild(
          tagged,
          notice(
            title  = "rmarkdown unavailable",
            detail = paste0(
              "Install rmarkdown to generate reports: ",
              "`install.packages(\"rmarkdown\")`. ",
              "The Generate buttons below are disabled until it is available."
            ),
            kind   = "warn"
          )
        )
      }
      proj <- current_project()
      if (is.null(proj)) {
        tagged <- htmltools::tagAppendChild(
          tagged,
          notice(
            title  = "Import data first",
            detail = paste0(
              "Generate a report after importing data and running at ",
              "least one analysis (QC, differential, enrichment, or ",
              "integration)."
            ),
            kind   = "info"
          )
        )
      } else {
        bnd <- proj$bundles
        n_bundles <- length(bnd)
        if (n_bundles == 0L) {
          tagged <- htmltools::tagAppendChild(
            tagged,
            notice(
              title  = "No analysis bundles yet",
              detail = paste0(
                "Run QC, differential, enrichment, or integration ",
                "before generating a report."
              ),
              kind   = "info"
            )
          )
        }
      }
      tagged
    })

    # ---- bundle status cards ------------------------------------------
    output$bundle_cards <- shiny::renderUI({
      proj <- current_project()

      make_row <- function(ix, title, desc, role, ready) {
        schema_row(
          ix         = ix,
          title      = title,
          desc       = desc,
          role       = role,
          confidence = if (ready) 1.0 else 0.0,
          accent     = if (ready) "ok" else "warn",
          role_kind  = if (ready) "ok" else "warn"
        )
      }

      if (is.null(proj)) {
        # Demo fixture: show example counts.
        qc  <- example_qc_bundle()
        dif <- example_diff_bundle()
        enr <- example_enrich_table()
        itg <- example_integration_tables()
        qc_n  <- nrow(qc$results$qc_summary$missingness$sample_metrics)
        dif_n <- nrow(dif$results$diff_result_df)
        enr_n <- nrow(enr)
        itg_n <- nrow(itg$concordance_df)
        rows <- list(
          make_row("1", "QC", sprintf("%d samples", qc_n),
                   "demo", TRUE),
          make_row("2", "Differential", sprintf("%d features", dif_n),
                   "demo", TRUE),
          make_row("3", "Enrichment", sprintf("%d pathways", enr_n),
                   "demo", TRUE),
          make_row("4", "Integration", sprintf("%d paired features", itg_n),
                   "demo", TRUE)
        )
      } else {
        bnd <- proj$bundles
        qc_ready <- !is.null(bnd$qc)
        diff_ready <- !is.null(bnd$diff)
        enrich_ready <- !is.null(bnd$enrich)
        int_ready <- !is.null(bnd$integration)

        qc_desc <- if (qc_ready) {
          sm <- bnd$qc$results$qc_summary$missingness$sample_metrics
          sprintf("%d samples \u00B7 %d features passed",
                  nrow(sm), nrow(bnd$qc$results$qc_summary$missingness$feature_metrics))
        } else "Not yet run"

        diff_desc <- if (diff_ready) {
          sprintf("%d features \u00B7 %s vs %s",
                  nrow(bnd$diff$results$diff_result_df),
                  bnd$diff$params$case_group %||% "case",
                  bnd$diff$params$control_group %||% "control")
        } else "Not yet run"

        enrich_desc <- if (enrich_ready) {
          sprintf("%s \u00B7 %d pathways",
                  toupper(bnd$enrich$params$type %||% "ora"),
                  nrow(bnd$enrich$results$enrich_result_df))
        } else "Not yet run"

        int_desc <- if (int_ready) {
          sprintf("%d paired features",
                  nrow(bnd$integration$results$integration_df))
        } else "Not yet run"

        rows <- list(
          make_row("1", "QC", qc_desc,
                   if (qc_ready) "QC" else "pending", qc_ready),
          make_row("2", "Differential", diff_desc,
                   if (diff_ready) bnd$diff$params$method %||% "auto" else "pending",
                   diff_ready),
          make_row("3", "Enrichment", enrich_desc,
                   if (enrich_ready) toupper(bnd$enrich$params$type %||% "ora") else "pending",
                   enrich_ready),
          make_row("4", "Integration", int_desc,
                   if (int_ready) "concordance" else "pending",
                   int_ready)
        )
      }

      bslib::card(
        bslib::card_header(
          htmltools::tags$h3(class = "card-title", "Bundles available"),
          htmltools::tags$span(class = "card-sub",
                               "ready for the report template")
        ),
        bslib::card_body(rows)
      )
    })

    # ---- download handlers --------------------------------------------
    output$download_html <- shiny::downloadHandler(
      filename = function() {
        proj <- current_project()
        nm <- if (!is.null(proj)) proj$name %||% "omics" else "omics"
        sprintf("%s_report.html", gsub("[^A-Za-z0-9_.-]+", "_", nm))
      },
      content = function(file) {
        shiny::req(have_rmd)
        proj <- current_project()
        shiny::req(proj)
        omicsCore::export_report(proj, file, format = "html",
                                 overwrite = TRUE)
      }
    )

    output$download_pdf <- shiny::downloadHandler(
      filename = function() {
        proj <- current_project()
        nm <- if (!is.null(proj)) proj$name %||% "omics" else "omics"
        sprintf("%s_report.pdf", gsub("[^A-Za-z0-9_.-]+", "_", nm))
      },
      content = function(file) {
        shiny::req(have_rmd)
        proj <- current_project()
        shiny::req(proj)
        omicsCore::export_report(proj, file, format = "pdf",
                                 overwrite = TRUE)
      }
    )
  })
}
