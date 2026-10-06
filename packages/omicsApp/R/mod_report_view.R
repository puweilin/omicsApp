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
    shiny::uiOutput(ns("bundle_cards")),
    shiny::uiOutput(ns("script_card"))
  )
}

#' @rdname report_view_ui
#' @keywords internal
#' @noRd
report_view_server <- function(id, current_project = shiny::reactiveVal(NULL)) {
  shiny::moduleServer(id, function(input, output, session) {

    have_rmd <- has_pkg("rmarkdown")
    ns <- session$ns

    # ---- header -------------------------------------------------------
    output$header <- shiny::renderUI({
      proj <- current_project()
      subtitle <- if (is.null(proj)) {
        "Demo data \u00B7 rmarkdown"
      } else {
        sprintf("%s \u00B7 %d layer(s)",
                proj$name %||% "Project",
                length(proj$experiments))
      }
      # There is nothing to report on until a project exists, and
      # do_export() stops on req(proj) when there is not -- silently, as
      # req() does. A live-looking button that produces no file and no
      # message is indistinguishable from a broken one, and that is what
      # it looked like.
      ready <- have_rmd && !is.null(proj)
      dl <- function(id, label, primary, extra_ok = TRUE,
                     extra_why = NULL) {
        ok <- ready && extra_ok
        btn <- shiny::downloadButton(
          ns(id), label,
          class = paste("btn", if (ok && primary) "btn-primary" else "btn-ghost"))
        if (ok) return(btn)
        # `disabled` on the anchor, plus the pointer-events guard,
        # because a download link is an <a> and browsers do not honour
        # the attribute on its own.
        # And out of the keyboard's reach: Tab + Enter still started a
        # PDF export that failed on the missing LaTeX.
        htmltools::tagAppendAttributes(
          btn, disabled = NA, tabindex = "-1", `aria-disabled` = "true",
          style = "pointer-events:none;opacity:0.55",
          title = if (!have_rmd) "rmarkdown is not installed"
                  else if (is.null(proj)) "Import a project first"
                  else extra_why)
      }
      view_header(
        title    = "Report",
        subtitle = subtitle,
        actions  = htmltools::tagList(
          dl("download_html", "Generate HTML", TRUE),
          # PDF needs LaTeX, which most servers do not have; the button
          # used to be live and fail after a long render.
          dl("download_pdf", "Generate PDF", FALSE, extra_ok = have_latex(),
             extra_why = "PDF needs a LaTeX installation (e.g. tinytex::install_tinytex())")
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
              title  = "No results yet",
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
          cmps <- tryCatch(omicsCore::diff_comparisons(bnd$diff),
                           error = function(e) character(0))
          if (length(cmps) > 1L) {
            sprintf("%d comparisons: %s", length(cmps),
                    paste(gsub("_vs_", " vs ", cmps), collapse = ", "))
          } else {
            sprintf("%d features \u00B7 %s",
                    nrow(bnd$diff$results$diff_result_df),
                    gsub("_vs_", " vs ", bnd$diff$params$comparison %||% "case vs control"))
          }
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
                   if (int_ready) bnd$integration$params$method %||% "integration" else "pending",
                   int_ready)
        )
      }

      bslib::card(
        bslib::card_header(
          htmltools::tags$h3(class = "card-title", "Results in the report"),
          htmltools::tags$span(class = "card-sub",
                               "ready for the report template")
        ),
        bslib::card_body(rows)
      )
    })

    # ---- download handlers --------------------------------------------
    # `export_report` is wrapped in tryCatch + showNotification so a
    # failure inside the rmarkdown render surfaces as a banner instead
    # of an opaque Shiny error overlay.
    do_export <- function(file, format) {
      shiny::req(have_rmd)
      proj <- current_project()
      shiny::req(proj)
      if (identical(format, "pdf") && !have_latex()) {
        shiny::showNotification(
          "PDF needs a LaTeX installation on the server; use Generate HTML instead.",
          type = "warning", duration = 8)
        shiny::req(FALSE)
      }
      tryCatch(
        omicsCore::export_report(proj, file, format = format,
                                 overwrite = TRUE),
        error = function(e) {
          shiny::showNotification(
            sprintf("Report export failed (%s): %s",
                    format, conditionMessage(e)),
            type = "error",
            duration = 8
          )
          stop(e)
        }
      )
    }

    output$download_html <- shiny::downloadHandler(
      filename = function() {
        proj <- current_project()
        nm <- if (!is.null(proj)) proj$name %||% "omics" else "omics"
        sprintf("%s_report.html", download_stem(nm))
      },
      content = function(file) do_export(file, "html")
    )

    output$download_pdf <- shiny::downloadHandler(
      filename = function() {
        proj <- current_project()
        nm <- if (!is.null(proj)) proj$name %||% "omics" else "omics"
        sprintf("%s_report.pdf", download_stem(nm))
      },
      content = function(file) do_export(file, "pdf")
    )

    # ---- reproducibility script ---------------------------------------
    # Shown inline rather than offered only as a download: the point of
    # this is that a reader can check what ran, and a file they have to
    # download first is a file most people will not open.
    script_lines <- shiny::reactive({
      proj <- current_project()
      if (is.null(proj)) return(NULL)
      tryCatch(
        omicsCore::export_script(proj),
        error = function(e) paste("# Could not build the script:",
                                  conditionMessage(e))
      )
    })

    output$script_card <- shiny::renderUI({
      lines <- script_lines()
      if (is.null(lines)) {
        return(bslib::card(
          bslib::card_header(
            htmltools::tags$h3(class = "card-title", "Analysis code")
          ),
          bslib::card_body(
            htmltools::tags$div(
              class = "muted", style = "font-size:13px",
              "Import data and run an analysis to see the code that produced it."
            )
          )
        ))
      }
      bslib::card(
        bslib::card_header(
          htmltools::tags$h3(class = "card-title", "Analysis code"),
          htmltools::tags$span(
            class = "card-sub",
            "The calls that produced this project, as they were resolved"
          )
        ),
        bslib::card_body(
          htmltools::tags$div(
            style = "display:flex;justify-content:flex-end;padding-bottom:8px",
            shiny::downloadButton(session$ns("download_script"), "Download .R",
                                  class = "btn btn-sm btn-primary"),
            # The script reads raw/<file>: without the files it points
            # at, a downloaded script could not run anywhere.
            shiny::downloadButton(session$ns("download_bundle"), "Script + data files",
                                  class = "btn btn-sm btn-ghost")
          ),
          htmltools::tags$pre(
            class = "text-mono",
            style = paste("max-height:460px;overflow:auto;font-size:12px",
                          "background:var(--bg);padding:12px;border:1px solid",
                          "var(--border);border-radius:6px", sep = ";"),
            paste(lines, collapse = "\n")
          )
        )
      )
    })

    # The script, the archived uploads it reads (raw/...) and a README,
    # as one archive that runs where it is unpacked: setwd() to the
    # folder and source() the script.
    output$download_bundle <- shiny::downloadHandler(
      filename = function() {
        proj <- current_project()
        nm <- if (!is.null(proj)) proj$name %||% "omics" else "omics"
        sprintf("%s_analysis.%s", download_stem(nm), if (nzchar(Sys.which("zip"))) "zip" else "tar.gz")
      },
      content = function(file) {
        proj <- current_project()
        lines <- script_lines()
        shiny::req(proj, lines)
        dir <- file.path(tempfile("omicsapp-bundle-"), "analysis")
        dir.create(file.path(dir, "raw"), recursive = TRUE)
        on.exit(unlink(dirname(dir), recursive = TRUE), add = TRUE)
        writeLines(lines, file.path(dir, "analysis.R"))
        srcs <- unique(stats::na.omit(unlist(lapply(proj$experiments, function(e)
          c(e$source_path %||% NA_character_, e$sample_sheet_path %||% NA_character_)))))
        srcs <- srcs[file.exists(srcs)]
        file.copy(srcs, file.path(dir, "raw", basename(srcs)))
        writeLines(c(
          "Reproducing this analysis",
          "",
          "1. Install R and the omicsCore package.",
          "2. In R: setwd(\"<this folder>\"); source(\"analysis.R\")",
          "",
          if (length(srcs)) c("Data files in raw/:", paste0("  ", basename(srcs)))
          else "No data file was archived with this project; edit the read_omics() path(s) in analysis.R."),
          file.path(dir, "README.txt"))
        old <- setwd(dirname(dir))
        on.exit(setwd(old), add = TRUE)
        if (nzchar(Sys.which("zip"))) {
          utils::zip(file, "analysis", flags = "-rq9X", zip = Sys.which("zip"))
        } else {
          utils::tar(file, "analysis", compression = "gzip", tar = "internal")
        }
      }
    )

    output$download_script <- shiny::downloadHandler(
      filename = function() {
        proj <- current_project()
        nm <- if (!is.null(proj)) proj$name %||% "omics" else "omics"
        sprintf("%s_analysis.R", download_stem(nm))
      },
      content = function(file) {
        lines <- script_lines()
        shiny::req(lines)
        writeLines(lines, file)
      }
    )
  })
}

# Whether a PDF report can be built here.
have_latex <- function() {
  nzchar(Sys.which("pdflatex")) || nzchar(Sys.which("xelatex")) ||
    # Looked up rather than called with `::`: tinytex is not a dependency.
    (has_pkg("tinytex") && isTRUE(tryCatch(
      getExportedValue("tinytex", "is_tinytex")(), error = function(e) FALSE)))
}

# A file name stem for downloads: the store's slug, which keeps Chinese
# and other non-ASCII names ("\u4e2d\u6587\u9879\u76ee" came out as "__report.html").
download_stem <- function(nm) {
  st <- project_slug(nm %||% "")
  if (is.na(st) || !nzchar(st)) "project" else st
}
