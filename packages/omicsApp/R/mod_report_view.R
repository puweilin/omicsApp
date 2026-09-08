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
    shiny::uiOutput(ns("grid_card")),
    shiny::uiOutput(ns("sheet_card")),
    shiny::uiOutput(ns("bundle_cards")),
    shiny::uiOutput(ns("script_card"))
  )
}

#' @rdname report_view_ui
#' @keywords internal
#' @noRd
report_view_server <- function(id, current_project = shiny::reactiveVal(NULL),
                               diff_thresholds = shiny::reactive(NULL)) {
  shiny::moduleServer(id, function(input, output, session) {

    # One workbook of tables, not a rendered document.
    #
    # What a reader does next with a finished project is open its tables:
    # sort the differential by p, filter the pathways, paste a column into
    # a figure. HTML and PDF answered "what did this run do" in a form
    # nobody could carry into the next step, and the PDF half of that cost
    # a LaTeX toolchain, a CJK font stack and three rounds of debugging.
    have_xlsx <- has_pkg("openxlsx")
    ns <- session$ns

    # ---- header -------------------------------------------------------
    output$header <- shiny::renderUI({
      proj <- current_project()
      n_sheets <- length(sheet_names())
      subtitle <- if (is.null(proj)) {
        "Demo fixture \u00B7 one sheet per analysis"
      } else {
        sprintf("%s \u00B7 %d experiment(s) \u00B7 %d sheet(s)",
                proj$name %||% "Project",
                length(proj$experiments),
                n_sheets)
      }
      # There is nothing to export until a project exists, and do_export()
      # stops on req(proj) when there is not -- silently, as req() does. A
      # live-looking button that produces no file and no message is
      # indistinguishable from a broken one, and that is what it looked
      # like. A project with no analyses is the same case: the workbook
      # would have nothing in it, and export_workbook() says so.
      why <- if (!have_xlsx) "openxlsx is not installed"
             else if (is.null(proj)) "Import a project first"
             else if (n_sheets == 0L) "Run an analysis first"
             else NULL
      btn <- shiny::downloadButton(
        ns("download_xlsx"), "Generate xlsx",
        class = paste("btn", if (is.null(why)) "btn-primary" else "btn-ghost"))
      if (!is.null(why)) {
        # `disabled` on the anchor, plus the pointer-events guard, because
        # a download link is an <a> and browsers do not honour the
        # attribute on its own.
        btn <- htmltools::tagAppendAttributes(
          btn, disabled = NA,
          style = "pointer-events:none;opacity:0.55", title = why)
      }
      view_header(title = "Report", subtitle = subtitle, actions = btn)
    })

    # ---- notices ------------------------------------------------------
    # ---- the whole enrichment grid --------------------------------------
    #
    # The Enrich view runs one combination at a time, which is what you
    # want while looking at a result and not what you want when you are
    # finished: everything ORA and GSEA have to say about a differential
    # took fourteen visits to that view. This runs the grid once.
    #
    # A button rather than something the download does on its own. It is
    # tens of seconds of real work -- forty to sixty here, longer on a
    # single-core container -- and a download that blocks that long with
    # nothing on screen is indistinguishable from one that has hung.
    # Through run_async() like every other analysis, so it gets the
    # progress bar and the other sessions stay responsive.
    grid_bundle <- shiny::reactiveVal(NULL)
    grid_error <- shiny::reactiveVal(NULL)

    # A grid the project already carries -- restored, or run earlier this
    # session -- so the button does not have to be pressed twice for one
    # answer.
    shiny::observe({
      proj <- current_project()
      if (is.null(proj)) return()
      saved <- proj$bundles$enrich_all
      if (is.null(saved) || is.null(saved$provenance_binding)) return()
      if (!omicsCore::result_inputs_current(saved$provenance_binding, proj)) {
        return()
      }
      grid_bundle(saved)
    }, priority = -20L)

    shiny::observeEvent(input$run_grid, {
      proj <- current_project()
      diff <- proj$bundles$diff
      if (is.null(diff)) {
        grid_error("Run a differential first: enrichment reads its results.")
        return()
      }
      launch <- capture_launch(
        proj, "enrich_all",
        tags = diff$input_info$layer_tag %||% character(0),
        upstream = diff$provenance_binding$analysis_id %||% character(0),
        params = list(grid = "all"),
        activation = NA_integer_)
      # The path is chosen here and handed over, not chosen in the
      # worker: under `multisession` the worker has its own tempdir(),
      # so a name it picked would be one this side could not find.
      progress_path <- tempfile("omicsApp-grid-", fileext = ".progress")
      # The writer is built here and travels as a value. Built inside the
      # worker instead it would have to carry progress_writer() across,
      # and detached_call() re-parents to baseenv() -- where a package
      # internal is not.
      report_progress <- progress_writer(progress_path)
      # The cut the Enrich view is running at, so "all of them" means all
      # of them at the settings on screen rather than at whatever the
      # defaults happen to be. `output_p_cutoff = 1` for the same reason
      # the single run passes it: which genes go *in* and how much of the
      # pathway table comes back are two questions, and answering both
      # with one number stores only what the display threshold admitted.
      thr <- tryCatch(diff_thresholds(), error = function(e) NULL)
      grid_args <- list(output_p_cutoff = 1)
      if (!is.null(thr)) {
        grid_args$p_cutoff <- thr$p_cutoff
        grid_args$p_preference <- thr$p_preference
        grid_args$effect_cutoff <- thr$effect_cutoff
      }
      run_async(
        detached_call(
          function() {
            on.exit(unlink(progress_path), add = TRUE)
            do.call(omicsCore::run_enrichment_grid,
                    c(list(diff_bundle = bundle,
                           on_progress = report_progress), grid_args))
          },
          bundle = bundle_for_worker(diff),
          grid_args = grid_args,
          progress_path = progress_path,
          report_progress = report_progress),
        on_success = function(result) {
          unlink(progress_path)
          grid_error(NULL)
          grid_bundle(stamp_from_launch(result, launch))
        },
        on_error = function(msg) {
          unlink(progress_path)
          grid_bundle(NULL)
          grid_error(msg)
        },
        progress_poll = progress_reader(progress_path),
        message = "Running every enrichment..."
      )
    })

    # The sheets the workbook would hold, asked of omicsCore rather than
    # worked out again here. The header counts them, the card lists them
    # and the download writes them; three routes to the same answer is
    # three chances for the page to describe a file it did not produce.
    sheet_names <- shiny::reactive({
      proj <- current_project()
      if (is.null(proj)) return(character(0))
      tryCatch(names(omicsCore::workbook_sheets(proj)),
               error = function(e) character(0))
    })

    output$notices <- shiny::renderUI({
      tagged <- htmltools::tagList()
      if (!have_xlsx) {
        tagged <- htmltools::tagAppendChild(
          tagged,
          notice(
            title  = "openxlsx unavailable",
            detail = paste0(
              "Install openxlsx to export results: ",
              "`install.packages(\"openxlsx\")`. ",
              "The Generate button below is disabled until it is available."
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
    output$grid_card <- shiny::renderUI({
      proj <- current_project()
      have_diff <- !is.null(proj) && !is.null(proj$bundles$diff)
      g <- grid_bundle()
      err <- grid_error()

      status <- if (!is.null(err)) {
        htmltools::tags$div(class = "muted", style = "color:var(--danger)", err)
      } else if (!is.null(g)) {
        lg <- g$results$grid_log
        htmltools::tags$div(
          class = "muted", style = "font-size:12px",
          sprintf("%s combinations run · %s returned rows · %s empty · %s failed",
                  nrow(lg), sum(lg$status == "ok"),
                  sum(lg$status == "empty"), sum(lg$status == "failed")))
      } else {
        htmltools::tags$div(
          class = "muted", style = "font-size:12px",
          "ORA in three directions and GSEA, across every database. ",
          "Tens of seconds; the workbook picks the results up when it ",
          "finishes.")
      }

      btn <- shiny::actionButton(
        session$ns("run_grid"),
        if (is.null(g)) "Run every enrichment" else "Run them again",
        class = "btn btn-sm btn-primary")
      if (!have_diff) {
        btn <- htmltools::tagAppendAttributes(
          btn, disabled = NA, style = "pointer-events:none;opacity:0.55",
          title = "Run a differential first")
      }

      bslib::card(
        bslib::card_header(
          htmltools::tags$h3(class = "card-title", "Enrichment, all of it"),
          htmltools::tags$span(class = "card-sub",
                               "one pass instead of fourteen visits")
        ),
        bslib::card_body(btn, status)
      )
    })

    # What the workbook will contain, listed from the workbook itself.
    #
    # The card used to describe the four analyses and the download wrote
    # whatever it wrote; the two agreed only as long as someone kept them
    # agreeing. This asks omicsCore the same question the download asks,
    # so a sheet that appears here is a sheet in the file.
    output$sheet_card <- shiny::renderUI({
      proj <- current_project()
      sheets <- if (is.null(proj)) NULL else
        tryCatch(omicsCore::workbook_sheets(proj), error = function(e) NULL)

      body <- if (is.null(proj)) {
        htmltools::tags$div(
          class = "muted", style = "font-size:13px",
          "Import data and run an analysis to see what the workbook ",
          "would contain.")
      } else if (length(sheets) == 0L) {
        htmltools::tags$div(
          class = "muted", style = "font-size:13px",
          "Nothing to export yet. Run QC, a differential, an enrichment ",
          "or an integration.")
      } else {
        htmltools::tagList(
          lapply(seq_along(sheets), function(i) {
            df <- sheets[[i]]
            schema_row(
              ix    = as.character(i),
              title = names(sheets)[[i]],
              desc  = sprintf("%s rows · %d columns",
                              format(nrow(df), big.mark = ","), ncol(df)),
              role  = "sheet",
              confidence = 1.0, accent = "ok", role_kind = "ok")
          })
        )
      }
      bslib::card(
        bslib::card_header(
          htmltools::tags$h3(class = "card-title", "Workbook contents"),
          htmltools::tags$span(class = "card-sub",
                               "one sheet per analysis")
        ),
        bslib::card_body(body)
      )
    })

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
    # `export_report` is wrapped in tryCatch + showNotification so a
    # failure inside the rmarkdown render surfaces as a banner instead
    # of an opaque Shiny error overlay.
    output$download_xlsx <- shiny::downloadHandler(
      filename = function() {
        proj <- current_project()
        nm <- if (!is.null(proj)) proj$name %||% "omics" else "omics"
        sprintf("%s_results.xlsx", gsub("[^A-Za-z0-9_.-]+", "_", nm))
      },
      content = function(file) {
        shiny::req(have_xlsx)
        proj <- current_project()
        shiny::req(proj)
        tryCatch(
          omicsCore::export_workbook(proj, file, overwrite = TRUE),
          error = function(e) {
            shiny::showNotification(
              sprintf("Export failed: %s", conditionMessage(e)),
              type = "error", duration = 8)
            stop(e)
          }
        )
      }
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
                                  class = "btn btn-sm btn-primary")
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

    output$download_script <- shiny::downloadHandler(
      filename = function() {
        proj <- current_project()
        nm <- if (!is.null(proj)) proj$name %||% "omics" else "omics"
        sprintf("%s_analysis.R", gsub("[^A-Za-z0-9_.-]+", "_", nm))
      },
      content = function(file) {
        lines <- script_lines()
        shiny::req(lines)
        writeLines(lines, file)
      }
    )

    # Returned so app_server() can file it into the project, the same way
    # every other view's result is filed. The workbook then reads it from
    # the project rather than from this module, so a restored project
    # exports its grid without anyone pressing anything.
    shiny::reactive(grid_bundle())
  })
}
