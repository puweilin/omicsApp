#' Import view module
#'
#' Phase 3 slice 3A: real upload + smart-parse. The view accepts a
#' single Excel / CSV / TSV / RDS file via `shiny::fileInput()`,
#' hands the path to [omicsCore::read_omics()], renders the resulting
#' `ImportReport` (per-sheet classifier table + warnings strip), and
#' exposes an [omicsCore::omics_input()] via the module's return
#' value once the user clicks "Confirm".
#'
#' Slice 3A scope:
#'   * single-file upload only (no multi-experiment merge yet — that's
#'     slice 3E once Integration needs two layers);
#'   * `omics_type` chosen via a radio (proteomics / rnaseq) since the
#'     classifier needs it to build the actual `omics_input`;
#'   * read-only schema table (no Edit / Re-detect override UI);
#'   * the Confirm button stores the input in a module-scoped
#'     `reactiveVal` and the module returns a reactive over it; the
#'     parent currently ignores the return value. Slice 3B wires it
#'     to the app-level `current_project`.
#'
#' Reference markup: `omicsApp/mockup/index.html:565-683`.
#'
#' @param id Module namespace id.
#' @keywords internal
#' @noRd
import_view_ui <- function(id) {
  ns <- shiny::NS(id)

  htmltools::tagList(
    view_header(
      title    = "Import data",
      subtitle = "Auto-detect expression matrix, sample metadata, and feature annotation"
    ),
    shiny::uiOutput(ns("steps_strip")),
    htmltools::tags$div(
      class = "row-grid r-4-8",
      import_upload_card(ns),
      import_schema_card(ns)
    )
  )
}

#' @rdname import_view_ui
#' @keywords internal
#' @noRd
import_view_server <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ---- reactive state -----------------------------------------------
    # `parsed` holds the most recent successful read_omics() return.
    # `confirmed_input` holds the input the user has signed off on.
    # The two are separate so the Confirm step is deliberate.
    parsed         <- shiny::reactiveVal(NULL)   # list(input, report) or NULL
    confirmed_input <- shiny::reactiveVal(NULL)  # omics_input or NULL

    # ---- parse on upload or omics-type change -------------------------
    # The classifier ignores omics_type, but the final omics_input()
    # construction does — so changing the radio after an upload should
    # rebuild the candidate input. Observe both.
    do_parse <- function() {
      f <- input$file
      shiny::req(f)
      omics_type <- input$omics_type %||% "proteomics"
      assay_type <- if (omics_type == "rnaseq") "raw_count" else "intensity"

      out <- tryCatch(
        omicsCore::read_omics(
          f$datapath,
          omics_type = omics_type,
          assay_type = assay_type
        ),
        error = function(e) {
          list(
            input = NULL,
            report = omicsCore::new_import_report(
              warnings = paste0("read_omics() failed: ",
                                conditionMessage(e)),
              source = f$name
            )
          )
        }
      )
      # Stamp the user-visible source name so the schema card shows
      # the original filename, not the tempfile path Shiny gave us.
      out$report$source <- f$name
      parsed(out)
      # New upload (or radio change) always resets the confirmed state
      # so the user has to re-confirm against the rebuilt input.
      confirmed_input(NULL)
    }

    shiny::observeEvent(input$file,       do_parse())
    shiny::observeEvent(input$omics_type, do_parse(), ignoreInit = TRUE)

    # ---- derived state for the UI -------------------------------------
    has_file       <- shiny::reactive(!is.null(input$file))
    parse_ok       <- shiny::reactive(!is.null(parsed()) &&
                                        !is.null(parsed()$input))
    is_confirmed   <- shiny::reactive(!is.null(confirmed_input()))

    # ---- steps strip --------------------------------------------------
    output$steps_strip <- shiny::renderUI({
      step1 <- if (has_file())       "done"   else "active"
      step2 <- if (!has_file())      "pending"
               else if (parse_ok())  "done"
               else                  "active"
      step3 <- if (is_confirmed())   "done"
               else if (parse_ok())  "active"
               else                  "pending"

      desc1 <- if (has_file()) input$file$name else "Excel / CSV / TSV / RDS"
      desc2 <- if (has_file()) {
        rep <- parsed()$report
        if (parse_ok()) {
          sprintf("%d sheet%s detected",
                  nrow(rep$sheets),
                  if (nrow(rep$sheets) == 1L) "" else "s")
        } else {
          "needs review"
        }
      } else {
        "auto-detect roles"
      }
      desc3 <- if (is_confirmed()) "omics_input ready"
               else if (parse_ok()) "click Confirm"
               else "pending"

      htmltools::tags$div(
        class = "steps",
        step_item(1L, "Upload",                  desc1, state = step1),
        step_arrow(),
        step_item(2L, "Review inferred schema",  desc2, state = step2),
        step_arrow(),
        step_item(3L, "Confirm & import",        desc3, state = step3)
      )
    })

    # ---- upload card --------------------------------------------------
    output$upload_status <- shiny::renderUI({
      if (!has_file()) {
        return(htmltools::tags$div(
          class = "muted",
          style = "font-size:12px;margin-top:8px",
          "No file selected yet."
        ))
      }
      f <- input$file
      file_row(
        name = f$name,
        meta = sprintf("%s \u00B7 %s",
                       toupper(tools::file_ext(f$name)),
                       input$omics_type %||% "proteomics"),
        size = format_file_size(f$size)
      )
    })

    # ---- schema card --------------------------------------------------
    output$schema_table <- DT::renderDT({
      rep <- parsed()$report
      shiny::req(rep)
      df <- rep$sheets
      if (is.null(df) || nrow(df) == 0L) {
        return(DT::datatable(
          data.frame(message = "No sheets to display."),
          options = list(dom = "t"), rownames = FALSE
        ))
      }
      DT::datatable(
        df[, c("name", "role", "n_rows", "n_cols",
               "confidence", "orientation", "notes"), drop = FALSE],
        rownames  = FALSE,
        selection = "none",
        options   = list(
          pageLength = 8,
          dom        = "tip",
          scrollX    = TRUE,
          columnDefs = list(list(className = "dt-right",
                                 targets = c(2, 3, 4)))
        )
      ) |>
        DT::formatRound("confidence", 2)
    }, server = TRUE)

    output$schema_warnings <- shiny::renderUI({
      rep <- parsed()$report
      if (is.null(rep) || length(rep$warnings) == 0L) return(NULL)
      htmltools::tags$div(
        style = "margin-top:12px;display:flex;flex-direction:column;gap:6px",
        lapply(rep$warnings, function(w) notice(title = w, kind = "warn"))
      )
    })

    output$schema_summary <- shiny::renderUI({
      rep <- parsed()$report
      if (is.null(rep)) return(NULL)
      sug <- rep$suggested_input
      if (length(sug) == 0L) return(NULL)
      bits <- list()
      if (!is.null(sug$matrix_sheet))
        bits$matrix <- sprintf("matrix=%s", sug$matrix_sheet)
      if (!is.null(sug$metadata_sheet))
        bits$meta <- sprintf("meta=%s", sug$metadata_sheet)
      if (!is.null(sug$feature_sheet))
        bits$feat <- sprintf("features=%s", sug$feature_sheet)
      if (!is.null(sug$orientation))
        bits$orient <- sprintf("orient=%s", sug$orientation)
      if (length(bits) == 0L) return(NULL)
      htmltools::tags$div(
        class = "muted",
        style = "font-size:12px;margin-top:6px",
        paste(unlist(bits), collapse = " \u00B7 ")
      )
    })

    # ---- confirm button gating ----------------------------------------
    output$confirm_state <- shiny::renderUI({
      if (is_confirmed()) {
        return(htmltools::tags$div(
          style = "display:flex;align-items:center;gap:8px;justify-content:flex-end",
          pill("omics_input ready", kind = "ok"),
          htmltools::tags$span(
            class = "muted", style = "font-size:12px",
            sprintf("%d features \u00D7 %d samples",
                    nrow(confirmed_input()$expr_mat),
                    ncol(confirmed_input()$expr_mat))
          )
        ))
      }
      if (!parse_ok()) {
        return(htmltools::tags$div(
          class = "muted", style = "font-size:12px;text-align:right",
          "Upload a file and pick the omics type to enable Confirm."
        ))
      }
      NULL
    })

    shiny::observeEvent(input$confirm, {
      shiny::req(parse_ok())
      confirmed_input(parsed()$input)
    })

    # ---- module return ------------------------------------------------
    shiny::reactive(confirmed_input())
  })
}

# ---- internal helpers ------------------------------------------------

# Tiny `%||%` so the module doesn't pull rlang in just for one operator.
`%||%` <- function(a, b) if (is.null(a)) b else a

format_file_size <- function(bytes) {
  if (is.null(bytes) || !is.finite(bytes)) return("")
  if (bytes < 1024)        return(sprintf("%d B", as.integer(bytes)))
  if (bytes < 1024^2)      return(sprintf("%.1f KB", bytes / 1024))
  if (bytes < 1024^3)      return(sprintf("%.1f MB", bytes / 1024^2))
  sprintf("%.1f GB", bytes / 1024^3)
}

# Chevron between two steps. Inline SVG so we don't depend on a
# specific bsicon name for what is essentially a typographic glyph.
step_arrow <- function() {
  htmltools::tags$svg(
    class    = "icon step-arrow",
    viewBox  = "0 0 24 24",
    fill     = "none",
    stroke   = "currentColor",
    width    = "16",
    height   = "16",
    htmltools::tags$path(d = "M9 6l6 6-6 6")
  )
}

import_upload_card <- function(ns) {
  bslib::card(
    bslib::card_header(
      htmltools::tags$h3(class = "card-title", "Upload"),
      htmltools::tags$span(class = "card-sub",
                           "single file \u00B7 one omics layer")
    ),
    bslib::card_body(
      shiny::fileInput(
        ns("file"),
        label = NULL,
        multiple = FALSE,
        accept = c(".xlsx", ".xls", ".csv", ".tsv", ".txt", ".rds"),
        placeholder = "Drop or browse \u2026"
      ),
      shiny::radioButtons(
        ns("omics_type"),
        label  = "Omics layer",
        choices = c("Proteomics" = "proteomics", "RNA-seq" = "rnaseq"),
        selected = "proteomics",
        inline = TRUE
      ),
      shiny::uiOutput(ns("upload_status"))
    )
  )
}

import_schema_card <- function(ns) {
  bslib::card(
    bslib::card_header(
      htmltools::tags$h3(class = "card-title", "Inferred schema"),
      htmltools::tags$span(class = "card-sub",
                           "per-sheet classification \u00B7 read-only in slice 3A")
    ),
    bslib::card_body(
      DT::DTOutput(ns("schema_table")),
      shiny::uiOutput(ns("schema_summary")),
      shiny::uiOutput(ns("schema_warnings")),
      htmltools::tags$div(
        style = "display:flex;gap:8px;justify-content:flex-end;align-items:center;margin-top:18px",
        shiny::uiOutput(ns("confirm_state"), inline = TRUE),
        shiny::actionButton(
          ns("confirm"),
          "Confirm & build omics_input",
          class = "btn btn-primary"
        )
      )
    )
  )
}
