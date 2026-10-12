# The upload side of the import view: the value-scale picker, the steps
# strip, and the upload and schema cards.
#
# Split out of mod_import_view.R. These are plain functions called from
# inside import_view_server()'s moduleServer(), not modules, so every
# input and output keeps its id ("import-assay_type", "import-schema_table").

import_scale_server <- function(input, output, session, parsed, confirmed_input, parse_ok) {
  ns <- session$ns

  # ---- assay type: inferred, then owned by the user -----------------
  # The picker is rendered only once there is a matrix to infer from, so
  # the default it shows is a statement about this file rather than a
  # blanket guess.
  output$assay_type_picker <- shiny::renderUI({
    if (!parse_ok()) return(NULL)
    omics_type <- input$omics_type %||% "proteomics"
    choices <- omicsCore::SUPPORTED_ASSAY_TYPES[[omics_type]]
    if (is.null(choices)) return(NULL)

    chosen <- parsed()$input$assay_type
    guess <- parsed()$assay_guess
    htmltools::tagList(
      shiny::selectInput(
        ns("assay_type"),
        label = "Value scale",
        choices = stats::setNames(choices, gsub("_", " ", choices)),
        selected = chosen
      ),
      # Why this was preselected, so the user can judge the guess rather
      # than take it on trust; and, for RNA-seq, what it decides.
      htmltools::tags$div(
        class = "muted", style = "font-size:12px;margin:-8px 0 10px",
        if (!is.null(guess$reason)) {
          if (identical(chosen, guess$assay_type)) {
            paste("Guessed from the values:", guess$reason)
          } else {
            sprintf("Guessed '%s' from the values; you changed it.",
                    gsub("_", " ", guess$assay_type))
          }
        },
        if (identical(omics_type, "rnaseq") && !identical(chosen, "raw_count")) {
          paste(" DESeq2 and edgeR need read counts, so Differential will offer",
                "limma, the t-test and lm for this layer.")
        }
      )
    )
  })

  # Relabelling does not re-read the file; it rewrites the field the
  # analysis backends read, and the fingerprint that decides whether a
  # re-import counts as new data.
  shiny::observeEvent(input$assay_type, {
    cand <- parsed()
    shiny::req(cand, cand$input)
    if (identical(cand$input$assay_type, input$assay_type)) return()

    cand$input$assay_type <- input$assay_type
    f <- input$file
    if (!is.null(f)) {
      cand$input$source_fingerprint <- input_fingerprint(
        f$datapath, cand$input$omics_type, input$assay_type)
    }
    parsed(cand)
    confirmed_input(NULL)
  }, ignoreInit = TRUE)

  # omicsCore signals a scale mismatch with warning(), which never reaches a
  # Shiny user. Surfaced here, because getting this wrong is silent
  # everywhere else: limma would run on untransformed intensities and still
  # return a full result table.
  output$scale_notice <- shiny::renderUI({
    cand <- parsed()
    if (!parse_ok()) return(NULL)
    chosen <- input$assay_type %||% cand$input$assay_type
    if (is.null(chosen)) return(NULL)

    probe <- cand$input
    probe$assay_type <- chosen
    msg <- NULL
    withCallingHandlers(
      omicsCore::check_assay_scale(probe),
      warning = function(w) {
        msg <<- conditionMessage(w)
        invokeRestart("muffleWarning")
      }
    )
    if (is.null(msg)) return(NULL)
    htmltools::tags$div(
      style = "margin-top:-8px;margin-bottom:8px",
      notice(title = msg, kind = "warn")
    )
  })
  invisible()
}

import_status_server <- function(input, output, session, parsed, has_file, parse_ok,
                                 is_confirmed) {
  # ---- steps strip --------------------------------------------------
  output$steps_strip <- shiny::renderUI({
    step1 <- if (has_file())       "done"   else "active"
    step2 <- if (!has_file())      "pending"
             else if (parse_ok())  "done"
             else                  "active"
    step3 <- if (is_confirmed())   "done"
             else if (parse_ok())  "active"
             else                  "pending"

    desc1 <- if (!has_file()) "Excel / CSV / TSV / RDS"
             else if (length(input$file$name) > 1L) sprintf("%d files", length(input$file$name))
             else input$file$name
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
    desc3 <- if (is_confirmed()) "layer imported"
             else if (parse_ok()) "click Import this layer"
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
    if (length(f$name) > 1L) {
      return(file_row(
        name = sprintf("%d files: %s%s", length(f$name),
                       paste(utils::head(f$name, 3L), collapse = ", "),
                       if (length(f$name) > 3L) ", \u2026" else ""),
        meta = sprintf("one sample each \u00B7 %s", input$omics_type %||% "proteomics"),
        size = format_file_size(sum(f$size))
      ))
    }
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
  invisible()
}

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
                           "one omics layer \u00B7 data file + optional sample sheet")
    ),
    bslib::card_body(
      # Several files only for per-sample quantification output (Salmon's
      # quant.sf, RSEM's .results, kallisto's abundance.tsv), which is
      # merged into one layer; anything else is one file.
      shiny::fileInput(
        ns("file"),
        label = NULL,
        multiple = TRUE,
        accept = c(".xlsx", ".xlsm", ".xls", ".csv", ".tsv", ".txt", ".gz", ".rds",
                   ".sf", ".results"),
        placeholder = "Drop or browse \u2026"
      ),
      htmltools::tags$div(
        class = "muted", style = "font-size:12px;margin:-8px 0 10px",
        "Salmon, RSEM or kallisto output: select every sample's file at once."
      ),
      # Optional: for a matrix that carries no sample information, such as
      # a featureCounts table or a CSV export. Without it there were no
      # groups to compare and nowhere to add them.
      shiny::fileInput(
        ns("sample_file"),
        label = "Sample sheet (optional)",
        multiple = FALSE,
        accept = c(".xlsx", ".xlsm", ".xls", ".csv", ".tsv", ".txt", ".gz"),
        placeholder = "samples.csv: one row per sample, with a group column"
      ),
      shiny::radioButtons(
        ns("omics_type"),
        label  = "Omics layer",
        choices = c("Proteomics" = "proteomics", "RNA-seq" = "rnaseq"),
        selected = "proteomics",
        inline = TRUE
      ),
      # What scale the numbers are on is not recoverable from the file, and
      # every analysis backend reads it without re-deriving anything. The
      # guess is filled in from the data; this is where it gets corrected.
      shiny::uiOutput(ns("assay_type_picker")),
      shiny::uiOutput(ns("scale_notice")),
      shiny::uiOutput(ns("omics_hint")),
      shiny::uiOutput(ns("normalize_controls")),
      shiny::uiOutput(ns("upload_status")),

      # The classifier reads whatever a vendor sent; it cannot invent the
      # grouping, and without that there is nothing for Differential to
      # compare. A file that already parses is a better starting point
      # than a description of one.
      htmltools::tags$hr(style = "margin:14px 0 10px"),
      htmltools::tags$div(
        class = "muted", style = "font-size:12px;margin-bottom:6px",
        "Not sure of the format? Start from a template."
      ),
      htmltools::tags$div(
        style = "display:flex;gap:8px;flex-wrap:wrap",
        shiny::downloadButton(ns("template_proteomics"), "Proteomics template",
                              class = "btn btn-ghost btn-sm"),
        shiny::downloadButton(ns("template_rnaseq"), "RNA-seq template",
                              class = "btn btn-ghost btn-sm")
      )
    )
  )
}

import_schema_card <- function(ns) {
  bslib::card(
    bslib::card_header(
      htmltools::tags$h3(class = "card-title", "Inferred schema"),
      htmltools::tags$span(class = "card-sub",
                           "per-sheet classification \u00B7 correct it below")
    ),
    bslib::card_body(
      DT::DTOutput(ns("schema_table")),
      shiny::uiOutput(ns("schema_summary")),
      shiny::uiOutput(ns("schema_warnings")),
      htmltools::tags$div(
        style = "display:flex;gap:8px;justify-content:flex-end;align-items:center;margin-top:18px",
        shiny::uiOutput(ns("confirm_state"), inline = TRUE),
        shiny::uiOutput(ns("layer_name_ui"), inline = TRUE),
        shiny::actionButton(
          ns("confirm"),
          "Import this layer",
          class = "btn btn-primary"
        )
      )
    )
  )
}
