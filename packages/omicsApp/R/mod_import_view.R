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
    ),
    import_confirm_card(ns)
  )
}

#' @rdname import_view_ui
#' @param current_project Reactive yielding the live `omics_project` or
#'   `NULL`. Read only, to tell a first import from one that would
#'   replace a layer other analyses were computed on.
#' @keywords internal
#' @noRd
import_view_server <- function(id,
                               current_project = shiny::reactiveVal(NULL)) {
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
    # Roles the user has overridden, as read_omics() wants them. Kept apart
    # from `parsed` because they have to survive the re-parse they trigger.
    role_overrides <- shiny::reactiveVal(NULL)

    # Bumped on every new file. Shiny keeps an input's value across re-renders
    # of the control, so without this the role dropdowns would still hold the
    # previous workbook's answers and the observer below would apply them to
    # the new one the moment it parsed. Folding the generation into the input
    # ids means a new file starts with genuinely empty controls.
    parse_gen <- shiny::reactiveVal(0L)

    do_parse <- function() {
      f <- input$file
      shiny::req(f)
      omics_type <- input$omics_type %||% "proteomics"
      # Parse first with the modality default, then re-label from the values
      # once there is a matrix to look at. read_omics() needs *an* assay_type,
      # and the data it would be inferred from does not exist until it returns.
      assay_type <- if (omics_type == "rnaseq") "raw_count" else "raw_intensity"

      out <- tryCatch(
        omicsCore::read_omics(
          f$datapath,
          omics_type = omics_type,
          assay_type = assay_type,
          sheet_roles = role_overrides()
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
      if (!is.null(out$input)) {
        inferred <- omicsCore::infer_assay_type(out$input$expr_mat, omics_type)
        if (!is.na(inferred)) {
          out$input$assay_type <- inferred
          assay_type <- inferred
        }
        # Fingerprint the upload so Confirm can tell a genuinely new
        # dataset from the same file picked twice.
        out$input$source_fingerprint <-
          input_fingerprint(f$datapath, omics_type, assay_type)
      }
      parsed(out)
      # New upload (or radio change) always resets the confirmed state
      # so the user has to re-confirm against the rebuilt input.
      confirmed_input(NULL)
    }

    shiny::observeEvent(input$file, {
      # A new file makes the previous sheet assignment meaningless
      role_overrides(NULL)
      parse_gen(shiny::isolate(parse_gen()) + 1L)
      do_parse()
    })
    shiny::observeEvent(input$omics_type, do_parse(), ignoreInit = TRUE)

    # ---- assay type: inferred, then owned by the user -----------------
    # The picker is rendered only once there is a matrix to infer from, so
    # the default it shows is a statement about this file rather than a
    # blanket guess.
    output$assay_type_picker <- shiny::renderUI({
      if (!parse_ok()) return(NULL)
      omics_type <- input$omics_type %||% "proteomics"
      choices <- omicsCore::SUPPORTED_ASSAY_TYPES[[omics_type]]
      if (is.null(choices)) return(NULL)

      shiny::selectInput(
        ns("assay_type"),
        label = "Value scale",
        choices = stats::setNames(choices, gsub("_", " ", choices)),
        selected = parsed()$input$assay_type
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

    # ---- confirmation card --------------------------------------------
    output$confirm_shape <- shiny::renderUI({
      cand <- parsed()
      if (is.null(cand)) return(NULL)
      confirm_shape_ui(cand$input, cand$report)
    })

    output$confirm_roles <- shiny::renderUI({
      cand <- parsed()
      if (is.null(cand)) return(NULL)
      confirm_roles_ui(ns, cand$report, gen = parse_gen())
    })

    # One observer per sheet row, created after the report is known. The
    # dropdowns are rendered by confirm_roles_ui(), so the number of them is
    # data-dependent; observers are registered once and read whatever exists.
    shiny::observe({
      cand <- parsed()
      shiny::req(cand, cand$report$sheets)
      sheets <- cand$report$sheets

      gen <- parse_gen()
      chosen <- vapply(seq_len(nrow(sheets)), function(i) {
        val <- input[[role_input_id(gen, i)]]
        if (is.null(val)) NA_character_ else val
      }, character(1))
      if (all(is.na(chosen))) return()

      current <- sheets$role
      changed <- !is.na(chosen) & chosen != current
      if (!any(changed)) return()

      # Carry every explicit choice, not just the changed one: re-parsing
      # rebuilds the table from the classifier, and an earlier override would
      # otherwise be undone by the next one.
      overrides <- stats::setNames(chosen[!is.na(chosen)],
                                   sheets$name[!is.na(chosen)])
      role_overrides(overrides)
      do_parse()
    })

    output$confirm_matrix_preview <- shiny::renderTable({
      cand <- parsed()
      shiny::req(cand, cand$input)
      preview_matrix(cand$input$expr_mat)
    }, striped = TRUE, spacing = "xs", width = "100%", digits = 2)

    output$confirm_meta_preview <- shiny::renderTable({
      cand <- parsed()
      shiny::req(cand, cand$input)
      preview_metadata(cand$input$meta_df)
    }, striped = TRUE, spacing = "xs", width = "100%")

    # ---- normalization, applied on commit -----------------------------
    # This belongs to import rather than QC because it is what turns a file
    # into something the analysis backends can read: limma applies no
    # transform of its own, so an un-normalized layer means limma runs on raw
    # instrument output. The legacy framework normalized in its
    # data-input layer for the same reason; that layer is what did not survive
    # the port into this package.
    #
    # RNA-seq is deliberately excluded: DESeq2 and edgeR model raw counts
    # directly, and the t-test / lm backends log-transform "raw_count"
    # themselves.
    normalizable <- shiny::reactive({
      if (!parse_ok()) return(FALSE)
      chosen <- input$assay_type %||% parsed()$input$assay_type
      identical(parsed()$input$omics_type, "proteomics") &&
        !is.null(chosen) &&
        !chosen %in% omicsCore::LOG_SCALE_ASSAY_TYPES
    })

    output$normalize_controls <- shiny::renderUI({
      if (!parse_ok()) return(NULL)
      if (!normalizable()) {
        if (!identical(parsed()$input$omics_type, "proteomics")) return(NULL)
        return(htmltools::tags$div(
          class = "muted",
          style = "font-size:12px;margin-bottom:10px",
          "Already on a transformed scale \u2014 nothing to normalize."
        ))
      }
      htmltools::tagList(
        shiny::checkboxInput(
          ns("normalize"),
          label = "Normalize on import",
          value = TRUE
        ),
        shiny::conditionalPanel(
          condition = "input.normalize",
          ns = ns,
          shiny::selectInput(
            ns("normalize_method"),
            label = "Method",
            choices = c("vsn (variance stabilising)" = "vsn", "log2" = "log2"),
            selected = "vsn"
          )
        )
      )
    })

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

    # ---- confirm, guarded when it would replace live data -------------
    # Committing an input whose omics_type already exists in the project
    # replaces that layer, and every analysis computed on it becomes
    # meaningless. Three cases, only the last of which is destructive:
    #
    #   no such layer yet     -> commit straight away
    #   same file re-selected -> nothing changed, keep the results
    #   different data        -> ask first, then commit and clear
    #
    # The middle case is why the fingerprint is worth carrying: a
    # mis-click that re-picks the same file must not cost the user an
    # afternoon of analysis.
    layer_being_replaced <- function(cand) {
      proj <- current_project()
      if (is.null(proj) || is.null(cand)) return(NULL)
      tag <- cand$omics_type %||% "experiment"
      if (!tag %in% names(proj$experiments)) return(NULL)
      proj$experiments[[tag]]
    }

    # Archive the upload alongside the parsed input. Only on commit:
    # parsing happens on every file pick and radio change, most of which
    # the user never confirms. Archiving failing must not stop the
    # import, so its outcome is surfaced but not acted on.
    # Which normalization the current controls would apply. "none" when the
    # layer is not normalizable or the box is unticked -- both have to be
    # distinguishable in the fingerprint from an actual method.
    pending_normalize_method <- function() {
      if (normalizable() && isTRUE(input$normalize %||% TRUE)) {
        input$normalize_method %||% "vsn"
      } else {
        "none"
      }
    }

    # The identity the data would have once committed. Computed before the
    # replace check as well as inside commit(), so that re-picking the same
    # file with the same settings still reads as "nothing changed" rather than
    # as new data.
    stamp_fingerprint <- function(cand) {
      f <- input$file
      if (is.null(f)) return(cand)
      cand$source_fingerprint <- input_fingerprint(
        f$datapath, cand$omics_type, cand$assay_type,
        normalize = pending_normalize_method())
      cand
    }

    commit <- function(cand) {
      method <- pending_normalize_method()
      do_normalize <- !identical(method, "none")
      cand <- stamp_fingerprint(cand)
      f <- input$file

      if (do_normalize) {
        normalized <- tryCatch(
          suppressMessages(omicsCore::normalize_omics(cand, method = method)),
          error = function(e) e
        )
        if (inherits(normalized, "error")) {
          # Importing raw and telling the user beats importing something the
          # backends will silently mistreat
          shiny::showNotification(
            paste0("Normalization failed, importing unnormalized: ",
                   conditionMessage(normalized)),
            type = "error", duration = 12
          )
        } else {
          # normalize_omics() keeps the pre-normalization matrix in raw_mat
          normalized$source_fingerprint <- cand$source_fingerprint
          cand <- normalized
          shiny::showNotification(
            sprintf("Normalized with %s; values are now '%s'.",
                    method, cand$assay_type),
            type = "message", duration = 6
          )
        }
      }

      if (!is.null(f)) {
        res <- store_raw_upload(f$datapath, f$name, cand$source_fingerprint)
        if (isTRUE(res$ok)) {
          # Recorded so `omicsCore::export_script()` can point its
          # read_omics() line at a file that exists, which is the
          # difference between a script that runs and one that only
          # documents what was run.
          cand$source_path <- res$path
        } else if (grepl("quota", res$message, fixed = TRUE)) {
          shiny::showNotification(res$message, type = "warning", duration = 8)
        }
      }
      confirmed_input(cand)
    }

    shiny::observeEvent(input$confirm, {
      shiny::req(parse_ok())
      cand <- stamp_fingerprint(parsed()$input)
      existing <- layer_being_replaced(cand)

      if (is.null(existing)) {
        commit(cand)
        return()
      }
      if (fingerprints_match(existing, cand)) {
        shiny::showNotification(
          "That is the file already loaded \u2014 nothing to re-import.",
          type = "message"
        )
        return()
      }
      shiny::showModal(
        replace_layer_modal(ns, cand$omics_type %||% "experiment",
                            current_project())
      )
    })

    shiny::observeEvent(input$confirm_replace, {
      shiny::removeModal()
      shiny::req(parse_ok())
      commit(parsed()$input)
    })

    # ---- module return ------------------------------------------------
    shiny::reactive(confirmed_input())
  })
}

# ---- internal helpers ------------------------------------------------

# Tiny `%||%` so the module doesn't pull rlang in just for one operator.
`%||%` <- function(a, b) if (is.null(a)) b else a

# Identity of an upload. The file digest alone is not enough: the same
# workbook imported as proteomics and as RNA-seq yields two different
# inputs, so the parse settings are part of what makes it "the same".
# `normalize` is in here for the same reason -- re-importing the same file
# with normalization turned off produces different numbers, and that has to
# read as new data rather than as "nothing changed".
input_fingerprint <- function(path, omics_type, assay_type, normalize = NULL) {
  digest <- unname(tools::md5sum(path))
  if (is.na(digest)) return(NULL)
  paste(digest, omics_type %||% "", assay_type %||% "",
        if (is.null(normalize)) "" else as.character(normalize), sep = ":")
}

# Both sides must carry a fingerprint for a match to mean anything. An
# input built straight from matrices has none, and two missing values
# are not evidence of sameness — treat that as "assume it changed",
# which costs a confirmation click rather than an afternoon of results.
fingerprints_match <- function(existing, candidate) {
  a <- existing$source_fingerprint
  b <- candidate$source_fingerprint
  !is.null(a) && !is.null(b) && identical(a, b)
}

# Human-readable names for the bundles about to be discarded, so the
# dialog names what is at stake rather than saying "analyses".
BUNDLE_LABELS <- c(
  qc          = "Quality control",
  diff        = "Differential analysis",
  enrich      = "Pathway enrichment",
  integration = "Multi-omics integration"
)

replace_layer_modal <- function(ns, tag, project) {
  bundles <- names(project$bundles %||% list())
  losing <- if (length(bundles) == 0L) {
    htmltools::tags$p(
      class = "muted",
      "No analyses have been run on the current data yet."
    )
  } else {
    htmltools::tagList(
      htmltools::tags$p("These results were computed on the current data ",
                        "and will be cleared:"),
      htmltools::tags$ul(
        lapply(bundles, function(b) {
          htmltools::tags$li(unname(BUNDLE_LABELS[b]) %||% b)
        })
      )
    )
  }
  shiny::modalDialog(
    title = sprintf("Replace the %s layer?", tag),
    losing,
    htmltools::tags$p(
      class = "muted",
      style = "font-size:12px",
      "The uploaded file differs from the one currently loaded. Keeping ",
      "results computed on the previous data would misreport them as ",
      "belonging to the new data."
    ),
    easyClose = FALSE,
    footer = htmltools::tagList(
      shiny::modalButton("Cancel"),
      shiny::actionButton(ns("confirm_replace"), "Replace and clear",
                          class = "btn btn-danger")
    )
  )
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
      # What scale the numbers are on is not recoverable from the file, and
      # every analysis backend reads it without re-deriving anything. The
      # guess is filled in from the data; this is where it gets corrected.
      shiny::uiOutput(ns("assay_type_picker")),
      shiny::uiOutput(ns("scale_notice")),
      shiny::uiOutput(ns("normalize_controls")),
      shiny::uiOutput(ns("upload_status"))
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
        shiny::actionButton(
          ns("confirm"),
          "Confirm & build omics_input",
          class = "btn btn-primary"
        )
      )
    )
  )
}

# ---- confirmation card -----------------------------------------------

# The schema card answers "what did the classifier decide". This one answers
# "what does that decision mean for my data", which is the question a user can
# actually check. A sheet assignment can be wrong at high confidence and still
# yield an omics_input that analyses cleanly -- metadata read as the matrix
# gives numbers, dimensions, and a full result table, all meaningless. Nothing
# downstream errors on it, so this is the last point where it is catchable.
import_confirm_card <- function(ns) {
  bslib::card(
    bslib::card_header(
      htmltools::tags$h3(class = "card-title", "What will be imported"),
      htmltools::tags$span(
        class = "card-sub",
        "check this against what you know about the file"
      )
    ),
    bslib::card_body(
      shiny::uiOutput(ns("confirm_shape")),
      shiny::uiOutput(ns("confirm_roles")),
      htmltools::tags$div(
        class = "row-grid r-6-6",
        htmltools::tags$div(
          htmltools::tags$h5("Expression matrix"),
          htmltools::tags$div(class = "muted",
                              style = "font-size:12px;margin-bottom:6px",
                              "first rows and columns, as parsed"),
          shiny::tableOutput(ns("confirm_matrix_preview"))
        ),
        htmltools::tags$div(
          htmltools::tags$h5("Sample metadata"),
          htmltools::tags$div(class = "muted",
                              style = "font-size:12px;margin-bottom:6px",
                              "columns available for grouping and covariates"),
          shiny::tableOutput(ns("confirm_meta_preview"))
        )
      )
    )
  )
}

# Numbers first, because "3 features x 240 samples" on a file the user knows
# has 240 features is the fastest way to catch a transposed matrix.
confirm_shape_ui <- function(input_obj, report) {
  if (is.null(input_obj)) {
    return(notice(
      title  = "Nothing to import yet",
      detail = "Upload a file, or correct the sheet roles above if the classifier could not find an expression matrix.",
      kind   = "warn"
    ))
  }
  mat <- input_obj$expr_mat
  n_missing <- sum(is.na(mat))
  orientation <- report$suggested_input$orientation %||% "features_in_rows"

  htmltools::tags$div(
    class = "stat-grid",
    style = "margin-bottom:16px",
    stat_card(
      label = "Features", value = format(nrow(mat), big.mark = ","),
      trend = "rows of the matrix", mono = TRUE
    ),
    stat_card(
      label = "Samples", value = format(ncol(mat), big.mark = ","),
      trend = "columns of the matrix", mono = TRUE
    ),
    stat_card(
      label = "Missing", value = sprintf("%.1f%%", 100 * n_missing / length(mat)),
      trend = sprintf("%s cells", format(n_missing, big.mark = ",")),
      accent = if (n_missing / length(mat) > 0.5) "warn" else "ok"
    ),
    stat_card(
      label = "Orientation",
      value = if (identical(orientation, "features_in_rows")) "features in rows"
              else "samples in rows",
      trend = "swap the role below if reversed"
    )
  )
}

# Which sheet became what, with a dropdown to say otherwise. Sheets the
# classifier could not place are worth showing too: an "unknown" sheet is
# often the metadata, and silently dropping it is how a grouping column goes
# missing later.
# Input id for one sheet's role dropdown. The generation is part of the id so
# a new upload gets fresh controls rather than inheriting the last file's.
role_input_id <- function(gen, i) paste0("role_", gen, "_", i)

confirm_roles_ui <- function(ns, report, gen = 0L) {
  sheets <- report$sheets
  if (is.null(sheets) || nrow(sheets) == 0L) return(NULL)

  choices <- c("expression matrix" = "matrix",
               "sample metadata" = "metadata",
               "feature annotation" = "feature_annot",
               "ignore" = "unknown")

  rows <- lapply(seq_len(nrow(sheets)), function(i) {
    nm <- sheets$name[i]
    low_conf <- !is.na(sheets$confidence[i]) && sheets$confidence[i] < 0.5
    htmltools::tags$div(
      style = "display:flex;align-items:center;gap:10px;margin-bottom:6px",
      htmltools::tags$code(style = "min-width:150px", nm),
      htmltools::tags$span(
        class = "muted", style = "font-size:12px;min-width:110px",
        sprintf("%s x %s", sheets$n_rows[i], sheets$n_cols[i])
      ),
      shiny::selectInput(
        ns(role_input_id(gen, i)), label = NULL, choices = choices,
        selected = sheets$role[i], width = "180px"
      ),
      if (low_conf) {
        pill("low confidence", kind = "warn")
      } else if (identical(sheets$notes[i], "role set by user")) {
        pill("you set this", kind = "ok")
      } else NULL
    )
  })

  htmltools::tagList(
    htmltools::tags$h5("Sheet roles"),
    htmltools::tags$div(class = "muted",
                        style = "font-size:12px;margin-bottom:8px",
                        "Changing a role re-reads the file with your assignment."),
    rows,
    htmltools::tags$hr(style = "margin:14px 0")
  )
}

# A corner of the matrix. Seeing the actual values is what catches a header
# row parsed as data, or an ID column read as a sample.
preview_matrix <- function(mat, n_row = 5L, n_col = 4L) {
  if (is.null(mat) || nrow(mat) == 0L) return(NULL)
  sub <- mat[seq_len(min(n_row, nrow(mat))),
             seq_len(min(n_col, ncol(mat))), drop = FALSE]
  df <- as.data.frame(round(sub, 2))
  df <- cbind(feature = rownames(sub), df)
  rownames(df) <- NULL
  df
}

preview_metadata <- function(meta, n_row = 5L) {
  if (is.null(meta) || nrow(meta) == 0L) return(NULL)
  df <- meta[seq_len(min(n_row, nrow(meta))), , drop = FALSE]
  cbind(sample = rownames(df), as.data.frame(df, stringsAsFactors = FALSE))
}
