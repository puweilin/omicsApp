#' Import view module
#'
#' Phase 3 slice 3A: real upload + smart-parse. The view accepts a
#' single Excel / CSV / TSV / RDS file via `shiny::fileInput()`,
#' hands the path to [omicsCore::read_omics()] (or several per-sample
#' Salmon / RSEM / kallisto files to [omicsCore::read_quant_files()],
#' merged into one layer), renders the resulting
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
                               current_project = shiny::reactiveVal(NULL),
                               navigate = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    if (is.function(navigate)) {
      shiny::observeEvent(input$go_qc, navigate("qc"))
    }

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
    # The matrix orientation the user chose, when the guess was wrong.
    orientation_override <- shiny::reactiveVal(NULL)

    # Bumped on every new file. Shiny keeps an input's value across re-renders
    # of the control, so without this the role dropdowns would still hold the
    # previous workbook's answers and the observer below would apply them to
    # the new one the moment it parsed. Folding the generation into the input
    # ids means a new file starts with genuinely empty controls.
    parse_gen <- shiny::reactiveVal(0L)

    # A new file picked while the previous one is still being read wins;
    # the older parse is dropped when it lands (see run_epoch()).
    parse_epoch <- run_epoch()

    # ---- derived state for the UI -------------------------------------
    has_file       <- shiny::reactive(!is.null(input$file))
    parse_ok       <- shiny::reactive(!is.null(parsed()) &&
                                        !is.null(parsed()$input))
    is_confirmed   <- shiny::reactive(!is.null(confirmed_input()))

    do_parse <- function() {
      f <- input$file
      shiny::req(f)
      omics_type <- input$omics_type %||% "proteomics"
      my_parse <- parse_epoch$start()
      # Off the main thread: a large workbook takes tens of seconds to
      # read, and read in the observer it froze every view of the session
      # -- and, under one R process per container, the whole app -- for
      # that long.
      parsed(NULL)
      confirmed_input(NULL)
      run_async(
        detached_call(
          function() {
            # Parse first with the modality default, then re-label from
            # the values once there is a matrix to look at. read_omics()
            # needs *an* assay_type, and the data it would be inferred
            # from does not exist until it returns.
            assay_type <- if (omics_type == "rnaseq") "raw_count" else "raw_intensity"
            # Several files are one sample each -- Salmon, RSEM or
            # kallisto output -- and are merged into one layer. Named by
            # the names the browser sent, not by their temporary copies.
            several <- length(datapath) > 1L
            read_quant <- function() {
              omicsCore::read_quant_files(datapath, file_names = name,
                                          sample_sheet = sample_sheet)
            }
            out <- tryCatch(
              # The scale check inside omics_input() is muffled for this
              # one call and nothing else: the label handed in here is
              # the modality's default, and it is replaced from the
              # values a few lines down. Warning about a label that is
              # about to be corrected only trains people to ignore the
              # warning that matters, the one at confirm time.
              withCallingHandlers(
                if (several) {
                  read_quant()
                } else {
                  res <- omicsCore::read_omics(
                    datapath,
                    omics_type = omics_type,
                    assay_type = assay_type,
                    sheet_roles = roles,
                    orientation = orientation,
                    sample_sheet = sample_sheet
                  )
                  # One quantification file: read again so its sample is
                  # named after the uploaded file, not "0".
                  if (!is.null(res$report$suggested_input$quant_format)) {
                    single <- read_quant()
                    single$report$warnings <- unique(c(
                      single$report$warnings,
                      grep("read as RNA-seq, not", res$report$warnings, value = TRUE)))
                    res <- single
                  }
                  res
                },
                warning = function(w) {
                  if (grepl("implies (linear|log-scale) values", conditionMessage(w))) {
                    invokeRestart("muffleWarning")
                  }
                }
              ),
              error = function(e) {
                list(
                  input = NULL,
                  report = omicsCore::new_import_report(
                    warnings = paste0(
                      if (several) paste("Several files can be imported together only",
                                         "when each is one sample's Salmon, RSEM or",
                                         "kallisto quantification: ")
                      else "The file could not be read: ",
                      conditionMessage(e)),
                    source = paste(name, collapse = ", ")
                  )
                )
              }
            )
            if (!is.null(out$input)) {
              if (!is.null(out$report$suggested_input$quant_format)) {
                # Estimated reads with their lengths, whatever their totals.
                out$assay_guess <- list(
                  assay_type = "raw_count",
                  reason = "Estimated read counts from the quantification files.")
              } else {
                # The layer's own modality, which is the radio's except for
                # a SummarizedExperiment's read counts: RNA-seq whatever
                # the radio said.
                guess <- omicsCore::infer_assay_type(out$input$expr_mat,
                                                     out$input$omics_type,
                                                     explain = TRUE)
                if (!is.na(guess$assay_type)) {
                  out$input$assay_type <- guess$assay_type
                  out$assay_guess <- guess
                }
              }
            }
            out
          },
          datapath = f$datapath, name = f$name, omics_type = omics_type,
          roles = role_overrides(), orientation = orientation_override(),
          sample_sheet = input$sample_file$datapath
        ),
        on_success = function(out) {
          if (!parse_epoch$is_current(my_parse)) return(invisible())
          # Stamp the user-visible source name so the schema card shows
          # the original filename, not the tempfile path Shiny gave us.
          out$report$source <- paste(f$name, collapse = ", ")
          # Quantification files are RNA-seq whatever the radio said; the
          # radio follows, so the scale choices and the layer name agree.
          # So are a SummarizedExperiment's read counts.
          forced_rnaseq <- !is.null(out$report$suggested_input$quant_format) ||
            (!is.null(out$report$suggested_input$se_class) &&
               identical(out$input$omics_type, "rnaseq"))
          if (forced_rnaseq && !identical(omics_type, "rnaseq")) {
            shiny::updateRadioButtons(session, "omics_type", selected = "rnaseq")
          }
          if (!is.null(out$input)) {
            # Fingerprint the upload so Confirm can tell a genuinely new
            # dataset from the same file picked twice.
            out$input$source_fingerprint <-
              input_fingerprint(f$datapath, omics_type, out$input$assay_type)
          }
          parsed(out)
          # New upload (or radio change) always resets the confirmed
          # state so the user has to re-confirm against the rebuilt input.
          confirmed_input(NULL)
        },
        on_error = function(msg) {
          if (!parse_epoch$is_current(my_parse)) return(invisible())
          parsed(list(input = NULL, report = omicsCore::new_import_report(
            warnings = paste0("The file could not be read: ", msg),
            source = paste(f$name, collapse = ", "))))
        },
        message = if (length(f$name) > 1L) sprintf("Reading %d files...", length(f$name))
                  else sprintf("Reading %s...", f$name)
      )
    }

    # Both templates carry the same donors under different sample ids --
    # the pairing Integration needs, shown rather than described.
    output$template_proteomics <- shiny::downloadHandler(
      filename = function() "omicsapp_template_proteomics.xlsx",
      content = function(file) write_import_template(file, "proteomics")
    )
    output$template_rnaseq <- shiny::downloadHandler(
      filename = function() "omicsapp_template_rnaseq.xlsx",
      content = function(file) write_import_template(file, "rnaseq")
    )

    shiny::observeEvent(input$file, {
      # A new file makes the previous sheet assignment meaningless
      role_overrides(NULL)
      orientation_override(NULL)
      parse_gen(shiny::isolate(parse_gen()) + 1L)
      do_parse()
    })
    shiny::observeEvent(input$omics_type, do_parse(), ignoreInit = TRUE)
    # A sample sheet for a matrix file that has none (counts.txt beside
    # samples.csv): the two are read together.
    shiny::observeEvent(input$sample_file, {
      if (!is.null(input$file)) do_parse()
    }, ignoreInit = TRUE)

    # The value-scale picker and its warning (mod_import_upload.R).
    import_scale_server(input, output, session, parsed, confirmed_input, parse_ok)

    # The "What will be imported" card: shape, sheet roles, gene symbols,
    # study design and normalization (mod_import_confirm.R).
    confirm <- import_confirm_server(input, output, session, parsed, parse_gen,
                                     orientation_override, role_overrides, do_parse,
                                     parse_ok)
    apply_design <- confirm$apply_design
    normalizable <- confirm$normalizable

    # The steps strip, the upload card and the schema card
    # (mod_import_upload.R).
    import_status_server(input, output, session, parsed, has_file, parse_ok, is_confirmed)

    # ---- confirm button gating ----------------------------------------
    # Disabled until there is something to import. It used to be live
    # from the start and do nothing when pressed before a parse had
    # succeeded -- req() stopped it silently -- which read as broken.
    shiny::observe({
      ok <- isTRUE(parse_ok())
      tryCatch(shinyjs::toggleState("confirm", condition = ok),
               error = function(e) NULL)
    })

    output$confirm_state <- shiny::renderUI({
      if (is_confirmed()) {
        return(htmltools::tags$div(
          style = "display:flex;align-items:center;gap:8px;justify-content:flex-end",
          pill("layer imported", kind = "ok"),
          htmltools::tags$span(
            class = "muted", style = "font-size:12px",
            sprintf("%d features \u00D7 %d samples",
                    nrow(confirmed_input()$expr_mat),
                    ncol(confirmed_input()$expr_mat))
          ),
          # The import is done; say what comes next rather than leaving
          # the user to find it in the sidebar.
          if (is.function(navigate)) {
            shiny::actionButton(ns("go_qc"), "Next: Quality control \u2192",
                                class = "btn btn-sm btn-primary")
          }
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

    # Confirm: name the layer, normalize, archive the upload, and ask
    # before replacing live data (mod_import_commit.R).
    import_commit_server(input, output, session, current_project, parsed, parse_ok,
                         normalizable, apply_design, confirmed_input)

    # ---- module return ------------------------------------------------
    shiny::reactive(confirmed_input())
  })
}

# ---- internal helpers ------------------------------------------------

# The normalisations offered on import. "log2, samples aligned on their
# medians" is log2 followed by subtracting each sample's median and adding
# back the overall median: it takes out a sample that was simply loaded
# with more material, and keeps the values on a log2-intensity scale.
NORMALIZE_CHOICES <- c(
  "vsn (variance stabilising)" = "vsn",
  "log2" = "log2",
  "log2, samples aligned on their medians" = "log2_median"
)

# A choice above as the arguments normalize_omics() takes.
normalize_args <- function(choice) {
  switch(choice,
         log2_median = list(method = "log2", center = "median"),
         list(method = choice))
}

# Identity of an upload. The file digest alone is not enough: the same
# workbook imported as proteomics and as RNA-seq yields two different
# inputs, so the parse settings are part of what makes it "the same".
# `normalize` is in here for the same reason -- re-importing the same file
# with normalization turned off produces different numbers, and that has to
# read as new data rather than as "nothing changed".
input_fingerprint <- function(path, omics_type, assay_type, normalize = NULL) {
  digest <- unname(tools::md5sum(path))
  if (!length(digest) || anyNA(digest)) return(NULL)
  # Several quantification files make one layer, and all of them are its
  # identity.
  digest <- paste(digest, collapse = "+")
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

# Bound by detached_call() in do_parse(), not visible to R CMD check.
utils::globalVariables(c("datapath", "name", "orientation", "roles", "sample_sheet"))
