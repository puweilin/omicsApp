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

    # ---- confirmation card --------------------------------------------
    output$confirm_shape <- shiny::renderUI({
      cand <- parsed()
      if (is.null(cand)) return(NULL)
      htmltools::tagList(
        confirm_shape_ui(cand$input, cand$report),
        orientation_picker_ui(ns, cand$report)
      )
    })

    # Samples as columns or rows. Asked outright rather than left to the
    # role dropdowns: a transposed matrix imports cleanly and analyses
    # features as samples, and the only sign is a count the user has to
    # notice.
    shiny::observeEvent(input$orientation_pick, {
      cand <- parsed()
      shiny::req(cand)
      current <- cand$report$suggested_input$orientation
      if (identical(input$orientation_pick, current)) return()
      orientation_override(input$orientation_pick)
      do_parse()
    }, ignoreInit = TRUE)

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

    # Which column the gene symbol came from, or that none was found.
    #
    # The import matches it by heading, and a heading can be wrong or
    # absent. Getting it wrong relabels every feature; getting nothing
    # leaves feature_symbol as the accession, and enrichment then
    # returns empty rather than failing -- clusterProfiler answers "No
    # gene can be mapped" and hands back NULL, which arrives as a result
    # with nothing in it. Either way the user cannot tell from the
    # outcome, so it is said here, before anything is run.
    output$confirm_symbol_source <- shiny::renderUI({
      cand <- parsed()
      if (is.null(cand) || is.null(cand$input)) return(NULL)
      src <- attr(cand$input$feature_df, "symbol_column") %||% NA_character_
      syms <- cand$input$feature_df$feature_symbol
      share <- if (length(syms)) mean(is_gene_symbol(syms)) else 0
      # Judged on the values, not on whether a column was named like a
      # gene column: ids that are gene symbols, or Ensembl ids mapped to
      # symbols, were told enrichment "will return nothing".
      if (is.na(src) && share >= 0.5) {
        example <- utils::head(syms[is_gene_symbol(syms)], 3L)
        return(htmltools::tags$div(
          class = "muted", style = "font-size:12.5px;margin:4px 0 12px",
          htmltools::tags$strong("Gene symbol"), htmltools::HTML(" &middot; "),
          sprintf("%.0f%% of features carry a gene symbol ", 100 * share),
          htmltools::tags$span(class = "text-mono", paste(example, collapse = ", "))))
      }

      if (is.na(src)) {
        return(htmltools::tags$div(
          class = "muted",
          style = "font-size:12.5px;margin:4px 0 12px",
          htmltools::tags$strong("Gene symbol"),
          htmltools::HTML(" &middot; "),
          "no gene column found \u2014 feature IDs will be used instead. ",
          htmltools::tags$span(
            style = "color:var(--warn)",
            "Enrichment needs gene symbols, so it will return nothing."
          )
        ))
      }
      example <- utils::head(
        stats::na.omit(cand$input$feature_df$feature_symbol), 3L)
      htmltools::tags$div(
        class = "muted",
        style = "font-size:12.5px;margin:4px 0 12px",
        htmltools::tags$strong("Gene symbol"),
        htmltools::HTML(" &middot; "),
        "taken from ",
        htmltools::tags$span(class = "text-mono", sprintf("\u201c%s\u201d", src)),
        if (length(example)) {
          htmltools::tagList(
            htmltools::HTML(" &middot; "),
            htmltools::tags$span(class = "text-mono",
                                 paste(example, collapse = ", "))
          )
        }
      )
    })

    # Integer values with gene ids imported as proteomics: almost always
    # RNA-seq counts on the wrong radio button, which then go through vsn
    # and limma as intensities.
    output$omics_hint <- shiny::renderUI({
      cand <- parsed()
      if (is.null(cand) || is.null(cand$input)) return(NULL)
      if (!identical(input$omics_type %||% "proteomics", "proteomics")) return(NULL)
      m <- cand$input$expr_mat
      v <- m[!is.na(m)]
      v <- v[seq_len(min(length(v), 20000L))]
      ids <- utils::head(rownames(m), 200L)
      if (length(v) && all(v >= 0) && all(v == round(v)) && max(v) > 100 &&
          mean(grepl("^ENS[A-Z]*G[0-9]+", ids)) > 0.5) {
        notice(title = "This looks like RNA-seq counts",
               detail = "Whole-number values on Ensembl gene ids. If these are read counts, choose RNA-seq above so they are modelled as counts.",
               kind = "warn")
      }
    })

    output$confirm_matrix_preview <- shiny::renderTable({
      cand <- parsed()
      shiny::req(cand, cand$input)
      preview_matrix(cand$input$expr_mat)
    }, striped = TRUE, spacing = "xs", width = "100%", digits = 2)

    # ---- study design ---------------------------------------------------
    output$confirm_design <- shiny::renderUI({
      cand <- parsed()
      shiny::req(cand, cand$input)
      meta <- cand$input$meta_df
      cands <- grouping_candidates(meta)
      if (!length(cands)) {
        info_cols <- setdiff(names(meta), "sample_id")
        msg <- if (!length(info_cols)) {
          paste("This file has no sample information, so there are no groups to",
                "compare. Add a sample sheet above (one row per sample, a column",
                "naming its group), or put it in a second sheet of the workbook.")
        } else {
          counts <- lapply(meta[info_cols], function(x) table(as.character(x)))
          single <- unique(unlist(lapply(counts, function(t) names(t)[t < 2L])))
          paste0("No column splits the samples into groups of two or more",
                 if (length(single)) sprintf(" (only one sample in: %s)",
                                             paste(utils::head(single, 5L), collapse = ", ")),
                 ". A group needs at least two samples to be compared.")
        }
        return(htmltools::tags$div(
          class = "notice notice-warn", style = "font-size:12.5px;margin:6px 0 10px",
          msg))
      }
      sel <- shiny::isolate(input$design_group)
      if (is.null(sel) || !sel %in% c(cands, "")) sel <- cands[[1L]]
      htmltools::tags$div(
        class = "design-picker",
        htmltools::tags$h5("Study design",
                           info_tip(paste("Which column holds the groups, and which group",
                                          "is the control. QC colours by it and the",
                                          "Differential view compares against the control",
                                          "by default; both can still be changed there."))),
        htmltools::tags$div(
          class = "row-grid r-6-6",
          shiny::selectInput(ns("design_group"), "Group column",
                             choices = c(cands, "(none)" = ""), selected = sel),
          shiny::uiOutput(ns("design_reference_ui"))
        )
      )
    })

    output$design_reference_ui <- shiny::renderUI({
      cand <- parsed()
      gc <- input$design_group
      shiny::req(cand, cand$input, gc, nzchar(gc), gc %in% names(cand$input$meta_df))
      lv <- sort(unique(as.character(stats::na.omit(cand$input$meta_df[[gc]]))))
      shiny::selectInput(ns("design_reference"), "Control (reference) group",
                         choices = lv, selected = default_control_level(lv))
    })

    # The design as chosen, checked against the metadata it names.
    apply_design <- function(inp) {
      gc <- input$design_group
      if (is.null(gc) || !nzchar(gc) || !gc %in% names(inp$meta_df)) return(inp)
      ref <- input$design_reference
      lv <- unique(as.character(stats::na.omit(inp$meta_df[[gc]])))
      if (is.null(ref) || !ref %in% lv) ref <- default_control_level(sort(lv))
      tryCatch(omicsCore::set_study_design(inp, gc, ref), error = function(e) inp)
    }

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
            choices = NORMALIZE_CHOICES,
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
      tag <- target_tag(cand)
      if (!tag %in% names(proj$experiments)) return(NULL)
      proj$experiments[[tag]]
    }

    # The name the layer goes in under. A project can hold several layers
    # of one omics type (two proteomics batches), so the name is the
    # user's: it defaults to the omics type, and a name already in the
    # project means "replace that layer".
    target_tag <- function(cand) {
      clean_layer_name(input$layer_name, cand$omics_type %||% "experiment")
    }
    output$layer_name_ui <- shiny::renderUI({
      shiny::req(parse_ok())
      default <- parsed()$input$omics_type %||% "experiment"
      htmltools::tags$div(
        class = "inline-control", style = "max-width:200px",
        shiny::textInput(ns("layer_name"), label = "Layer name",
                         value = shiny::isolate(input$layer_name) %||% default,
                         placeholder = default)
      )
    })
    # The omics type is the default name, and changes with the type
    # control until the user writes a name of their own.
    shiny::observeEvent(parsed()$input$omics_type, {
      typ <- parsed()$input$omics_type
      cur <- input$layer_name
      if (is.null(cur) || !nzchar(trimws(cur)) || cur %in% c("proteomics", "rnaseq")) {
        shiny::updateTextInput(session, "layer_name", value = typ)
      }
    }, ignoreInit = TRUE)

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
      # A different sample sheet is different data.
      sf <- input$sample_file
      if (!is.null(sf)) {
        cand$source_fingerprint <- paste0(cand$source_fingerprint, ":sheet=",
                                          unname(tools::md5sum(sf$datapath)))
      }
      cand
    }

    commit <- function(cand, tag = target_tag(cand)) {
      method <- pending_normalize_method()
      do_normalize <- !identical(method, "none")
      cand <- stamp_fingerprint(cand)
      f <- input$file

      if (do_normalize) {
        normalized <- tryCatch(
          suppressMessages(do.call(omicsCore::normalize_omics,
                                   c(list(cand), normalize_args(method)))),
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
                    names(NORMALIZE_CHOICES)[NORMALIZE_CHOICES == method],
                    cand$assay_type),
            type = "message", duration = 6
          )
        }
      }

      cand <- apply_design(cand)

      # Quantification files are archived one by one, with the names
      # they were uploaded under: export_script() reads them back with
      # read_quant_files(), and the names are the sample names (the
      # archived copies are renamed to stay unique).
      is_quant <- !is.null(parsed()$report$suggested_input$quant_format)
      if (!is.null(f) && is_quant) {
        res <- lapply(seq_along(f$datapath), function(i) {
          store_raw_upload(f$datapath[[i]], f$name[[i]], unname(tools::md5sum(f$datapath[[i]])))
        })
        if (all(vapply(res, function(r) isTRUE(r$ok), logical(1)))) {
          cand$quant_source <- list(paths = vapply(res, `[[`, character(1), "path"),
                                    names = as.character(f$name))
        } else if (any(grepl("quota", vapply(res, function(r) r$message %||% "", character(1)),
                             fixed = TRUE))) {
          shiny::showNotification(res[[1L]]$message, type = "warning", duration = 8)
        }
      } else if (!is.null(f) && length(f$datapath) == 1L) {
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
      sf <- input$sample_file
      if (!is.null(sf) && !is.null(cand$sample_sheet_path)) {
        res <- store_raw_upload(sf$datapath, sf$name,
                                unname(tools::md5sum(sf$datapath)))
        cand$sample_sheet_path <- if (isTRUE(res$ok)) res$path else NULL
      }
      cand$layer_tag <- tag
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
      tag <- target_tag(cand)
      shiny::showModal(
        replace_layer_modal(ns, tag, current_project(),
                            keep_both = next_free_tag(tag, names(current_project()$experiments)))
      )
    })

    shiny::observeEvent(input$confirm_replace, {
      shiny::removeModal()
      shiny::req(parse_ok())
      commit(parsed()$input)
    })
    # Both: the new data as a layer of its own beside the old one, whose
    # results stay.
    shiny::observeEvent(input$confirm_keep_both, {
      shiny::removeModal()
      shiny::req(parse_ok())
      cand <- parsed()$input
      tag <- next_free_tag(target_tag(cand), names(current_project()$experiments))
      shiny::updateTextInput(session, "layer_name", value = tag)
      commit(cand, tag = tag)
    })

    # ---- module return ------------------------------------------------
    shiny::reactive(confirmed_input())
  })
}

# ---- internal helpers ------------------------------------------------

# Tiny `%||%` so the module doesn't pull rlang in just for one operator.
`%||%` <- function(a, b) if (is.null(a)) b else a

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

# Human-readable names for the bundles about to be discarded, so the
# dialog names what is at stake rather than saying "analyses".
BUNDLE_LABELS <- c(
  qc          = "Quality control",
  diff        = "Differential analysis",
  enrich      = "Pathway enrichment",
  integration = "Multi-omics integration"
)

# A layer name as typed, made safe to file under: letters, digits and
# `_ - .` (a space becomes `_`), at most 40 characters, the omics type
# when nothing usable is left.
clean_layer_name <- function(x, default) {
  x <- trimws(as.character(x %||% ""))
  if (!length(x) || !nzchar(x[[1L]])) return(default)
  x <- gsub("[^[:alnum:]_.-]+", "_", x[[1L]], perl = TRUE)
  x <- gsub("^_+|_+$", "", substr(x, 1L, 40L))
  if (nzchar(x)) x else default
}

# `tag`, or `tag_2`, `tag_3`, ... -- the first name no layer has.
next_free_tag <- function(tag, taken) {
  if (!tag %in% taken) return(tag)
  i <- 2L
  while (paste0(tag, "_", i) %in% taken) i <- i + 1L
  paste0(tag, "_", i)
}

replace_layer_modal <- function(ns, tag, project, keep_both = NULL) {
  # The results on this layer, not every result in the project: another
  # layer's analyses are not touched by replacing this one.
  all_b <- project$bundles %||% list()
  bundles <- names(all_b)[vapply(all_b, function(b) {
    !omicsCore::is_analysis_bundle(b) ||
      tag %in% (b$params$experiments %||% character(0)) ||
      identical(omicsCore::bundle_layer(project, b), tag)
  }, logical(1))]
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
      if (!is.null(keep_both)) {
        shiny::actionButton(ns("confirm_keep_both"),
                            sprintf("Keep both (add as '%s')", keep_both),
                            class = "btn btn-outline-primary")
      },
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
      shiny::uiOutput(ns("confirm_symbol_source")),
      # The groups and the control, stated once here by the person who
      # knows the study, and read by QC, Differential and Integration.
      shiny::uiOutput(ns("confirm_design")),
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
      trend = if (identical(report$suggested_input$orientation_source, "user")) "as you set it"
              else sprintf("detected (confidence %.2f)",
                           report$suggested_input$orientation_confidence %||% NA_real_),
      accent = if ((report$suggested_input$orientation_confidence %||% 1) < 0.6) "warn" else "ok"
    )
  )
}

# The control that corrects the orientation, prominent when it was a guess.
orientation_picker_ui <- function(ns, report) {
  sug <- report$suggested_input
  if (is.null(sug$orientation) || identical(sug$matrix_sheet, "rds")) return(NULL)
  guessed <- (sug$orientation_confidence %||% 1) < 0.6 &&
    !identical(sug$orientation_source, "user")
  htmltools::tags$div(
    style = "margin:-6px 0 14px",
    if (guessed) {
      notice(title = "Check the orientation",
             detail = paste("The file did not make it clear whether samples are",
                            "columns or rows. If the Samples count above is really",
                            "the number of features, switch it here."),
             kind = "warn")
    },
    shiny::radioButtons(
      ns("orientation_pick"), label = "Samples are the matrix's",
      choices = c("columns (features in rows)" = "features_in_rows",
                  "rows (samples in rows, e.g. Olink NPX)" = "samples_in_rows"),
      selected = sug$orientation, inline = TRUE)
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
        ns(role_input_id(gen, i)),
        label = htmltools::tags$span(class = "visually-hidden-label", sprintf("Role of sheet %s", nm)),
        choices = choices,
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

# Bound by detached_call() in do_parse(), not visible to R CMD check.
utils::globalVariables(c("datapath", "name", "orientation", "roles", "sample_sheet"))
