#' Differential analysis view module
#'
#' Slice 3D: replaces the inert mockup controls with a live design
#' panel. The Method dropdown lists every backend supported by
#' `omicsCore::run_diff()`; engines that need an absent
#' Bioconductor Suggest are kept in the list but disabled via a
#' notice strip. Group column / Control / Case / Covariates are
#' populated from the active experiment's `meta_df`. The Re-run
#' button is gated by `bindEvent`; failures surface in a notice
#' strip instead of Shiny's red overlay.
#'
#' With several treatment groups the Case control takes more than one
#' level: every chosen group is fitted against the Control in one model
#' (see [omicsCore::run_diff()]), a "Showing" selector picks which
#' contrast the volcano, the hit table and everything downstream read,
#' and a summary card compares the contrasts side by side.
#'
#' Nothing runs until the user presses Run.
#'
#' Reference markup: `omicsApp/mockup/index.html:753-915`.
#'
#' @param id Module namespace id.
#' @keywords internal
#' @noRd
diff_view_ui <- function(id) {
  ns <- shiny::NS(id)

  htmltools::tagList(
    shiny::uiOutput(ns("header")),
    shiny::uiOutput(ns("notices")),
    shiny::uiOutput(ns("stats")),
    htmltools::tags$div(
      class = "row-grid r-3-9",
      diff_params_card(ns),
      htmltools::tags$div(
        shiny::uiOutput(ns("contrast_summary")),
        diff_volcano_card(ns),
        diff_hits_card(ns),
        shiny::uiOutput(ns("anova_card"))
      )
    )
  )
}

#' @rdname diff_view_ui
#' @param current_project Reactive (or reactiveVal) yielding the
#'   live `omics_project` or `NULL`.
#' @keywords internal
#' @noRd
diff_view_server <- function(id, current_project = shiny::reactiveVal(NULL),
                             invalidate = shiny::reactiveVal(0L),
                             navigate = NULL) {
  shiny::moduleServer(id, function(input, output, session) {

    # Active experiment: whichever layer the user picked, else
    # proteomics, else the first of any kind.
    #
    # The layer is a control rather than a fixed choice because the
    # engines available depend on it: applicable_diff_methods() offers
    # deseq2 and edger only for rnaseq raw counts. Pinned to proteomics,
    # as this was, those two could never be reached from the interface
    # at all -- the gate was right and there was no way to the side of
    # it where it opens.
    #
    # The demo resolves through example_project() for the same reason,
    # rather than being handed a fixed proteomics input: its rnaseq
    # layer is raw counts, so it is the one place a user can see the
    # method list change without importing anything.
    # The layer a restored result was computed on, until the control exists.
    preferred_layer <- shiny::reactiveVal(NULL)
    skip_layer_clear <- FALSE

    active <- shiny::reactive({
      proj <- current_project()
      is_demo <- is.null(proj)
      if (is_demo) proj <- example_project()
      exps <- proj$experiments
      if (length(exps) == 0L) {
        return(list(input = NULL, tag = NULL, is_demo = TRUE))
      }
      want <- input$layer %||% preferred_layer()
      tag <- if (!is.null(want) && want %in% names(exps)) {
        want
      } else {
        default_layer_tag(exps)
      }
      list(input = exps[[tag]], tag = tag, is_demo = is_demo)
    })

    output$ui_layer <- shiny::renderUI({
      proj <- current_project() %||% example_project()
      tags_avail <- names(proj$experiments)
      if (length(tags_avail) == 0L) return(NULL)
      # isolate(), because active() reads input$layer and this output
      # writes it. Reading it here closed the loop: re-rendering the
      # control re-sent its value, which invalidated active(), which
      # re-rendered the control.
      sel <- shiny::isolate(input$layer) %||% shiny::isolate(preferred_layer())
      if (is.null(sel) || !sel %in% tags_avail) {
        sel <- default_layer_tag(proj$experiments)
      }
      shiny::selectInput(
        session$ns("layer"), label = "Omics layer",
        choices = tags_avail, selected = sel
      )
    })

    # Method dropdown, restricted to the engines whose assumptions the
    # active layer meets. DESeq2 handed continuous intensities does not
    # error: it rounds them to integers and reports p-values for a
    # negative-binomial model the data never fitted. Nothing downstream
    # can tell that apart from a real result, so the guard has to be
    # here, at the point of choosing.
    output$ui_method <- shiny::renderUI({
      a <- active()
      inp <- a$input
      shiny::req(inp)   # a project with no layers has no methods to offer
      choices <- omicsCore::applicable_diff_methods(inp)
      # Keep the user's pick. active() changes whenever the project does
      # -- including when a finished run is attached to it -- so the
      # control re-rendered after every run and snapped back to "auto".
      sel <- shiny::isolate(input$method)
      if (is.null(sel) || !sel %in% choices) sel <- "auto"
      shiny::selectInput(session$ns("method"),
                         label = htmltools::tags$span(class = "visually-hidden-label", "Method"),
                         choices = choices, selected = sel)
    })

    output$method_note <- shiny::renderUI({
      a <- active()
      inp <- a$input
      shiny::req(inp)
      dropped <- setdiff(omicsCore::SUPPORTED_DIFF_METHODS,
                         omicsCore::applicable_diff_methods(inp))
      if (length(dropped) == 0L) return(NULL)
      htmltools::tags$div(
        class = "muted", style = "font-size:11.5px;padding-top:4px",
        sprintf("%s hidden: not valid for %s data.",
                paste(dropped, collapse = ", "),
                inp$assay_type %||% inp$omics_type)
      )
    })

    # The contrast the controls will default to, computed from the data
    # rather than read back off them.
    #
    # do_run() fires once as soon as active() settles, and at that
    # moment the controls do not exist yet: they are renderUI output,
    # and Shiny has not been round the loop. Reading input$group_col
    # there gets NULL, which used to surface as "Pick a group column
    # with distinct Control and Case levels" on a view the user had
    # only just opened. Deriving the default in one place means the
    # first run and the rendered controls cannot disagree about it.
    default_contrast <- shiny::reactive({
      meta <- active()$input$meta_df
      if (is.null(meta) || !ncol(meta)) {
        return(list(group_col = NULL, control = NULL, case = NULL,
                    levels = character(0), candidates = character(0)))
      }
      cands <- grouping_candidates(meta)
      # No column with two replicated groups. Offer the columns that at
      # least repeat a value -- never one that names each sample once,
      # which is the sample ID and cannot be compared.
      if (length(cands) == 0L) cands <- grouping_candidates(meta, min_per_level = 1L,
                                                            replicated = TRUE)
      # What was stated at import comes first; the guess is the fallback.
      design <- tryCatch(omicsCore::study_design(active()$input), error = function(e) NULL)
      if (!is.null(design)) cands <- c(design$group_col, setdiff(cands, design$group_col))
      # NULL, not cands[1L]: that is NA on an empty vector, and the
      # "no column splits the samples" message below never fired.
      gc <- if (length(cands)) cands[1L] else NULL
      if (is.null(gc)) {
        return(list(group_col = NULL, control = NULL, case = NULL,
                    levels = character(0), candidates = character(0)))
      }
      lv <- sort(unique(as.character(stats::na.omit(meta[[gc]]))))
      ctrl <- if (!is.null(design$reference) && identical(gc, design$group_col)) {
        design$reference
      } else default_control_level(lv)
      n_per <- table(as.character(meta[[gc]]))
      list(
        group_col  = gc,
        control    = ctrl,
        # Every other group against the control: with one control and
        # several treatments that is the design people mean, and with two
        # groups it is the one contrast there is. A group of one sample
        # cannot be tested and is left out of the default.
        case       = if (length(lv) >= 2L) setdiff(lv[n_per[lv] >= 2L], ctrl) else NULL,
        levels     = lv,
        candidates = cands
      )
    })

    # Groups or a continuous variable (dose, time, age). Numeric columns
    # with many values had no way in: they are not groups, and the view
    # offered nothing else.
    continuous_cols <- shiny::reactive({
      meta <- active()$input$meta_df
      if (is.null(meta)) return(character(0))
      names(meta)[vapply(meta, function(col) {
        v <- suppressWarnings(as.numeric(as.character(col)))
        sum(!is.na(v)) >= 4L && length(unique(stats::na.omit(v))) >= 3L
      }, logical(1))]
    })
    design_mode <- shiny::reactive(input$design_mode %||% "groups")
    # A trend's effect is a slope, not a fold change; the cutoff says so.
    shiny::observeEvent(design_mode(), {
      shiny::updateNumericInput(
        session, "fc_cut",
        label = if (identical(design_mode(), "continuous")) "|slope| cutoff" else "|log2FC| cutoff")
    }, ignoreInit = TRUE)

    output$ui_design_mode <- shiny::renderUI({
      if (!length(continuous_cols())) return(NULL)
      shiny::radioButtons(session$ns("design_mode"), label = "Compare",
                          choices = c("groups" = "groups",
                                      "a continuous variable" = "continuous"),
                          selected = shiny::isolate(input$design_mode) %||% "groups",
                          inline = TRUE)
    })

    output$ui_continuous <- shiny::renderUI({
      if (!identical(design_mode(), "continuous")) return(NULL)
      cols <- continuous_cols()
      sel <- shiny::isolate(input$continuous_col)
      if (is.null(sel) || !sel %in% cols) sel <- cols[[1L]]
      htmltools::tagList(
        shiny::selectInput(session$ns("continuous_col"), "Variable (dose, time, age...)",
                           choices = cols, selected = sel),
        shiny::radioButtons(session$ns("continuous_model"), label = "Model",
                            choices = c("linear trend" = "linear",
                                        "any smooth change (spline)" = "spline"),
                            selected = shiny::isolate(input$continuous_model) %||% "linear",
                            inline = TRUE),
        htmltools::tags$div(
          class = "muted", style = "font-size:11.5px;margin-top:-6px",
          paste("Each feature is tested for a trend with the variable. The",
                "effect shown is Spearman's rho; a spline asks whether it",
                "changes at all, without a direction."))
      )
    })

    # Columns that pair samples across the groups: one sample per group
    # for each patient, donor, mouse. Offered whenever such a column
    # exists, and chosen for the user when its name says it is one --
    # a before/after design analysed unpaired found nothing at all.
    pairing_candidates <- shiny::reactive({
      meta <- active()$input$meta_df
      gc <- input$group_col %||% default_contrast()$group_col
      if (is.null(meta) || is.null(gc) || !gc %in% names(meta)) return(character(0))
      g <- as.character(meta[[gc]])
      names(meta)[vapply(names(meta), function(nm) {
        if (identical(nm, gc)) return(FALSE)
        b <- as.character(meta[[nm]])
        ok <- !is.na(b) & !is.na(g)
        if (sum(ok) < 4L) return(FALSE)
        per_block <- table(b[ok])
        nb <- length(per_block)
        # Repeated, but not one level per sample and not the group again.
        nb >= 2L && nb < sum(ok) && all(per_block >= 2L) &&
          # The same block appears in more than one group.
          mean(tapply(g[ok], b[ok], function(x) length(unique(x))) >= 2L) >= 0.5
      }, logical(1))]
    })

    output$ui_paired <- shiny::renderUI({
      if (identical(design_mode(), "continuous")) return(NULL)
      cands <- pairing_candidates()
      if (!length(cands)) return(NULL)
      sel <- shiny::isolate(input$paired_col)
      if (is.null(sel) || !sel %in% c("", cands)) {
        hinted <- cands[grepl(PAIRING_COL_HINTS, tolower(cands))]
        sel <- if (length(hinted)) hinted[[1L]] else ""
      }
      htmltools::tagList(
        shiny::selectInput(session$ns("paired_col"),
                           label = "Paired by (same patient / donor / animal)",
                           choices = c("not paired" = "", cands), selected = sel),
        htmltools::tags$div(
          class = "muted", style = "font-size:11.5px;margin-top:-6px",
          paste("Each pair is compared with itself, which removes the",
                "differences between individuals."))
      )
    })

    # Group column dropdown: any meta_df column with >= 2 unique
    # non-NA values (continuous columns like `age` are excluded
    # for the simple "control vs case" UI in this slice).
    output$ui_group_col <- shiny::renderUI({
      d <- default_contrast()
      if (!length(d$candidates)) return(NULL)
      # Keep what the user picked. This output re-renders whenever the
      # layer or the project changes, and re-rendering with the default
      # threw away their choice -- they selected condition / G1 / G2,
      # ran it, and the control came back reading `label`.
      sel <- shiny::isolate(input$group_col)
      if (is.null(sel) || !sel %in% d$candidates) sel <- d$group_col
      shiny::selectInput(session$ns("group_col"),
                         label    = "Group column",
                         choices  = d$candidates,
                         selected = sel)
    })

    # Reactive level set for the chosen group column.
    levels_ <- shiny::reactive({
      a <- active()
      meta <- a$input$meta_df
      gc <- input$group_col %||% default_contrast()$group_col
      if (is.null(gc) || !(gc %in% names(meta))) return(character(0))
      sort(unique(as.character(stats::na.omit(meta[[gc]]))))
    })

    output$ui_contrast <- shiny::renderUI({
      if (identical(design_mode(), "continuous")) return(NULL)
      lv <- levels_()
      if (length(lv) < 2L) {
        return(htmltools::tags$div(
          class = "muted",
          style = "font-size:12px",
          "Pick a group column with at least two levels."
        ))
      }
      keep <- function(current, fallback) {
        current <- shiny::isolate(current)
        current <- current[current %in% lv]
        if (!length(current)) fallback else current
      }
      d <- default_contrast()
      meta <- active()$input$meta_df
      gc_now <- input$group_col %||% d$group_col
      n_per <- table(as.character(meta[[gc_now]]))
      small <- lv[n_per[lv] < 3L]
      ctrl <- keep(input$control,
                   if (!is.null(d$control) && d$control %in% lv &&
                       identical(input$group_col %||% d$group_col, d$group_col))
                     d$control else default_control_level(lv))[1L]
      mode <- shiny::isolate(input$contrast_mode) %||% "control"
      htmltools::tagList(
        # Three designs, one fit each: every treatment against the control
        # (the common case), every pair of groups, or comparisons written
        # out -- "(TreatA + TreatB)/2 - Control", "TreatB - TreatA".
        if (length(lv) > 2L) {
          shiny::radioButtons(
            session$ns("contrast_mode"), label = "Compare",
            choices = c("each vs control" = "control",
                        "all pairs" = "pairwise",
                        "custom" = "custom"),
            selected = mode, inline = TRUE)
        },
        shiny::selectInput(session$ns("control"),
                           label = "Control (reference)", choices = lv,
                           selected = ctrl),
        shiny::conditionalPanel(
          sprintf("!input['%s'] || input['%s'] == 'control'",
                  session$ns("contrast_mode"), session$ns("contrast_mode")),
          shiny::selectizeInput(
            session$ns("case"),
            label = if (length(lv) > 2L) "Compare against control (one or more)" else "Case",
            # The control is not something to compare with itself.
            choices = setdiff(lv, ctrl), multiple = TRUE,
            selected = setdiff(keep(input$case, setdiff(lv[n_per[lv] >= 2L], ctrl)), ctrl),
            options = list(plugins = list("remove_button")))
        ),
        if (length(small)) {
          htmltools::tags$div(
            class = "muted", style = "font-size:11.5px;color:var(--warn)",
            sprintf("Few samples: %s. A group needs at least 2 samples to be tested and 3 for a reliable estimate.",
                    paste(sprintf("%s (n = %d)", small, as.integer(n_per[small])), collapse = ", ")))
        },
        if (length(lv) > 2L) shiny::conditionalPanel(
          sprintf("input['%s'] == 'custom'", session$ns("contrast_mode")),
          shiny::textAreaInput(
            session$ns("custom_contrasts"),
            label = "Comparisons, one per line",
            value = shiny::isolate(input$custom_contrasts) %||% "",
            placeholder = paste(
              c(paste(contrast_token(lv[3L]), "-", contrast_token(lv[2L])),
                sprintf("(%s + %s)/2 - %s", contrast_token(lv[2L]),
                        contrast_token(lv[3L]), contrast_token(lv[1L]))),
              collapse = "\n"),
            rows = 3),
          htmltools::tags$div(
            class = "muted", style = "font-size:11.5px;margin-top:-6px",
            "Group names that are not plain words go in backticks: `Drug A` - Control.")
        ),
        if (length(lv) > 2L) {
          htmltools::tags$div(
            class = "muted", style = "font-size:11.5px;margin-top:-6px",
            paste("Every group compared is fitted in one model, so each",
                  "comparison borrows strength from all of their samples."))
        }
      )
    })

    # The control is not a case: keep the case list free of it when the
    # reference changes.
    shiny::observeEvent(input$control, {
      lv <- levels_()
      shiny::updateSelectizeInput(session, "case", choices = setdiff(lv, input$control),
                                  selected = setdiff(input$case, input$control))
    }, ignoreInit = TRUE)

    output$ui_covariates <- shiny::renderUI({
      a <- active()
      meta <- a$input$meta_df
      if (identical(input$method, "ttest")) {
        return(htmltools::tags$div(class = "muted", style = "font-size:12px",
                                   "The t-test does not adjust for covariates; choose limma or lm to add them."))
      }
      gc <- if (identical(design_mode(), "continuous")) input$continuous_col %||% ""
            else input$group_col %||% ""
      cands <- setdiff(names(meta), c(gc, "sample_id", input$paired_col))
      # Kept across re-renders, for the reason given at ui_method.
      sel <- intersect(shiny::isolate(input$covariates), cands)
      shiny::selectizeInput(
        session$ns("covariates"),
        label    = NULL,
        choices  = cands,
        multiple = TRUE,
        selected = if (length(sel)) sel,
        options  = list(placeholder = "optional, e.g. age")
      )
    })

    # Diff bundle: demo fallback when no project; otherwise gated
    # behind the Re-run button. We also auto-run once on first
    # mount when a real project is present, so the user lands on
    # a populated volcano without a Re-run click. After that,
    # changes only take effect on Re-run.
    diff_bundle <- shiny::reactiveVal(NULL)
    diff_error  <- shiny::reactiveVal(NULL)
    running     <- shiny::reactiveVal(FALSE)
    # Every run and every reset takes a new number. A result that comes
    # back under an older number is about a layer or a project that is
    # no longer on screen, and is dropped instead of shown as current.
    diff_epoch  <- run_epoch()
    anova_epoch <- run_epoch()

    # The button says what it will do, and cannot be pressed twice while
    # a run is in flight.
    set_busy <- function(busy) {
      running(busy)
      tryCatch(
        if (busy) shinyjs::disable("rerun") else shinyjs::enable("rerun"),
        error = function(e) NULL)
    }

    # ---- which contrast is on screen ---------------------------------
    comparisons <- shiny::reactive({
      b <- diff_bundle()
      if (!omicsCore::is_analysis_bundle(b)) return(character(0))
      omicsCore::diff_comparisons(b)
    })

    output$ui_comparison <- shiny::renderUI({
      cmp <- comparisons()
      if (length(cmp) < 2L) return(NULL)
      sel <- shiny::isolate(input$comparison)
      if (is.null(sel) || !sel %in% cmp) {
        sel <- diff_bundle()$params$shown_comparison %||% cmp[[1L]]
        if (!sel %in% cmp) sel <- cmp[[1L]]
      }
      shiny::selectInput(session$ns("comparison"), label = "Showing",
                         choices = stats::setNames(cmp, gsub("_vs_", " vs ", cmp)),
                         selected = sel)
    })

    # One contrast, shaped exactly like a single-contrast run. Everything
    # below -- the volcano, the table, Enrichment, Integration, the
    # report -- reads this, so none of them needs to know that several
    # contrasts were fitted.
    shown_bundle <- shiny::reactive({
      b <- diff_bundle()
      if (!omicsCore::is_analysis_bundle(b)) return(b)
      cmp <- comparisons()
      if (length(cmp) < 2L) return(b)
      sel <- input$comparison
      if (is.null(sel) || !sel %in% cmp) sel <- cmp[[1L]]
      omicsCore::select_comparison(b, sel)
    })

    # The layer this result was computed on has been replaced, so the
    # result is no longer about anything in the project. NULL is the
    # module's own start-up state, so this only rewinds it.
    shiny::observeEvent(invalidate(), {
      diff_epoch$bump()
      anova_epoch$bump()
      diff_bundle(NULL)
      diff_error(NULL)
      anova_bundle(NULL)
    }, ignoreInit = TRUE)

    do_run <- function() {
      a <- active()
      if (is.null(a$input)) {
        diff_error("This project has no layers to analyse.")
        return(invisible())
      }
      d <- default_contrast()
      method     <- input$method %||% "auto"
      # Fall back to the derived default, not to nothing: on the first
      # run the controls have not been rendered yet, and refusing then
      # put an error on a view the user had only just opened.
      group_col  <- input$group_col %||% d$group_col
      control    <- input$control   %||% d$control
      case       <- input$case      %||% d$case
      case       <- setdiff(case, control)
      covariates <- input$covariates
      mode       <- input$contrast_mode %||% "control"
      contrasts  <- switch(
        mode,
        pairwise = "pairwise",
        custom   = {
          lines <- trimws(strsplit(input$custom_contrasts %||% "", "\n")[[1L]])
          lines[nzchar(lines)]
        },
        NULL)
      if (identical(design_mode(), "continuous")) {
        return(run_continuous(a, method, covariates))
      }
      paired_col <- input$paired_col
      if (!length(paired_col) || !nzchar(paired_col)) paired_col <- NULL
      if (identical(mode, "custom") && !length(contrasts)) {
        diff_error("Write at least one comparison, e.g. \"TreatB - TreatA\".")
        return(invisible())
      }
      if (is.null(group_col) && !length(d$candidates)) {
        diff_error(paste("No column of the sample metadata splits the samples",
                         "into groups. Add one (e.g. 'group' = Control / Treated)",
                         "to the sample sheet and import again."))
        return(invisible())
      }
      if (is.null(contrasts) &&
          (is.null(group_col) || is.null(control) || !length(case))) {
        diff_error("Pick a group column, a control group, and at least one group distinct from the control to compare with it.")
        return(invisible())
      }
      set_busy(TRUE)
      my_run <- diff_epoch$start()
      run_async(
        # Detached, so the worker receives the input and the six
        # parameters rather than this module's whole scope. Defined
        # inline it carried the previous bundle, the project and the
        # demo fixtures with it -- 527 MB of "globals", which future
        # refused to export, from inside an observer, which greyed the
        # page.
        detached_call(
          function() {
            omicsCore::run_diff(
              input         = inp,
              method        = method,
              analysis_type = "group",
              group_col     = group_col,
              control_group = control,
              case_group    = if (is.null(contrasts)) case,
              covariates    = covariates,
              paired_col    = paired_col,
              contrasts     = contrasts
            )
          },
          inp        = a$input,
          method     = method,
          group_col  = group_col,
          control    = control,
          case       = case,
          covariates = if (length(covariates)) covariates else NULL,
          paired_col = paired_col,
          contrasts  = contrasts
        ),
        on_success = function(bundle) {
          if (diff_epoch$is_last_started(my_run)) set_busy(FALSE)
          if (!diff_epoch$is_current(my_run)) return(invisible())
          diff_error(NULL)
          diff_bundle(bundle)
        },
        on_error = function(msg) {
          if (diff_epoch$is_last_started(my_run)) set_busy(FALSE)
          if (!diff_epoch$is_current(my_run)) return(invisible())
          # The previous result stays, labelled as the previous one: a
          # mistyped contrast should not lose it, and unlabelled beside
          # the error it read as the answer to the new settings.
          diff_error(msg)
        },
        message = "Running differential analysis..."
      )
    }

    run_continuous <- function(a, method, covariates) {
      col <- input$continuous_col %||% continuous_cols()[1L]
      if (is.null(col) || is.na(col)) {
        diff_error("Pick a numeric column to test for a trend.")
        return(invisible())
      }
      model <- input$continuous_model %||% "linear"
      set_busy(TRUE)
      my_run <- diff_epoch$start()
      run_async(
        detached_call(
          function() {
            args <- list(input = inp, method = method, analysis_type = "continuous",
                         continuous_col = col, covariates = covariates)
            # The spline is limma's; a linear trend runs on any engine.
            if (identical(model, "spline")) {
              args$method <- "limma"
              args$model <- "spline"
            }
            do.call(omicsCore::run_diff, args)
          },
          inp = a$input, method = method, col = col, model = model,
          covariates = if (length(covariates)) covariates else NULL
        ),
        on_success = function(bundle) {
          if (diff_epoch$is_last_started(my_run)) set_busy(FALSE)
          if (!diff_epoch$is_current(my_run)) return(invisible())
          diff_error(NULL)
          diff_bundle(bundle)
        },
        on_error = function(msg) {
          if (diff_epoch$is_last_started(my_run)) set_busy(FALSE)
          if (!diff_epoch$is_current(my_run)) return(invisible())
          diff_error(msg)
        },
        message = "Testing for trends..."
      )
    }

    # Run once, when the view first settles, so the user lands on a
    # populated volcano rather than an empty panel.
    #
    # After that a change of layer *clears* the result instead of
    # recomputing it. Re-running would be worse than either extreme:
    # picking a layer is the first of several decisions -- method,
    # contrast, covariates -- and spending a DESeq2 run on the state
    # halfway through them is work nobody asked for, on settings nobody
    # has finished choosing. Keeping the old result would be worse
    # still: it is about the previous layer, and nothing on screen
    # would say so.
    # Nothing runs until Re-run is pressed.
    #
    # The view used to analyse once on arrival so it opened populated.
    # On the demo that is milliseconds; on a real workbook it is limma
    # or DESeq2 over thousands of features, started by walking into the
    # view, on a contrast the user has not looked at yet. Opening a tab
    # is not a request to compute.
    #
    # Changing layer clears the result for the same reason it is not
    # re-run: the result describes the previous layer, and nothing on
    # screen would say so. Keyed on input$layer rather than on active(),
    # which also invalidates when the control below it re-renders.
    shiny::observeEvent(input$layer, {
      # The switch a restore made itself, to the restored result's layer.
      if (isTRUE(skip_layer_clear)) {
        skip_layer_clear <<- FALSE
        return()
      }
      # The control's first value, or a switch back to the layer the
      # result is about: nothing to clear.
      b <- diff_bundle()
      exps <- (current_project() %||% example_project())$experiments
      if (!is.null(b) && identical(exps[[input$layer]]$omics_type %||% NA,
                                   b$input_info$omics_type)) return()
      diff_epoch$bump()
      anova_epoch$bump()
      diff_bundle(NULL)
      diff_error(NULL)
      anova_bundle(NULL)
    }, ignoreInit = TRUE)

    # A project opened or restored brings its results: the view shows the
    # saved comparison (and global test) on the layer it was computed
    # on, instead of "no result yet" beside a Workflow card that says
    # Differential is done. After the generation bump (priority -5), so
    # the clearing it causes does not undo this.
    shiny::observeEvent(current_project(), {
      proj <- current_project()
      b <- proj$bundles$diff
      if (is.null(proj) || is.null(b) || !is.null(diff_bundle())) return()
      types <- vapply(proj$experiments, function(e) e$omics_type %||% "", character(1))
      tag <- names(proj$experiments)[types == (b$input_info$omics_type %||% "")][1L]
      if (is.na(tag)) return()
      preferred_layer(tag)
      if (!is.null(input$layer) && !identical(input$layer, tag)) {
        skip_layer_clear <<- TRUE
        shiny::updateSelectInput(session, "layer", selected = tag)
      }
      diff_epoch$bump()
      diff_error(NULL)
      diff_bundle(b)
      a <- proj$bundles$anova
      if (!is.null(a) && identical(a$input_info$omics_type, b$input_info$omics_type)) {
        anova_bundle(a)
      }
    }, priority = -5, ignoreNULL = TRUE)

    # Re-run button is the user-driven path. bindEvent semantics
    # via observeEvent: any change to the controls *not* gated on
    # rerun is ignored except for refreshing the contrast UI
    # populated above.
    shiny::observeEvent(input$rerun, {
      do_run()
    })

    # Slider-derived significance mask. Shared by stat cards,
    # volcano, and top-hits table; recomputed on slider change
    # without re-running the full diff. The thresholds are
    # debounced so a slider drag fires one mask update instead of
    # one per pixel.
    fdr_cut_d <- shiny::debounce(shiny::reactive(input$fdr_cut %||% 0.05), 250)
    # numericInput hands back NA while the box is empty mid-typing, and
    # NA here would mark every feature non-significant with no
    # explanation. Fall back to the default rather than to nothing.
    fc_cut_d  <- shiny::debounce(shiny::reactive({
      v <- input$fc_cut
      if (is.null(v) || !is.finite(v) || v < 0) round(log2(1.2), 3) else v
    }), 250)
    # Which column "significant" is read from. The label follows it, so
    # a figure never says adj.P over a raw-p mask.
    p_col   <- shiny::reactive(
      if (identical(input$p_kind %||% "adj", "raw")) "p_value" else "adj_p_value")
    p_label <- shiny::reactive(
      if (identical(input$p_kind %||% "adj", "raw")) "p" else "adjusted p")

    marked <- shiny::reactive({
      shiny::req(shown_bundle())
      df <- shown_bundle()$results$diff_result_df
      pv <- df[[p_col()]]
      df$is_significant <- !is.na(pv) &
                           !is.na(df$effect) &
                           pv < fdr_cut_d() &
                           abs(df$effect) >= fc_cut_d()
      df
    })

    output$header <- shiny::renderUI({
      a <- active()
      b <- shown_bundle()
      source_note <- htmltools::tags$span(
        class = "muted",
        if (a$is_demo) "demo project (built-in)"
        else sprintf("layer = %s", a$tag)
      )

      # Before the first result there is nothing to summarise. Three em
      # dashes separated by middots said that in a way that read as
      # damage rather than as an empty state, which is how it was
      # reported.
      if (is.null(b)) {
        return(view_header(
          title    = "Differential",
          subtitle = htmltools::tagList(
            htmltools::tags$span(class = "muted", "no result yet"),
            htmltools::HTML(" &middot; "),
            source_note
          )
        ))
      }

      stale <- !is.null(diff_error()) && !identical(diff_error(), CANCELLED_MESSAGE)
      comparison <- if (length(b$params$comparison)) gsub("_vs_", " vs ", b$params$comparison) else
        sprintf("%s vs %s",
                b$params$case_group %||% input$case %||% "case",
                b$params$control_group %||% input$control %||% "control")
      n_cmp <- length(comparisons())
      view_header(
        title    = "Differential",
        actions  = if (is.function(navigate)) {
          shiny::actionButton(session$ns("go_enrich"), "Next: Enrichment \u2192",
                              class = "btn btn-ghost")
        },
        subtitle = htmltools::tagList(
          diff_omics_display(b$input_info$omics_type),
          if (n_cmp > 1L) htmltools::tagList(
            htmltools::HTML(" &middot; "),
            sprintf("%d comparisons", n_cmp)),
          htmltools::HTML(" &middot; "),
          comparison,
          htmltools::HTML(" &middot; "),
          b$params$method,
          htmltools::HTML(" &middot; "),
          source_note,
          if (stale) htmltools::tagList(
            htmltools::HTML(" &middot; "),
            htmltools::tags$span(style = "color:var(--warn);font-weight:600",
                                 "previous result \u2014 the latest run failed"))
        )
      )
    })

    output$notices <- shiny::renderUI({
      err <- diff_error()
      missing_engines <- diff_missing_engines()
      tagged <- htmltools::tagList()
      if (identical(err, CANCELLED_MESSAGE)) {
        tagged <- htmltools::tagAppendChild(tagged, notice(
          "Cancelled", "Press Run analysis to start again.", kind = "info"))
      } else if (!is.null(err)) {
        # The engine's own sentence when it is one: hiding "'control'
        # is not a group. Groups: Control, TreatA" behind "See the
        # technical details below" made users open a fold to read it.
        plain <- !looks_internal_error(err)
        tagged <- htmltools::tagAppendChild(
          tagged,
          notice(title  = if (is.null(diff_bundle())) "The differential analysis could not run"
                          else "The latest run failed; the result below is from the previous run",
                 detail = if (plain) err else diff_error_hint(err),
                 kind   = "error",
                 technical = if (!plain) err)
        )
      }
      # What the engine said about the result: an ignored covariate, a
      # scale conversion, genes too low to test.
      warns <- diff_bundle()$warnings
      if (length(warns)) {
        tagged <- htmltools::tagAppendChild(tagged, notice(
          title = "Notes on this result",
          detail = htmltools::tags$ul(lapply(unique(warns), htmltools::tags$li)),
          kind = "info"))
      }
      if (is.null(diff_bundle()) && is.null(err)) {
        tagged <- htmltools::tagAppendChild(
          tagged,
          notice(title  = "No result yet",
                 detail = paste("Check the layer, the control group and the groups to",
                                "compare in the Parameters card, then press Run analysis."),
                 kind   = "info")
        )
      }
      if (length(missing_engines)) {
        tagged <- htmltools::tagAppendChild(
          tagged,
          notice(
            title  = "Some engines are unavailable",
            detail = sprintf(
              "Not installed: %s. Install with `omicsCore::install_optional()`.",
              paste(missing_engines, collapse = ", ")
            ),
            kind = "info"
          )
        )
      }
      tagged
    })

    output$stats <- shiny::renderUI({
      b <- shown_bundle()
      if (is.null(b)) return(NULL)
      # A weighted contrast has no single case group; its "up" is the
      # direction of the contrast as written.
      case_lbl <- b$params$case_group %||% "contrast"
      if (length(case_lbl) > 1L) case_lbl <- "case"
      df <- marked()
      sig <- df[df$is_significant, , drop = FALSE]
      up_n   <- sum(sig$effect > 0, na.rm = TRUE)
      down_n <- sum(sig$effect < 0, na.rm = TRUE)
      top    <- if (nrow(sig) > 0L) sig[which.max(abs(sig$effect)), ] else NULL
      top_value <- if (is.null(top)) "\u2014" else as.character(top$feature_symbol[1L])
      top_trend <- if (is.null(top)) "no features pass thresholds"
                   else sprintf("%s %+.2f \u00B7 %s %.2g", omicsCore::effect_label(b),
                                top$effect[1L], p_label(), top[[p_col()]][1L])
      htmltools::tags$div(
        class = "stat-grid",
        {
          # Features with no p-value were not tested (too many missing
          # values, or set aside as too low to test); counting them as
          # tested overstated the screen.
          n_tested <- sum(!is.na(df$p_value))
          stat_card(
            label = "Tested features",
            value = if (n_tested < nrow(df))
                      sprintf("%s / %s", format(n_tested, big.mark = ","), format(nrow(df), big.mark = ","))
                    else format(nrow(df), big.mark = ","),
            trend = if (n_tested < nrow(df))
                      sprintf("%s not testable (missing values or too few counts)",
                              format(nrow(df) - n_tested, big.mark = ","))
                    else sprintf("%s, %s", b$params$method, b$params$comparison %||% "\u2014"),
            mono  = TRUE
          )
        },
        stat_card(
          label  = sprintf("Up in %s", case_lbl),
          value  = up_n,
          trend  = sprintf("%s > %.2f \u00B7 %s < %.3g", omicsCore::effect_label(b),
                           fc_cut_d(), p_label(), fdr_cut_d()),
          accent = "up"
        ),
        stat_card(
          label  = sprintf("Down in %s", case_lbl),
          value  = down_n,
          trend  = sprintf("%s < -%.2f \u00B7 %s < %.3g", omicsCore::effect_label(b),
                           fc_cut_d(), p_label(), fdr_cut_d()),
          accent = "down"
        ),
        stat_card(
          label = "Top hit",
          value = top_value,
          trend = top_trend,
          mono  = TRUE
        )
      )
    })

    output$volcano <- plotly::renderPlotly({
      b <- shown_bundle()
      shiny::validate(shiny::need(b, "Press Run analysis to draw the volcano."))
      # Deliberately not given the slider values. The volcano is drawn
      # at plot_volcano()'s own defaults, which is what an exported
      # report and an exported script also produce -- so the figure a
      # reader is shown is the figure they can reproduce, and a
      # screenshot does not depend on where a control happened to be.
      #
      # The sliders still drive the hit table and the stat cards, where
      # sweeping a threshold is the useful thing to do; the figure is
      # the stable reference next to them.
      p <- omicsCore::plot_volcano(
        b,
        top_n = if (isTRUE(input$label_top)) 20L else 0L
      )
      # WebGL rather than one SVG node per point: 60,000 genes painted in
      # 0.5 s instead of 4 s, with the same points and hover text.
      plotly::ggplotly(p, tooltip = "text") |>
        drop_hoveron() |>
        plotly::toWebGL() |>
        plotly::config(displaylogo = FALSE,
                       modeBarButtonsToRemove = c("lasso2d", "select2d"))
    })

    # An empty card before the first run read as a broken one.
    output$hits_empty <- shiny::renderUI({
      if (!is.null(diff_bundle())) return(NULL)
      htmltools::tags$p(class = "muted", style = "font-size:13px;margin:4px 0",
                        "Run the analysis to see the features that pass the thresholds.")
    })

    output$hits <- DT::renderDT({
      df <- marked()
      sig <- df[df$is_significant, , drop = FALSE]
      sig <- sig[order(-abs(sig$effect)), , drop = FALSE]
      out <- data.frame(
        Feature   = sig$feature_symbol,
        Effect    = round(sig$effect, 3),
        p         = signif(sig[[p_col()]], 3),
        Direction = ifelse(sig$effect > 0, "up", "down"),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      # The columns are named for what they hold: the p-value the mask
      # was read from, and the effect in the words the cards use.
      names(out)[2] <- omicsCore::effect_label(shown_bundle())
      names(out)[3] <- p_label()
      DT::datatable(
        out,
        rownames  = FALSE,
        selection = "single",
        options   = list(
          pageLength = 10,
          dom        = "ftip",
          language   = list(emptyTable = "No feature passes the current thresholds."),
          scrollX    = TRUE,
          columnDefs = list(list(className = "dt-right", targets = 1:2))
        )
      )
    }, server = TRUE)

    # The bundle *and* the thresholds it is read at. Enrichment used to
    # take only the bundle and apply its own defaults -- adjusted p at
    # 0.05, no fold-change bound -- so the gene list it enriched was a
    # different set from the hits shown here, and neither view said so.
    # A user with a hundred hits by raw p got an enrichment over however
    # many passed adj.P, which can be none.
    # ---- several contrasts side by side ------------------------------
    output$contrast_summary <- shiny::renderUI({
      if (length(comparisons()) < 2L) return(NULL)
      bslib::card(
        bslib::card_header(
          htmltools::tags$h3(class = "card-title", "Comparisons"),
          htmltools::tags$span(class = "card-sub",
                               "hits per comparison at the current thresholds")
        ),
        bslib::card_body(
          htmltools::tags$div(
            class = "row-grid r-6-6",
            shiny::plotOutput(session$ns("contrast_plot"),
                              height = paste0(160 + label_rows_px(comparisons(), 20L), "px")),
            DT::DTOutput(session$ns("contrast_table"))
          ),
          # Which comparisons share their hits. The table's "also in
          # another" says how many; this says with which.
          htmltools::tags$div(
            class = "inline-control",
            shiny::radioButtons(session$ns("overlap_dir"), label = "Overlap of",
                                choices = c("all hits" = "any", "up" = "up",
                                            "down" = "down"),
                                selected = "any", inline = TRUE)
          ),
          shiny::plotOutput(session$ns("overlap_plot"),
                            height = paste0(270 + label_rows_px(comparisons(), 10L), "px"))
        )
      )
    })

    output$overlap_plot <- shiny::renderPlot(res = PLOT_RES, alt = "Overlap of the significant features between comparisons", fit_to_width("overlap_plot", {
      b <- diff_bundle()
      shiny::req(omicsCore::is_analysis_bundle(b), length(comparisons()) > 1L)
      omicsCore::plot_diff_overlap(
        b, p_cutoff = fdr_cut_d(),
        p_preference = if (identical(input$p_kind %||% "adj", "raw")) "raw" else "adjusted",
        effect_cutoff = fc_cut_d(),
        direction = input$overlap_dir %||% "any")
    }))

    contrast_summary_df <- shiny::reactive({
      b <- diff_bundle()
      shiny::req(omicsCore::is_analysis_bundle(b), length(comparisons()) > 1L)
      omicsCore::summarize_diff_contrasts(
        b, p_cutoff = fdr_cut_d(),
        p_preference = if (identical(input$p_kind %||% "adj", "raw")) "raw" else "adjusted",
        effect_cutoff = fc_cut_d())
    })

    output$contrast_plot <- shiny::renderPlot(res = PLOT_RES, alt = "Number of up- and down-regulated features per comparison", fit_to_width("contrast_plot", {
      b <- diff_bundle()
      shiny::req(omicsCore::is_analysis_bundle(b), length(comparisons()) > 1L)
      omicsCore::plot_diff_contrasts(
        b, p_cutoff = fdr_cut_d(),
        p_preference = if (identical(input$p_kind %||% "adj", "raw")) "raw" else "adjusted",
        effect_cutoff = fc_cut_d())
    }))

    output$contrast_table <- DT::renderDT({
      s <- contrast_summary_df()
      out <- data.frame(
        Comparison = gsub("_vs_", " vs ", s$comparison),
        Up = s$n_up, Down = s$n_down,
        `Also in another` = s$n_shared,
        check.names = FALSE, stringsAsFactors = FALSE)
      DT::datatable(out, rownames = FALSE, selection = "none",
                    options = list(dom = "t", pageLength = 50))
    }, server = TRUE)

    # ---- global test across all groups -------------------------------
    # "Does this feature differ between any of the groups?" -- one test per
    # feature over every group at once, before (or instead of) reading the
    # comparisons one by one. Kept beside the comparisons rather than
    # replacing them: it has no direction and no fold change, so nothing
    # downstream (enrichment, integration) can use it.
    anova_bundle <- shiny::reactiveVal(NULL)
    anova_error <- shiny::reactiveVal(NULL)
    anova_running <- shiny::reactiveVal(FALSE)

    output$anova_card <- shiny::renderUI({
      if (length(levels_()) < 3L) return(NULL)
      b <- anova_bundle()
      bslib::card(
        bslib::card_header(
          htmltools::tags$h3(class = "card-title", "Any difference between groups"),
          htmltools::tags$span(class = "card-sub",
                               "one global test per feature (ANOVA / LRT)"),
          info_tip(paste("Tests all groups of the column at once: a small p says the",
                         "feature differs somewhere among them, not where. Use it to",
                         "screen, then read the comparisons for the direction."))
        ),
        bslib::card_body(
          htmltools::tags$div(
            style = "display:flex;gap:12px;align-items:center;flex-wrap:wrap",
            disabled_if(shiny::actionButton(session$ns("run_anova"),
                                            if (is.null(b)) "Run global test" else "Re-run global test",
                                            class = "btn btn-sm btn-outline-primary"),
                        shiny::isolate(anova_running())),
            shiny::uiOutput(session$ns("anova_summary"), inline = TRUE)
          ),
          if (!is.null(anova_error())) {
            notice("The global test could not run", kind = "error",
                   technical = anova_error())
          },
          if (!is.null(b)) DT::DTOutput(session$ns("anova_table"), fill = FALSE),
          if (!is.null(b)) shiny::downloadButton(session$ns("download_anova"),
                                                 "Download all (CSV)",
                                                 class = "btn btn-sm btn-ghost")
        )
      )
    })

    shiny::observeEvent(input$run_anova, {
      a <- active()
      shiny::req(a$input)
      d <- default_contrast()
      group_col <- input$group_col %||% d$group_col
      method <- input$method %||% "auto"
      if (!method %in% c("limma", "edger", "deseq2")) method <- "auto"
      covariates <- input$covariates
      my_run <- anova_epoch$start()
      set_button_busy("run_anova", TRUE, anova_running)
      run_async(
        detached_call(
          function() {
            omicsCore::run_diff(input = inp, method = method,
                                analysis_type = "anova", group_col = group_col,
                                covariates = covariates)
          },
          inp = a$input, method = method, group_col = group_col,
          covariates = if (length(covariates)) covariates else NULL
        ),
        on_success = function(bundle) {
          if (anova_epoch$is_last_started(my_run)) set_button_busy("run_anova", FALSE, anova_running)
          if (!anova_epoch$is_current(my_run)) return(invisible())
          anova_error(NULL)
          anova_bundle(bundle)
        },
        on_error = function(msg) {
          if (anova_epoch$is_last_started(my_run)) set_button_busy("run_anova", FALSE, anova_running)
          if (anova_epoch$is_current(my_run)) anova_error(msg)
        },
        message = "Running the global test..."
      )
    })

    anova_hits <- shiny::reactive({
      b <- anova_bundle()
      shiny::req(b)
      df <- b$results$diff_result_df
      df[order(df[[p_col()]], na.last = TRUE), , drop = FALSE]
    })

    output$anova_summary <- shiny::renderUI({
      df <- anova_hits()
      n <- sum(df[[p_col()]] < fdr_cut_d(), na.rm = TRUE)
      htmltools::tags$span(
        htmltools::tags$strong(format(n, big.mark = ",")),
        sprintf(" of %s features differ between the groups of '%s' (%s < %.3f)",
                format(nrow(df), big.mark = ","), anova_bundle()$params$group_col,
                p_label(), fdr_cut_d()))
    })

    output$anova_table <- DT::renderDT({
      df <- anova_hits()
      out <- data.frame(
        Feature = df$feature_symbol %||% df$feature_id,
        Statistic = signif(df$statistic, 3),
        p = signif(df[[p_col()]], 3),
        check.names = FALSE, stringsAsFactors = FALSE)
      names(out)[2] <- df$statistic_type[1] %||% "Statistic"
      names(out)[3] <- p_label()
      DT::datatable(out, rownames = FALSE, selection = "none",
                    options = list(pageLength = 10, dom = "ftip"))
    }, server = TRUE)

    # The whole result, every feature and every comparison, with the
    # significance at the current thresholds: the table a user takes to
    # a paper or a colleague. The hit list on screen is one page of one
    # comparison.
    full_table <- function() {
      b <- diff_bundle()
      shiny::req(b)
      df <- b$results$diff_result_df
      pv <- df[[p_col()]]
      df$significant <- !is.na(pv) & !is.na(df$effect) &
        pv < fdr_cut_d() & abs(df$effect) >= fc_cut_d()
      keep <- intersect(c("comparison", "feature_id", "feature_symbol", "effect",
                          "effect_type", "statistic", "p_value", "adj_p_value",
                          "base_mean", "significant", "method"), names(df))
      df[, keep, drop = FALSE]
    }
    table_name <- function(ext) {
      sprintf("differential_%s_%s.%s", active()$tag %||% "layer",
              format(Sys.Date(), "%Y%m%d"), ext)
    }
    output$download_csv <- shiny::downloadHandler(
      filename = function() table_name("csv"),
      content = function(file) utils::write.csv(full_table(), file, row.names = FALSE,
                                                fileEncoding = "UTF-8")
    )
    output$download_xlsx <- shiny::downloadHandler(
      filename = function() table_name("xlsx"),
      content = function(file) {
        df <- full_table()
        sheets <- split(df, df$comparison %||% "result")
        names(sheets) <- substr(gsub("[\\[\\]*?/:]", "_", names(sheets)), 1L, 31L)
        openxlsx::write.xlsx(sheets, file)
      }
    )
    output$download_anova <- shiny::downloadHandler(
      filename = function() sprintf("global_test_%s.csv", active()$tag %||% "layer"),
      content = function(file) {
        b <- anova_bundle()
        shiny::req(b)
        utils::write.csv(b$results$diff_result_df, file, row.names = FALSE, fileEncoding = "UTF-8")
      }
    )

    output$run_button <- shiny::renderUI({
      btn <- shiny::actionButton(
        session$ns("rerun"),
        if (is.null(diff_bundle())) "Run analysis" else "Re-run",
        icon = shiny::icon("play"),
        class = "btn btn-primary", style = "width:100%")
      # Re-rendered enabled while a run was still going (the label
      # changes when a result lands or is cleared), which let a second
      # run start alongside the first.
      if (shiny::isolate(running())) btn <- htmltools::tagAppendAttributes(btn, disabled = NA)
      btn
    })

    if (is.function(navigate)) {
      shiny::observeEvent(input$go_enrich, navigate("enrich"))
    }

    list(
      bundle = shown_bundle,
      # Every contrast of the last run, for anything that wants them all.
      all_bundle = shiny::reactive(diff_bundle()),
      # What the project keeps: every contrast, and which one the other
      # views were shown -- so the report covers them all and the script
      # can take the same one out for the steps that used it.
      project_bundle = shiny::reactive({
        b <- diff_bundle()
        if (!omicsCore::is_analysis_bundle(b)) return(b)
        # The thresholds the hits were read at, for the report: it said
        # "adjusted p < 0.05" whatever the user had chosen.
        b$params$display_thresholds <- list(
          p_cutoff = fdr_cut_d(),
          p_preference = if (identical(input$p_kind %||% "adj", "raw")) "raw" else "adjusted",
          effect_cutoff = fc_cut_d())
        if (length(comparisons()) >= 2L) b$params$shown_comparison <- shown_bundle()$params$comparison
        b
      }),
      # The global test, kept with the project like the comparisons.
      anova = shiny::reactive(anova_bundle()),
      # The layer this ran on. The bundle does not carry it -- a project
      # layer tag is an app concept, not something run_diff() knows --
      # and Enrichment needs it to say which layer its pathways came
      # from rather than leaving the reader to assume.
      layer = shiny::reactive(active()$tag),
      thresholds = shiny::reactive(list(
        p_cutoff      = fdr_cut_d(),
        p_preference  = if (identical(input$p_kind %||% "adj", "raw")) "raw"
                        else "adjusted",
        effect_cutoff = fc_cut_d()
      ))
    )
  })
}

# ---- internal helpers ------------------------------------------------

utils::globalVariables(".data")

`%||%` <- function(a, b) if (is.null(a)) b else a

# Which Bioconductor diff backends aren't installed in this R
# session. Returned as a character vector for the notices strip.
# DESeq2 / edgeR / limma are the ones run_diff() can dispatch to.
diff_missing_engines <- function() {
  engines <- c(limma = "limma", DESeq2 = "DESeq2", edgeR = "edgeR")
  missing <- vapply(engines, function(pkg) !has_pkg(pkg), logical(1))
  names(engines)[missing]
}

diff_omics_display <- function(t) {
  switch(t %||% "",
         proteomics = "Proteomics",
         rnaseq     = "RNA-seq",
         "\u2014")
}

diff_params_card <- function(ns) {
  bslib::card(
    class = "params-card",
    bslib::card_header(
      htmltools::tags$h3(class = "card-title", "Parameters"),
      htmltools::tags$span(class = "card-sub",
                           "design + thresholds")
    ),
    bslib::card_body(
      htmltools::tags$div(
        class = "param-stack",
        param_group(
          "Layer",
          # Which engines are on offer follows from this: deseq2 and
          # edger need rnaseq raw counts, so without a way to change
          # layer they were unreachable.
          shiny::uiOutput(ns("ui_layer"))
        ),
        param_group(
          "Method",
          help = paste("auto picks limma for proteomics and log-scale RNA-seq,",
                       "DESeq2 for raw counts. Engines that do not fit the",
                       "layer's data are hidden."),
          # Rendered server-side: which engines are valid depends on the
          # active layer's assay, and offering an invalid one produces a
          # complete, plausible, meaningless result table.
          shiny::uiOutput(ns("ui_method")),
          shiny::uiOutput(ns("method_note"))
        ),
        param_group(
          "Contrast",
          help = paste("The control is the reference every other group is compared",
                       "with. Choose several groups to fit them in one model;",
                       "'Showing' then picks the comparison on screen."),
          shiny::uiOutput(ns("ui_design_mode")),
          shiny::uiOutput(ns("ui_continuous")),
          shiny::uiOutput(ns("ui_group_col")),
          shiny::uiOutput(ns("ui_paired")),
          shiny::uiOutput(ns("ui_contrast")),
          shiny::uiOutput(ns("ui_comparison"))
        ),
        param_group(
          "Covariates",
          help = paste("Variables to adjust for (age, sex, batch). The group effect",
                       "is then estimated holding them constant. Leave empty for",
                       "an unadjusted comparison. For before/after or matched",
                       "samples use 'Paired by' above instead."),
          shiny::uiOutput(ns("ui_covariates"))
        ),
        param_group(
          "Thresholds",
          help = paste("Which features count as hits. They filter the table and the",
                       "counts, and are what Enrichment and Integration use; they",
                       "do not re-run the model."),
          # Which p to threshold on is the user's call, not ours. An
          # exploratory screen on 50 proteins and a confirmatory one on
          # 20,000 genes want different answers, and forcing adj.P made
          # the first look empty.
          shiny::radioButtons(
            ns("p_kind"), label = "Significance on",
            choices = c("adjusted p" = "adj", "raw p" = "raw"),
            selected = "adj", inline = TRUE
          ),
          # Ticks off: on a rail this narrow their labels ran together
          # ("0.020.04"); the handle shows the value.
          shiny::sliderInput(
            ns("fdr_cut"), label = "p cutoff",
            min = 0, max = 0.2, value = 0.05, step = 0.005, ticks = FALSE
          ),
          # A box, not a slider. The slider stepped 0.05, which cannot
          # express log2(1.2) = 0.263 -- so the fold change most often
          # wanted here was one of the few the control could not reach.
          shiny::numericInput(
            ns("fc_cut"), label = "|log2FC| cutoff",
            value = round(log2(1.2), 3), min = 0, max = 10, step = 0.05
          ),
          htmltools::tags$div(
            class = "muted", style = "font-size:11.5px;margin-top:-6px",
            sprintf("%.3f = %.2gx fold change", log2(1.2), 1.2)
          ),
          shinyWidgets::materialSwitch(
            ns("label_top"), label = "Label top 20",
            value = FALSE, status = "primary", right = TRUE
          )
        ),
        # Pinned to the bottom of the window while the rail scrolls: at
        # the end of the parameters it sat at y = 1,470 on a 900 px
        # screen, below the fold, under the controls it acts on.
        htmltools::tags$div(
          class = "run-sticky",
          shiny::uiOutput(ns("run_button"))
        )
      )
    )
  )
}

diff_volcano_card <- function(ns) {
  bslib::card(
    bslib::card_header(
      htmltools::tags$h3(class = "card-title", "Volcano"),
      # Says plainly that the thresholds do not reach this figure. The
      # cut it was drawn at is in the plot's own caption, so a
      # screenshot carries it too.
      htmltools::tags$span(
        class = "card-sub",
        "fixed thresholds \u00B7 sliders filter the table below")
    ),
    bslib::card_body(
      plotly::plotlyOutput(ns("volcano"), height = "360px"),
      htmltools::tags$div(
        class = "legend",
        legend_swatch("significant", omics_colors$up),
        legend_swatch("ns", omics_colors$ns)
      )
    )
  )
}

diff_hits_card <- function(ns) {
  bslib::card(
    bslib::card_header(
      htmltools::tags$h3(class = "card-title", "Top hits"),
      htmltools::tags$span(class = "card-sub",
                           "largest changes first, within current thresholds")
    ),
    bslib::card_body(
      shiny::uiOutput(ns("hits_empty")),
      DT::DTOutput(ns("hits"), fill = FALSE),
      htmltools::tags$div(
        style = "display:flex;gap:8px;margin-top:8px;flex-wrap:wrap",
        shiny::downloadButton(ns("download_csv"), "Download all results (CSV)",
                              class = "btn btn-sm btn-ghost"),
        shiny::downloadButton(ns("download_xlsx"), "Excel, one sheet per comparison",
                              class = "btn btn-sm btn-ghost")
      )
    )
  )
}

# A group name as a contrast expression needs it.
contrast_token <- function(x) {
  x <- as.character(x %||% "B")
  if (identical(make.names(x), x)) x else paste0("`", x, "`")
}

# The level a control group is usually called, so the default reference
# is the control rather than whichever label sorts first ("DMSO" before
# "Drug", but "Treated" after "Control" only by luck).
CONTROL_LEVEL_PATTERN <- paste0(
  "^(ctrl|control|controls|con|dmso|vehicle|veh|mock|sham|untreated|",
  "placebo|baseline|wt|wild[ _-]?type|normal|healthy|nc|neg|negative|",
  "blank|t0|0h|day0|d0|before|pre|pretreatment|pre[ _-]?treatment|",
  "pbs|saline|sinc|sictrl|shctrl|scramble|scrambled|",
  "\u5bf9\u7167|\u5bf9\u7167\u7ec4|\u7a7a\u767d|\u7a7a\u767d\u7ec4|",
  "\u6b63\u5e38|\u6b63\u5e38\u7ec4|\u9634\u6027\u5bf9\u7167)$")

default_control_level <- function(levels) {
  if (!length(levels)) return(NULL)
  hit <- grepl(CONTROL_LEVEL_PATTERN, trimws(tolower(levels)))
  if (!any(hit)) {
    hit <- grepl("(^|[ _-])(ctrl|control|vehicle|dmso|wt|before|pre|baseline)([ _-]|$)",
                 tolower(levels)) | grepl("\u5bf9\u7167", levels)
  }
  levels[if (any(hit)) which(hit)[1L] else 1L]
}

# An R-internal error (a subscript, a missing object) rather than a
# sentence written for the user.
looks_internal_error <- function(msg) {
  grepl(paste0("Error in|subscript|object '.*' not found|non-numeric argument|",
               "argument is of length zero|missing value where|unused argument|",
               "could not find function|values must be length|replacement has|",
               "non-conformable|NA/NaN/Inf|invalid 'type'|undefined columns"),
        msg %||% "")
}

# A plain-language reading of the errors a run most often ends in; the
# message itself stays available under "Technical details".
diff_error_hint <- function(msg) {
  msg <- msg %||% ""
  if (grepl("not a level of", msg, fixed = TRUE)) {
    return("One of the chosen groups is not in this layer's group column.")
  }
  if (grepl("confounded", msg, fixed = TRUE)) {
    return("A covariate cannot be separated from the groups being compared; remove it.")
  }
  if (grepl("missing values", msg, fixed = TRUE)) {
    return("A covariate or pairing column has empty cells for some samples.")
  }
  if (grepl("paired design", msg, fixed = TRUE)) {
    return("The pairing column does not pair the samples one-to-one across the groups.")
  }
  if (grepl("required|not installed", msg)) {
    return("The chosen method needs a package that is not installed; pick another method.")
  }
  "See the technical details below."
}

# Which layer the view lands on when the user has not chosen one:
# proteomics if present, else the first.
default_layer_tag <- function(experiments) {
  if (!length(experiments)) return(NULL)
  types <- vapply(experiments, function(e) e$omics_type %||% "", character(1))
  i <- which(types == "proteomics")
  names(experiments)[if (length(i)) i[1L] else 1L]
}

# Columns of `meta_df` that could name a contrast, best first.
#
# The old rule was "not numeric, and at least two distinct values",
# which a sample identifier satisfies perfectly: one level per sample.
# On a real workbook the first such column was `label`, whose values are
# the sample names, so the view defaulted to a contrast of one sample
# against one other. limma cannot fit that -- no residual degrees of
# freedom -- and reported it as "Partial NA coefficients for 2294
# probe(s)" and an empty result, which is not a sentence anyone can act
# on.
#
# The real requirement is replication: every level needs at least two
# samples, or the level cannot be tested. That single condition
# excludes identifiers exactly, without having to guess from names.
GROUP_COL_HINTS <- c("group", "condition", "treatment", "arm", "status",
                     "genotype", "diet", "disease", "phenotype", "diagnosis",
                     "cohort", "timepoint", "time_point", "visit",
                     "\u5206\u7ec4", "\u7ec4\u522b", "\u7ec4", "\u5904\u7406")
# Columns that describe samples but are rarely the comparison: ranked last.
NUISANCE_COL_RE <- paste0("^(batch|sex|gender|cage|plate|lane|run|replicate|rep|",
                          "patient|donor|subject|individual|animal|mouse|pair|",
                          "\u6279\u6b21|\u6027\u522b)([ _.-]?(id|no))?$")
# Columns that pair samples, preselected in 'Paired by'.
PAIRING_COL_HINTS <- paste0("^(patient|donor|subject|individual|animal|mouse|pair|",
                            "participant|case)([ _.-]?(id|no))?$|^(\u60a3\u8005|\u4e2a\u4f53)")

grouping_candidates <- function(meta, min_per_level = 2L, replicated = FALSE) {
  if (is.null(meta) || !ncol(meta)) return(character(0))
  usable <- vapply(names(meta), function(nm) {
    col <- meta[[nm]]
    # A numeric column is a grouping when it has a few repeated values:
    # 0/1 coding, dose levels. Age or BMI has as many values as samples.
    if (is.numeric(col) && length(unique(stats::na.omit(col))) > 6L) return(FALSE)
    counts <- table(as.character(col), useNA = "no")
    length(counts) >= 2L && min(counts) >= min_per_level &&
      (!replicated || max(counts) >= 2L)
  }, logical(1))
  cands <- names(meta)[usable]
  if (!length(cands)) return(character(0))

  # Fewest levels first: a two-level column is the contrast someone
  # almost always means. Conventional names win over the count, since a
  # column called `condition` is a stated intent and a level count is an
  # inference.
  n_levels <- vapply(cands, function(nm) {
    length(unique(stats::na.omit(as.character(meta[[nm]]))))
  }, integer(1))
  low <- tolower(cands)
  hinted <- low %in% GROUP_COL_HINTS | grepl("group|condition|treat", low)
  nuisance <- grepl(NUISANCE_COL_RE, low)
  cands[order(!hinted, nuisance, n_levels)]
}

# ggplotly() sets `hoveron` on its traces; scattergl has no such
# attribute, so once toWebGL() converts them plotly warns about it on
# every build. Dropped first, it is never there to warn about.
drop_hoveron <- function(fig) {
  fig$x$data <- lapply(fig$x$data, function(tr) {
    tr$hoveron <- NULL
    tr
  })
  fig
}
