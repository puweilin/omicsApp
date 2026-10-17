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
#' When the project is `NULL`, the view falls back to
#' `example_diff_bundle()` so the volcano + hits stay visible.
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
        diff_volcano_card(ns),
        diff_hits_card(ns)
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
                             activation = shiny::reactiveVal(0L)) {
  shiny::moduleServer(id, function(input, output, session) {

    # A restored result can arrive before the browser binds its layer
    # dropdown. Keep that handshake separate from a user changing layers.
    restoring_layer <- shiny::reactiveVal(NULL)
    diff_bundle <- shiny::reactiveVal(NULL)
    diff_error <- shiny::reactiveVal(NULL)
    restored_params <- shiny::reactive({
      b <- diff_bundle()
      if (is.null(b$provenance_binding)) return(list())
      tag <- restoring_layer() %||% input$layer
      if (!is.null(tag) && !identical(names(b$provenance_binding$inputs), tag)) return(list())
      b$params %||% list()
    })

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
    active <- shiny::reactive({
      proj <- current_project()
      is_demo <- is.null(proj)
      if (is_demo) proj <- example_project()
      exps <- proj$experiments
      if (length(exps) == 0L) {
        return(list(input = NULL, tag = NULL, is_demo = TRUE))
      }
      want <- restoring_layer() %||% input$layer
      tag <- if (!is.null(want) && want %in% names(exps)) {
        want
      } else {
        default_layer_tag(exps)
      }
      # Not exps[[tag]]. QC's filtering and imputation produced a cleaned
      # layer that this view never read, so a user who removed high-missing
      # features and imputed the rest watched the differential run on the
      # raw import anyway -- with the QC panel showing the cleaned numbers.
      # effective_input() is the one place that answers which data a layer
      # means, and every consumer goes through it.
      eff <- omicsCore::effective_input(proj, tag)
      list(input = eff$input, tag = tag, source = eff$source,
           version = eff$version, stale = eff$stale, is_demo = is_demo)
    })

    output$ui_layer <- shiny::renderUI({
      proj <- current_project() %||% example_project()
      tags_avail <- names(proj$experiments)
      if (length(tags_avail) == 0L) return(NULL)
      # isolate(), because active() reads input$layer and this output
      # writes it. Reading it here closed the loop: re-rendering the
      # control re-sent its value, which invalidated active(), which
      # re-rendered the control.
      sel <- restoring_layer() %||% shiny::isolate(input$layer)
      if (is.null(sel) || !sel %in% tags_avail) {
        sel <- default_layer_tag(proj$experiments)
      }
      shiny::selectInput(
        session$ns("layer"), label = "Experiment layer",
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
      choices <- omicsCore::applicable_diff_methods(inp)
      selected <- if (!is.null(restoring_layer())) restored_params()$method else shiny::isolate(input$method)
      selected <- selected %||% restored_params()$method %||% "auto"
      if (!is.character(selected) || length(selected) != 1L || is.na(selected) ||
          !selected %in% choices) selected <- "auto"
      shiny::selectInput(session$ns("method"), label = NULL,
                         choices = choices, selected = selected)
    })

    output$method_note <- shiny::renderUI({
      a <- active()
      inp <- a$input
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
      # No fallback to `names(meta)`. Offering every column when none
      # qualifies handed a metadata sheet whose only column was
      # `sample_id` back as a group column -- nine levels, one sample
      # each, an analysis that cannot run -- while the view claimed a
      # choice existed. Empty here is the truth, and ui_group_col below
      # says so out loud.
      cands <- grouping_candidates(meta)
      gc <- if (length(cands)) cands[1L] else NULL
      lv <- if (is.null(gc)) character(0)
            else sort(unique(as.character(stats::na.omit(meta[[gc]]))))
      list(
        group_col  = gc,
        control    = if (length(lv)) lv[1L] else NULL,
        case       = if (length(lv) >= 2L) lv[2L] else NULL,
        levels     = lv,
        candidates = cands
      )
    })

    # Group column dropdown: any meta_df column with >= 2 unique
    # non-NA values (continuous columns like `age` are excluded
    # for the simple "control vs case" UI in this slice).
    # What the sample names would give, if the metadata gives nothing.
    #
    # A suggestion, never applied on its own: deriving an experimental
    # design from a naming convention and acting on it unasked is
    # asserting what the experiment was.
    name_grouping <- shiny::reactive({
      a <- active()
      if (is.null(a$input)) return(NULL)
      if (length(default_contrast()$candidates)) return(NULL)
      tryCatch(omicsCore::suggest_grouping_from_names(a$input),
               error = function(e) NULL)
    })

    output$ui_group_col <- shiny::renderUI({
      d <- default_contrast()
      if (!length(d$candidates)) {
        # Nothing to pick, so do not ask. The view used to say "pick a
        # group column with at least two levels" to a user who had none --
        # a message pointing at an option that does not exist reads as
        # "you are holding it wrong", and one user spent an afternoon
        # believing the analysis was slow when it had never started.
        sg <- name_grouping()
        return(htmltools::tagList(
          htmltools::tags$div(
            class = "muted", style = "font-size:12px",
            paste0("This layer's sample metadata has no column that can ",
                   "group samples. A grouping column needs at least two ",
                   "values, with at least two samples in each. The import ",
                   "template calls it `condition`.")),
          if (!is.null(sg)) htmltools::tagList(
            htmltools::tags$div(
              class = "muted",
              style = "font-size:12px;margin-top:8px;color:var(--brand-600)",
              sprintf("The sample names look like %d groups: %s.",
                      sg$n_levels,
                      paste(utils::head(sg$levels, 6), collapse = ", "))),
            shiny::actionButton(
              session$ns("use_name_grouping"),
              "Use the sample names as groups",
              class = "btn btn-sm btn-primary",
              style = "margin-top:6px"))
        ))
      }
      # Keep what the user picked. This output re-renders whenever the
      # layer or the project changes, and re-rendering with the default
      # threw away their choice -- they selected condition / G1 / G2,
      # ran it, and the control came back reading `label`.
      sel <- if (!is.null(restoring_layer())) restored_params()$group_col else shiny::isolate(input$group_col)
      sel <- sel %||% restored_params()$group_col
      if (is.null(sel) || !sel %in% d$candidates) sel <- d$group_col
      shiny::selectInput(session$ns("group_col"),
                         label    = "Group column",
                         choices  = d$candidates,
                         selected = sel)
    })

    # Confirmed, so it goes into the layer -- and into the layer's
    # provenance, so the exported script rebuilds the same column rather
    # than reading a file that never had it.
    #
    # Writing the project from an analysis view is unusual and deliberate:
    # this *is* a data-preparation step, and the alternative on offer was
    # "go back to Excel, add a column, re-import, lose your QC".
    shiny::observeEvent(input$use_name_grouping, {
      proj <- current_project()
      a <- active()
      if (is.null(proj) || is.null(a$tag)) return()
      updated <- tryCatch(
        omicsCore::add_grouping_from_names(a$input, column = "group"),
        error = function(e) e)
      if (inherits(updated, "error")) {
        shiny::showNotification(conditionMessage(updated), type = "error",
                                duration = 8)
        return()
      }
      proj$experiments[[a$tag]] <- updated
      current_project(proj)
      shiny::showNotification(
        sprintf("Added a `group` column from the sample names: %d groups.",
                length(unique(updated$meta_df$group))),
        type = "message", duration = 6)
    })

    # Reactive level set for the chosen group column.
    levels_ <- shiny::reactive({
      a <- active()
      meta <- a$input$meta_df
      gc <- if (!is.null(restoring_layer())) restored_params()$group_col else input$group_col
      gc <- gc %||% restored_params()$group_col %||% default_contrast()$group_col
      if (is.null(gc) || !(gc %in% names(meta))) return(character(0))
      sort(unique(as.character(stats::na.omit(meta[[gc]]))))
    })

    output$ui_contrast <- shiny::renderUI({
      lv <- levels_()
      if (length(lv) < 2L) {
        # Only reachable with a column chosen: with none to choose,
        # ui_group_col above says so and offers what it can.
        if (!length(default_contrast()$candidates)) return(NULL)
        return(htmltools::tags$div(
          class = "muted",
          style = "font-size:12px",
          "Pick a group column with at least two levels."
        ))
      }
      keep <- function(current, saved, fallback) {
        current <- if (!is.null(restoring_layer())) saved else shiny::isolate(current)
        current <- current %||% saved
        if (is.null(current) || !current %in% lv) fallback else current
      }
      htmltools::tagList(
        shiny::selectInput(session$ns("control"),
                           label = "Control", choices = lv,
                           selected = keep(input$control, restored_params()$control_group, lv[1L])),
        shiny::selectInput(session$ns("case"),
                           label = "Case", choices = lv,
                           selected = keep(input$case, restored_params()$case_group, lv[min(2L, length(lv))]))
      )
    })

    # Which backend will actually run, so the controls can describe it. The
    # dropdown may say "auto", and auto resolves to something specific.
    resolved_method <- shiny::reactive({
      a <- active()
      m <- if (!is.null(restoring_layer())) restored_params()$method else input$method
      m <- m %||% restored_params()$method %||% "auto"
      if (!identical(m, "auto") || is.null(a$input)) return(m)
      # applicable_diff_methods() puts auto first; the one after it is what
      # auto would pick for this layer.
      setdiff(omicsCore::applicable_diff_methods(a$input), "auto")[[1L]] %||% m
    })

    output$ui_covariates <- shiny::renderUI({
      a <- active()
      # A control that is enabled and ignored is worse than one that is
      # absent: a t-test silently dropped whatever was selected here and
      # returned the unadjusted fit. run_diff() now refuses that outright,
      # so offering the control would only produce an error the user could
      # not have predicted from the interface.
      if (identical(resolved_method(), "ttest")) {
        return(htmltools::tags$div(
          class = "muted", style = "font-size:12px",
          "A t-test fits no covariates. Switch to limma or lm to adjust for one."
        ))
      }
      meta <- a$input$meta_df
      gc <- input$group_col %||% ""
      cands <- setdiff(names(meta), c(gc, "sample_id"))
      chosen <- shiny::isolate(input$covariates)
      if (!is.null(restoring_layer()) || !"covariates" %in% shiny::isolate(names(input))) {
        chosen <- restored_params()$covariates
      }
      shiny::selectizeInput(
        session$ns("covariates"),
        label    = NULL,
        choices  = cands,
        multiple = TRUE,
        selected = intersect(chosen %||% character(0), cands),
        options  = list(placeholder = "optional, e.g. age")
      )
    })

    # Diff bundle: demo fallback when no project; otherwise gated
    # behind the Re-run button. We also auto-run once on first
    # mount when a real project is present, so the user lands on
    # a populated volcano without a Re-run click. After that,
    # changes only take effect on Re-run.
    # The layer this result was computed on has been replaced, so the
    # result is no longer about anything in the project. NULL is the
    # module's own start-up state, so this only rewinds it.
    shiny::observeEvent(invalidate(), {
      restoring_layer(NULL)
      diff_bundle(NULL)
      diff_error(NULL)
    }, ignoreInit = TRUE)

    # A project arrives carrying results it was saved with. Clearing the
    # views on activation is right -- what they held belongs to the last
    # project -- but stopping there left the page empty in front of a
    # project whose differential is on disk and valid, so the user's only
    # move was to re-run it and hope the settings matched.
    #
    # Only results the project itself vouches for: a binding that names
    # this project and inputs that are still current. Anything else is
    # what retire_invalid() has already moved to history.
    shiny::observeEvent(activation(), {
      proj <- current_project()
      if (is.null(proj)) return()
      saved <- proj$bundles$diff
      if (is.null(saved) || is.null(saved$provenance_binding)) return()
      if (!omicsCore::result_inputs_current(saved$provenance_binding, proj)) return()
      diff_error(NULL)
      diff_bundle(saved)
      tag <- names(saved$provenance_binding$inputs)[[1L]]
      restoring_layer(tag)
      shiny::updateSelectInput(session, "layer", selected = tag)
      # And the controls, so the panel describes the result it is showing
      # rather than the defaults. A page that shows one contrast above
      # settings for another invites re-running the wrong thing.
      p <- saved$params %||% list()
      if (!is.null(p$method)) {
        shiny::updateSelectInput(session, "method", selected = p$method)
      }
      if (!is.null(p$group_col)) {
        shiny::updateSelectInput(session, "group_col", selected = p$group_col)
      }
      if (!is.null(p$control_group)) {
        shiny::updateSelectInput(session, "control", selected = p$control_group)
      }
      if (!is.null(p$case_group)) {
        shiny::updateSelectInput(session, "case", selected = p$case_group)
      }
      if (!is.null(p$covariates)) {
        shiny::updateSelectizeInput(session, "covariates", selected = p$covariates)
      }
    }, priority = -20L, ignoreInit = FALSE)

    # One number per view, bumped when a run starts and checked when it
    # comes back. Pressing Run again while the previous run is still out
    # used to mean whichever finished last won -- and the loser was often
    # the newer one, because the settings that made it slower were the
    # settings the user had just chosen.
    request_seq <- shiny::reactiveVal(0L)

    do_run <- function() {
      request_seq(shiny::isolate(request_seq()) + 1L)
      this_request <- shiny::isolate(request_seq())
      a <- active()
      if (is.null(a$input)) {
        diff_error("This project has no experiments to analyse.")
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
      covariates <- input$covariates
      # Shiny keeps an input's value when its control stops being rendered,
      # so switching to a t-test after selecting a covariate still sends
      # one. The control is hidden above; this is what makes it not apply.
      if (identical(method, "ttest") ||
          (identical(method, "auto") && identical(resolved_method(), "ttest"))) {
        covariates <- NULL
      }
      if (is.null(group_col) || is.null(control) || is.null(case) ||
          identical(control, case)) {
        # Two different situations, and one message for both told the
        # user with no grouping column at all to go and pick one.
        diff_error(if (!length(d$candidates)) {
          paste0("This layer's sample metadata has no column that can group ",
                 "samples, so there is nothing to compare. Add one (the ",
                 "import template calls it `condition`), or use the ",
                 "suggestion from the sample names if the view offers it.")
        } else {
          "Pick a group column with distinct Control and Case levels."
        })
        return(invisible())
      }
      # Captured before the run starts, and never re-read afterwards: the
      # first version of this took the versions from here and the project
      # identity from arrival, which produced a record that looked
      # internally consistent while describing a computation from another
      # project.
      launch <- capture_launch(current_project(), "diff", tags = a$tag,
                               params = list(method = method,
                                             group_col = group_col,
                                             control_group = control,
                                             case_group = case,
                                             covariates = covariates),
                               activation = activation())
      run_async(
        admissible = launch_admissible(
          launch,
          function() list(project = current_project(),
                          activation = activation(),
                          request_seq = request_seq()),
          this_request),
        on_discard = discard_notice(),
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
              case_group    = case,
              covariates    = covariates
            )
          },
          inp        = a$input,
          method     = method,
          group_col  = group_col,
          control    = control,
          case       = case,
          covariates = if (length(covariates)) covariates else NULL
        ),
        on_success = function(bundle) {
          # From the launch record. Nothing here is read from the live
          # application: if the project has moved on, the gate above has
          # already refused, and describing the result by whatever it
          # landed in is how a stale result comes to look current.
          bundle <- stamp_from_launch(bundle, launch)
          diff_error(NULL)
          diff_bundle(bundle)
        },
        on_error = function(msg) {
          diff_error(msg)
        },
        message = "Running differential analysis..."
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
      pending <- restoring_layer()
      if (!is.null(pending)) {
        if (identical(input$layer, pending)) restoring_layer(NULL)
        return()
      }
      b <- diff_bundle()
      proj <- current_project()
      if (!is.null(b$provenance_binding) && !is.null(proj) &&
          identical(names(b$provenance_binding$inputs), input$layer) &&
          omicsCore::result_inputs_current(b$provenance_binding, proj)) return()
      diff_bundle(NULL)
      diff_error(NULL)
    }, ignoreInit = TRUE)

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
      if (identical(input$p_kind %||% "adj", "raw")) "p" else "adj.P")

    marked <- shiny::reactive({
      shiny::req(diff_bundle())
      df <- diff_bundle()$results$diff_result_df
      pv <- df[[p_col()]]
      df$is_significant <- !is.na(pv) &
                           !is.na(df$effect) &
                           pv < fdr_cut_d() &
                           abs(df$effect) > fc_cut_d()
      df
    })

    output$header <- shiny::renderUI({
      a <- active()
      b <- diff_bundle()
      # Which version of the layer, not just which layer. "QC-cleaned" and
      # "raw import" are different data, and the whole point of routing this
      # view through effective_input() is lost if the page does not say
      # which one it got.
      source_note <- htmltools::tags$span(
        class = "muted",
        if (a$is_demo) {
          "demo project (built-in)"
        } else if (identical(a$source, "qc_cleaned")) {
          sprintf("layer = %s \u00B7 QC-cleaned", a$tag)
        } else if (isTRUE(a$stale)) {
          # Preprocessing exists but was built on a layer that has since
          # been replaced, so the baseline is what is running.
          sprintf("layer = %s \u00B7 raw import (QC applies to a replaced layer)",
                  a$tag)
        } else {
          sprintf("layer = %s \u00B7 raw import (QC not applied)", a$tag)
        }
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

      comparison <- b$params$comparison %||%
        sprintf("%s vs %s",
                b$params$case_group %||% input$case %||% "case",
                b$params$control_group %||% input$control %||% "control")
      view_header(
        title    = "Differential",
        subtitle = htmltools::tagList(
          diff_omics_display(b$input_info$omics_type),
          htmltools::HTML(" &middot; "),
          comparison,
          htmltools::HTML(" &middot; "),
          b$params$method,
          htmltools::HTML(" &middot; "),
          source_note
        )
      )
    })

    output$notices <- shiny::renderUI({
      err <- diff_error()
      missing_engines <- diff_missing_engines()
      tagged <- htmltools::tagList()
      if (!is.null(err)) {
        tagged <- htmltools::tagAppendChild(
          tagged,
          notice(title  = "run_diff failed",
                 detail = err,
                 kind   = "warn")
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
      b <- diff_bundle()
      if (is.null(b)) return(NULL)
      df <- marked()
      sig <- df[df$is_significant, , drop = FALSE]
      up_n   <- sum(sig$effect > 0, na.rm = TRUE)
      down_n <- sum(sig$effect < 0, na.rm = TRUE)
      top    <- if (nrow(sig) > 0L) sig[which.max(abs(sig$effect)), ] else NULL
      top_value <- if (is.null(top)) "\u2014" else as.character(top$feature_symbol[1L])
      top_trend <- if (is.null(top)) "no features pass thresholds"
                   else sprintf("effect %+.2f \u00B7 %s %.2g",
                                top$effect[1L], p_label(), top[[p_col()]][1L])
      htmltools::tags$div(
        class = "stat-grid",
        stat_card(
          label = "Tested features",
          value = format(nrow(df), big.mark = ","),
          trend = sprintf("%s, %s",
                          b$params$method,
                          b$params$comparison %||% "\u2014"),
          mono  = TRUE
        ),
        stat_card(
          label  = sprintf("Up in %s", input$case %||% "case"),
          value  = up_n,
          trend  = sprintf("effect > %.2f \u00B7 %s < %.3f",
                           fc_cut_d(), p_label(), fdr_cut_d()),
          accent = "up"
        ),
        stat_card(
          label  = sprintf("Down in %s", input$case %||% "case"),
          value  = down_n,
          trend  = sprintf("effect < -%.2f \u00B7 %s < %.3f",
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
      b <- diff_bundle()
      shiny::req(b)
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
      plotly::ggplotly(p) |>
        plotly::config(displaylogo = FALSE,
                       modeBarButtonsToRemove = c("lasso2d", "select2d"))
    })

    output$hits <- DT::renderDT({
      df <- marked()
      sig <- df[df$is_significant, , drop = FALSE]
      sig <- sig[order(-abs(sig$effect)), , drop = FALSE]
      out <- data.frame(
        Feature   = sig$feature_symbol,
        Effect    = round(sig$effect, 3),
        `adj.P`   = signif(sig$adj_p_value, 3),
        Direction = ifelse(sig$effect > 0, "up", "down"),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      DT::datatable(
        out,
        rownames  = FALSE,
        selection = "single",
        options   = list(
          pageLength = 10,
          dom        = "ftip",
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
    list(
      bundle = shiny::reactive(diff_bundle()),
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
          # Rendered server-side: which engines are valid depends on the
          # active layer's assay, and offering an invalid one produces a
          # complete, plausible, meaningless result table.
          shiny::uiOutput(ns("ui_method")),
          shiny::uiOutput(ns("method_note"))
        ),
        param_group(
          "Contrast",
          shiny::uiOutput(ns("ui_group_col")),
          shiny::uiOutput(ns("ui_contrast"))
        ),
        param_group(
          "Covariates",
          shiny::uiOutput(ns("ui_covariates"))
        ),
        param_group(
          "Thresholds",
          # Which p to threshold on is the user's call, not ours. An
          # exploratory screen on 50 proteins and a confirmatory one on
          # 20,000 genes want different answers, and forcing adj.P made
          # the first look empty.
          shiny::radioButtons(
            ns("p_kind"), label = "Significance on",
            choices = c("adjusted p" = "adj", "raw p" = "raw"),
            selected = "adj", inline = TRUE
          ),
          shiny::sliderInput(
            ns("fdr_cut"), label = "p cutoff",
            min = 0, max = 0.2, value = 0.05, step = 0.005
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
        htmltools::tags$div(
          style = "margin-top:8px",
          shiny::actionButton(
            ns("rerun"), "Re-run",
            class = "btn btn-primary",
            style = "width:100%"
          )
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
        legend_swatch("ns", omics_colors$ns),
        htmltools::tags$span(class = "muted",
                             style = "font-size:12px",
                             "dashed lines mirror the slider thresholds")
      )
    )
  )
}

diff_hits_card <- function(ns) {
  bslib::card(
    bslib::card_header(
      htmltools::tags$h3(class = "card-title", "Top hits"),
      htmltools::tags$span(class = "card-sub",
                           "ranked by |effect| within current thresholds")
    ),
    bslib::card_body(
      DT::DTOutput(ns("hits"))
    )
  )
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
GROUP_COL_HINTS <- c("group", "condition", "treatment", "arm", "status")

grouping_candidates <- function(meta, min_per_level = 2L) {
  if (is.null(meta) || !ncol(meta)) return(character(0))
  usable <- vapply(names(meta), function(nm) {
    col <- meta[[nm]]
    if (is.numeric(col)) return(FALSE)
    counts <- table(as.character(col), useNA = "no")
    length(counts) >= 2L && min(counts) >= min_per_level
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
  hinted <- tolower(cands) %in% GROUP_COL_HINTS
  cands[order(!hinted, n_levels)]
}
