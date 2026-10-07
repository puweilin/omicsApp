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

    # The layer a restored result was computed on, until the control exists.
    preferred_layer <- shiny::reactiveVal(NULL)
    skip_layer_clear <- FALSE

    # The Parameters card: layer, method, contrast, covariates
    # (mod_diff_params.R).
    params <- diff_params_server(input, output, session, current_project, preferred_layer)
    active             <- params$active
    default_contrast   <- params$default_contrast
    continuous_cols    <- params$continuous_cols
    design_mode        <- params$design_mode
    pairing_candidates <- params$pairing_candidates
    levels_            <- params$levels_

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

    # The global test's result, held here because a change of layer or
    # project resets it along with the comparisons (mod_diff_anova.R).
    anova_bundle <- shiny::reactiveVal(NULL)
    anova_error <- shiny::reactiveVal(NULL)
    anova_running <- shiny::reactiveVal(FALSE)

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

    # The settings a result answers. Read with the same fallbacks do_run()
    # uses, so a control that has not rendered yet counts as its default
    # and the comparison below does not flag a result as out of date
    # merely because the page is still drawing.
    settings_now <- shiny::reactive({
      d <- default_contrast()
      mode <- design_mode()
      common <- list(layer = input$layer %||% "", mode = mode,
                     method = input$method %||% "auto",
                     covariates = sort(as.character(input$covariates)))
      if (identical(mode, "continuous")) {
        return(c(common, list(
          col = input$continuous_col %||% continuous_cols()[1L] %||% "",
          model = input$continuous_model %||% "linear")))
      }
      control <- input$control %||% d$control
      cmode <- input$contrast_mode %||% "control"
      custom <- if (identical(cmode, "custom")) {
        x <- trimws(strsplit(input$custom_contrasts %||% "", "\n")[[1L]])
        x[nzchar(x)]
      }
      c(common, list(
        group_col = input$group_col %||% d$group_col %||% "",
        control = control %||% "",
        case = sort(setdiff(input$case %||% d$case, control)),
        paired = input$paired_col %||% "",
        contrast_mode = cmode, custom = custom))
    })
    # What the result on screen was computed with, or, for a result a
    # project brought back, the controls as they stood when it was shown.
    ran_with <- shiny::reactiveVal(NULL)
    # The controls moved since the result was computed: the volcano and
    # the table still answer the old question until Re-run is pressed,
    # and nothing on screen said so.
    settings_changed <- shiny::reactive({
      !is.null(diff_bundle()) && !is.null(ran_with()) &&
        !identical(settings_now(), ran_with())
    })

    # What Re-run does, by groups or along a continuous variable
    # (mod_diff_run.R).
    runner <- diff_run_server(input, output, session, active, default_contrast,
                              design_mode, continuous_cols, settings_now, set_busy,
                              diff_epoch, diff_error, ran_with, diff_bundle)
    do_run         <- runner$do_run
    run_continuous <- runner$run_continuous

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
      proj <- current_project() %||% example_project()
      if (!is.null(b) && identical(omicsCore::bundle_layer(proj, b), input$layer)) return()
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
      tag <- omicsCore::bundle_layer(proj, b)
      if (is.na(tag)) return()
      preferred_layer(tag)
      if (!is.null(input$layer) && !identical(input$layer, tag)) {
        skip_layer_clear <<- TRUE
        shiny::updateSelectInput(session, "layer", selected = tag)
      }
      diff_epoch$bump()
      diff_error(NULL)
      ran_with(shiny::isolate(settings_now()))
      diff_bundle(b)
      a <- proj$bundles$anova
      if (!is.null(a) && identical(omicsCore::bundle_layer(proj, a), tag)) {
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

    # The header, the notices, the stat cards, the volcano and the hit
    # table (mod_diff_results.R).
    diff_results_server(input, output, session, navigate, active, shown_bundle,
                        diff_bundle, diff_error, comparisons, settings_changed,
                        marked, p_col, p_label, fdr_cut_d, fc_cut_d)

    # ---- several contrasts side by side ------------------------------
    contrast_summary_df <- diff_contrasts_server(input, output, session, diff_bundle,
                                                 comparisons, fdr_cut_d, fc_cut_d)

    # ---- global test across all groups -------------------------------
    anova_hits <- diff_anova_server(input, output, session, active, default_contrast,
                                    levels_, anova_bundle, anova_error, anova_running,
                                    anova_epoch, p_col, p_label, fdr_cut_d)

    # Every feature of every comparison, as CSV or Excel (mod_diff_results.R).
    diff_downloads_server(input, output, session, active, diff_bundle, p_col,
                          fdr_cut_d, fc_cut_d)

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
      # The bundle *and* the thresholds it is read at. Enrichment used to
      # take only the bundle and apply its own defaults -- adjusted p at
      # 0.05, no fold-change bound -- so the gene list it enriched was a
      # different set from the hits shown here, and neither view said so.
      # A user with a hundred hits by raw p got an enrichment over however
      # many passed adj.P, which can be none.
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

# Results carry the name of the layer they were computed on: a project
# can hold two proteomics layers, and matching a result to its layer by
# omics type picked the first of them. The demo's layers are not the
# user's, so its results carry none.
with_layer <- function(bundle, active) {
  if (omicsCore::is_analysis_bundle(bundle) && !isTRUE(active$is_demo) &&
      length(active$tag) == 1L) {
    bundle$input_info$layer <- active$tag
  }
  bundle
}
