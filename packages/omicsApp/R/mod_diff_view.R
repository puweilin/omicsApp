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
    active <- shiny::reactive({
      proj <- current_project()
      is_demo <- is.null(proj)
      if (is_demo) proj <- example_project()
      exps <- proj$experiments
      if (length(exps) == 0L) {
        return(list(input = NULL, tag = NULL, is_demo = TRUE))
      }
      want <- input$layer
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
      sel <- shiny::isolate(input$layer)
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
      shiny::selectInput(session$ns("method"), label = NULL,
                         choices = choices, selected = "auto")
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
      cands <- grouping_candidates(meta)
      if (length(cands) == 0L) cands <- names(meta)
      gc <- cands[1L]
      lv <- sort(unique(as.character(stats::na.omit(meta[[gc]]))))
      ctrl <- default_control_level(lv)
      list(
        group_col  = gc,
        control    = ctrl,
        # Every other group against the control: with one control and
        # several treatments that is the design people mean, and with two
        # groups it is the one contrast there is.
        case       = if (length(lv) >= 2L) setdiff(lv, ctrl) else NULL,
        levels     = lv,
        candidates = cands
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
      ctrl <- keep(input$control, default_control_level(lv))[1L]
      htmltools::tagList(
        shiny::selectInput(session$ns("control"),
                           label = "Control (reference)", choices = lv,
                           selected = ctrl),
        shiny::selectizeInput(
          session$ns("case"),
          label = if (length(lv) > 2L) "Compare against control (one or more)" else "Case",
          choices = lv, multiple = TRUE,
          selected = keep(input$case, setdiff(lv, ctrl)),
          options = list(plugins = list("remove_button"))),
        if (length(lv) > 2L) {
          htmltools::tags$div(
            class = "muted", style = "font-size:11.5px;margin-top:-6px",
            paste("Several groups are fitted in one model, so each comparison",
                  "with the control borrows strength from all samples."))
        }
      )
    })

    output$ui_covariates <- shiny::renderUI({
      a <- active()
      meta <- a$input$meta_df
      gc <- input$group_col %||% ""
      cands <- setdiff(names(meta), c(gc, "sample_id"))
      shiny::selectizeInput(
        session$ns("covariates"),
        label    = NULL,
        choices  = cands,
        multiple = TRUE,
        selected = NULL,
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
      if (is.null(sel) || !sel %in% cmp) sel <- cmp[[1L]]
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
      diff_bundle(NULL)
      diff_error(NULL)
    }, ignoreInit = TRUE)

    do_run <- function() {
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
      case       <- setdiff(case, control)
      covariates <- input$covariates
      if (is.null(group_col) || is.null(control) || !length(case)) {
        diff_error("Pick a group column, a control group, and at least one group distinct from the control to compare with it.")
        return(invisible())
      }
      set_busy(TRUE)
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
          set_busy(FALSE)
          diff_error(NULL)
          diff_bundle(bundle)
        },
        on_error = function(msg) {
          set_busy(FALSE)
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
      shiny::req(shown_bundle())
      df <- shown_bundle()$results$diff_result_df
      pv <- df[[p_col()]]
      df$is_significant <- !is.na(pv) &
                           !is.na(df$effect) &
                           pv < fdr_cut_d() &
                           abs(df$effect) > fc_cut_d()
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
          notice(title  = "The differential analysis could not run",
                 detail = diff_error_hint(err),
                 kind   = "error",
                 technical = err)
        )
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
      case_lbl <- b$params$case_group %||% "case"
      if (length(case_lbl) > 1L) case_lbl <- "case"
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
          label  = sprintf("Up in %s", case_lbl),
          value  = up_n,
          trend  = sprintf("effect > %.2f \u00B7 %s < %.3f",
                           fc_cut_d(), p_label(), fdr_cut_d()),
          accent = "up"
        ),
        stat_card(
          label  = sprintf("Down in %s", case_lbl),
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
        p         = signif(sig[[p_col()]], 3),
        Direction = ifelse(sig$effect > 0, "up", "down"),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      # The column is named for the p-value the mask was read from.
      names(out)[3] <- p_label()
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
                              height = paste0(90 + 42 * length(comparisons()), "px")),
            DT::DTOutput(session$ns("contrast_table"))
          )
        )
      )
    })

    contrast_summary_df <- shiny::reactive({
      b <- diff_bundle()
      shiny::req(omicsCore::is_analysis_bundle(b), length(comparisons()) > 1L)
      omicsCore::summarize_diff_contrasts(
        b, p_cutoff = fdr_cut_d(),
        p_preference = if (identical(input$p_kind %||% "adj", "raw")) "raw" else "adjusted",
        effect_cutoff = fc_cut_d())
    })

    output$contrast_plot <- shiny::renderPlot({
      b <- diff_bundle()
      shiny::req(omicsCore::is_analysis_bundle(b), length(comparisons()) > 1L)
      omicsCore::plot_diff_contrasts(
        b, p_cutoff = fdr_cut_d(),
        p_preference = if (identical(input$p_kind %||% "adj", "raw")) "raw" else "adjusted",
        effect_cutoff = fc_cut_d())
    })

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

    output$run_button <- shiny::renderUI({
      shiny::actionButton(
        session$ns("rerun"),
        if (is.null(diff_bundle())) "Run analysis" else "Re-run",
        icon = shiny::icon("play"),
        class = "btn btn-primary", style = "width:100%")
    })

    if (is.function(navigate)) {
      shiny::observeEvent(input$go_enrich, navigate("enrich"))
    }

    list(
      bundle = shown_bundle,
      # Every contrast of the last run, for anything that wants them all.
      all_bundle = shiny::reactive(diff_bundle()),
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
          shiny::uiOutput(ns("ui_group_col")),
          shiny::uiOutput(ns("ui_contrast")),
          shiny::uiOutput(ns("ui_comparison"))
        ),
        param_group(
          "Covariates",
          help = paste("Variables to adjust for (age, sex, batch). The group effect",
                       "is then estimated holding them constant. Leave empty for",
                       "an unadjusted comparison."),
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
                           "ranked by |effect| within current thresholds")
    ),
    bslib::card_body(
      DT::DTOutput(ns("hits"))
    )
  )
}

# The level a control group is usually called, so the default reference
# is the control rather than whichever label sorts first ("DMSO" before
# "Drug", but "Treated" after "Control" only by luck).
CONTROL_LEVEL_PATTERN <- paste0(
  "^(ctrl|control|controls|con|dmso|vehicle|veh|mock|sham|untreated|",
  "placebo|baseline|wt|wild[ _-]?type|normal|healthy|nc|neg|negative|",
  "blank|t0|0h|day0|d0)$")

default_control_level <- function(levels) {
  if (!length(levels)) return(NULL)
  hit <- grepl(CONTROL_LEVEL_PATTERN, trimws(tolower(levels)))
  if (!any(hit)) {
    hit <- grepl("(^|[ _-])(ctrl|control|vehicle|dmso|wt)([ _-]|$)",
                 tolower(levels))
  }
  levels[if (any(hit)) which(hit)[1L] else 1L]
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
