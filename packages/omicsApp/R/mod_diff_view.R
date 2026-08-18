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
                             invalidate = shiny::reactiveVal(0L)) {
  shiny::moduleServer(id, function(input, output, session) {

    # Active experiment: prefer proteomics, then the first layer
    # of any kind, then NULL (= demo fallback).
    active <- shiny::reactive({
      proj <- current_project()
      if (is.null(proj)) return(list(input = NULL, tag = NULL, is_demo = TRUE))
      exps <- proj$experiments
      if (length(exps) == 0L) return(list(input = NULL, tag = NULL, is_demo = TRUE))
      types <- vapply(exps, function(e) e$omics_type %||% "", character(1))
      idx <- which(types == "proteomics")
      if (length(idx) == 0L) idx <- 1L
      list(
        input   = exps[[idx[1L]]],
        tag     = names(exps)[idx[1L]],
        is_demo = FALSE
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
      inp <- if (a$is_demo) example_proteomics_input() else a$input
      choices <- omicsCore::applicable_diff_methods(inp)
      shiny::selectInput(session$ns("method"), label = NULL,
                         choices = choices, selected = "auto")
    })

    output$method_note <- shiny::renderUI({
      a <- active()
      inp <- if (a$is_demo) example_proteomics_input() else a$input
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

    # Group column dropdown: any meta_df column with >= 2 unique
    # non-NA values (continuous columns like `age` are excluded
    # for the simple "control vs case" UI in this slice).
    output$ui_group_col <- shiny::renderUI({
      a <- active()
      meta <- if (a$is_demo) example_proteomics_input()$meta_df
              else a$input$meta_df
      cands <- names(meta)[vapply(meta, function(col) {
        u <- unique(stats::na.omit(col))
        length(u) >= 2L && !is.numeric(col)
      }, logical(1))]
      if (length(cands) == 0L) cands <- names(meta)
      default <- if ("group" %in% cands) "group" else cands[1L]
      shiny::selectInput(session$ns("group_col"),
                         label    = "Group column",
                         choices  = cands,
                         selected = default)
    })

    # Reactive level set for the chosen group column.
    levels_ <- shiny::reactive({
      a <- active()
      meta <- if (a$is_demo) example_proteomics_input()$meta_df
              else a$input$meta_df
      gc <- input$group_col
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
      htmltools::tagList(
        shiny::selectInput(session$ns("control"),
                           label = "Control", choices = lv,
                           selected = lv[1L]),
        shiny::selectInput(session$ns("case"),
                           label = "Case", choices = lv,
                           selected = lv[min(2L, length(lv))])
      )
    })

    output$ui_covariates <- shiny::renderUI({
      a <- active()
      meta <- if (a$is_demo) example_proteomics_input()$meta_df
              else a$input$meta_df
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

    # The layer this result was computed on has been replaced, so the
    # result is no longer about anything in the project. NULL is the
    # module's own start-up state, so this only rewinds it.
    shiny::observeEvent(invalidate(), {
      diff_bundle(NULL)
      diff_error(NULL)
    }, ignoreInit = TRUE)

    do_run <- function() {
      a <- active()
      if (a$is_demo) {
        diff_error(NULL)
        diff_bundle(example_diff_bundle())
        return(invisible())
      }
      method   <- input$method %||% "auto"
      group_col <- input$group_col
      control  <- input$control
      case     <- input$case
      covariates <- input$covariates
      if (is.null(group_col) || is.null(control) || is.null(case) ||
          identical(control, case)) {
        diff_error("Pick a group column with distinct Control and Case levels.")
        return(invisible())
      }
      run_async(
        function() {
          omicsCore::run_diff(
            input         = a$input,
            method        = method,
            analysis_type = "group",
            group_col     = group_col,
            control_group = control,
            case_group    = case,
            covariates    = if (length(covariates)) covariates else NULL
          )
        },
        on_success = function(bundle) {
          diff_error(NULL)
          diff_bundle(bundle)
        },
        on_error = function(msg) {
          diff_error(msg)
        },
        message = "Running differential analysis..."
      )
    }

    # Initial population: fire as soon as the active() reactive
    # settles. For demo we just emit the fixture; for live data
    # we attempt a default-contrast run.
    shiny::observeEvent(active(), {
      do_run()
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
    fc_cut_d  <- shiny::debounce(shiny::reactive(input$fc_cut  %||% 1),    250)

    marked <- shiny::reactive({
      shiny::req(diff_bundle())
      df <- diff_bundle()$results$diff_result_df
      df$is_significant <- !is.na(df$adj_p_value) &
                           !is.na(df$effect) &
                           df$adj_p_value < fdr_cut_d() &
                           abs(df$effect)  > fc_cut_d()
      df
    })

    output$header <- shiny::renderUI({
      a <- active()
      b <- diff_bundle()
      method <- if (is.null(b)) "\u2014" else b$params$method
      comparison <- if (is.null(b)) "\u2014"
                    else (b$params$comparison %||%
                          sprintf("%s vs %s",
                                  input$case %||% "case",
                                  input$control %||% "control"))
      omics  <- if (is.null(b)) "\u2014"
                else diff_omics_display(b$input_info$omics_type)
      view_header(
        title    = "Differential",
        subtitle = htmltools::tagList(
          omics,
          htmltools::HTML(" &middot; "),
          comparison,
          htmltools::HTML(" &middot; "),
          method,
          htmltools::HTML(" &middot; "),
          htmltools::tags$span(
            class = "muted",
            if (a$is_demo) "demo project (built-in)"
            else sprintf("layer = %s", a$tag)
          )
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
                   else sprintf("effect %+.2f \u00B7 adj.P %.2g",
                                top$effect[1L], top$adj_p_value[1L])
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
          trend  = sprintf("effect > %.2f \u00B7 adj.P < %.3f",
                           fc_cut_d(), fdr_cut_d()),
          accent = "up"
        ),
        stat_card(
          label  = sprintf("Down in %s", input$case %||% "case"),
          value  = down_n,
          trend  = sprintf("effect < -%.2f \u00B7 adj.P < %.3f",
                           fc_cut_d(), fdr_cut_d()),
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
      # Drawn by omicsCore so this figure, the one in an exported
      # report, and the one an exported script reproduces are the same
      # figure. The thresholds come from the sliders, and plot_volcano()
      # re-derives significance from them, so the points re-colour
      # rather than only the dashed rules moving.
      p <- omicsCore::plot_volcano(
        b,
        top_n            = if (isTRUE(input$label_top)) 20L else 0L,
        p_threshold      = fdr_cut_d(),
        effect_threshold = fc_cut_d()
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

    # Expose the bundle so siblings (Enrichment in slice 3E,
    # Report in 3F) can subscribe later. Returning here is
    # harmless when the parent ignores it.
    shiny::reactive(diff_bundle())
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
  missing <- vapply(engines, function(pkg) {
    !requireNamespace(pkg, quietly = TRUE)
  }, logical(1))
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
          shiny::sliderInput(
            ns("fdr_cut"), label = "adj.P cutoff",
            min = 0, max = 0.2, value = 0.05, step = 0.005
          ),
          shiny::sliderInput(
            ns("fc_cut"), label = "|log2FC| cutoff",
            min = 0, max = 4, value = 1, step = 0.05
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
      htmltools::tags$span(class = "card-sub",
                           "drag to zoom \u00B7 double-click to reset")
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
