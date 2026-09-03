#' Pathway enrichment view module
#'
#' Slice 3E: takes a live `diff_bundle` from the Differential
#' view and runs `omicsCore::run_enrichment()` against it
#' (ORA or GSEA, against an MSigDB database). The Re-run button
#' gates the call per slice-3 convention. When clusterProfiler is
#' missing we fall back to `example_enrich_table()` and surface a
#' notice with the install hint.
#'
#' Reference markup: `omicsApp/mockup/index.html:918-964`.
#'
#' @param id Module namespace id.
#' @keywords internal
#' @noRd
enrich_view_ui <- function(id) {
  ns <- shiny::NS(id)

  htmltools::tagList(
    shiny::uiOutput(ns("header")),
    shiny::uiOutput(ns("notices")),
    shiny::uiOutput(ns("input_summary")),
    # Parameters across the top rather than down a 3fr rail: three
    # radio groups and a dropdown do not need a quarter of the width,
    # and the dotplot -- whose y axis carries pathway names -- needed
    # every bit of it. The two result cards then split the full width.
    enrich_params_card(ns),
    htmltools::tags$div(
      class = "row-grid r-7-5",
      enrich_dot_card(ns),
      enrich_hits_card(ns)
    )
  )
}

#' @rdname enrich_view_ui
#' @param diff_bundle Reactive (or reactiveVal) yielding the
#'   live differential `analysis_bundle` from the Diff view, or
#'   `NULL`.
#' @keywords internal
#' @noRd
enrich_view_server <- function(id, diff_bundle = shiny::reactiveVal(NULL),
                               diff_thresholds = shiny::reactive(list(
                                 p_cutoff = 0.05, p_preference = "adjusted",
                                 effect_cutoff = NULL)),
                               invalidate = shiny::reactiveVal(0L)) {
  shiny::moduleServer(id, function(input, output, session) {

    have_cp <- requireNamespace("clusterProfiler", quietly = TRUE)

    enrich_bundle <- shiny::reactiveVal(NULL)
    enrich_error  <- shiny::reactiveVal(NULL)
    is_demo       <- shiny::reactiveVal(TRUE)

    # The layer the upstream diff was computed on has been replaced, so
    # this enrichment is no longer about anything in the project. Back
    # to the module's own start-up state.
    shiny::observeEvent(invalidate(), {
      enrich_bundle(NULL)
      enrich_error(NULL)
      is_demo(TRUE)
    }, ignoreInit = TRUE)

    do_run <- function() {
      db_arg <- input$database %||% "hallmark"
      type   <- input$type %||% "ora"
      dir_   <- input$direction %||% "both"
      thr    <- diff_thresholds()
      bundle <- diff_bundle()
      if (is.null(bundle) || !have_cp) {
        # Demo fallback: synthetic table re-shaped to match the
        # standardized enrich schema.
        enrich_error(if (!have_cp) {
          paste0("clusterProfiler is not installed; showing the demo ",
                 "fixture. Install with `omicsCore::install_optional",
                 "('enrichment')`.")
        } else NULL)
        is_demo(TRUE)
        enrich_bundle(NULL)
        return(invisible())
      }
      run_async(
        # Detached for the same reason as the Differential view: a
        # closure defined here carries this module's whole scope to the
        # worker, and an enrichment runs on a bundle that is already
        # large.
        detached_call(
          function() {
            # ORA takes the thresholds; GSEA does not.
            #
            # GSEA is a ranked-list method: it reads every feature,
            # ordered by effect, and pre-filtering the list is not a
            # stricter version of the analysis but a different and
            # invalid one. run_enrichment() also spends p_cutoff on
            # bounding the *pathway* table, so handing it a threshold
            # chosen for features would quietly narrow which pathways
            # are reported -- two different tests on two different
            # objects, sharing one number.
            if (identical(type, "ora")) {
              omicsCore::run_enrichment(
                diff_bundle   = bundle,
                type          = type,
                database      = db_arg,
                direction     = dir_,
                p_cutoff      = thr$p_cutoff,
                p_preference  = thr$p_preference,
                effect_cutoff = thr$effect_cutoff
              )
            } else {
              omicsCore::run_enrichment(
                diff_bundle = bundle,
                type        = type,
                database    = db_arg,
                direction   = dir_
              )
            }
          },
          bundle = bundle, type = type, db_arg = db_arg, dir_ = dir_,
          thr = thr
        ),
        on_success = function(result) {
          enrich_error(NULL)
          is_demo(FALSE)
          enrich_bundle(result)
        },
        on_error = function(msg) {
          enrich_error(msg)
          is_demo(TRUE)
          enrich_bundle(NULL)
        },
        message = "Running pathway enrichment..."
      )
    }

    # Auto-run once on first diff_bundle change so the view lands
    # populated; subsequent updates require Re-run.
    shiny::observeEvent(diff_bundle(), {
      do_run()
    }, ignoreInit = TRUE)
    shiny::observeEvent(input$rerun, do_run())

    # The table powering both the dot card and the hits card. In
    # demo mode this is the static fixture; in live mode it's
    # the enrichment bundle's standardized result data frame. We
    # guard the live path with `req()` so a transient NULL during
    # an async flush doesn't crash the reactive graph (downstream
    # renderers also `req()` on `nrow(df) > 0`).
    table_data <- shiny::reactive({
      if (isTRUE(is_demo())) return(example_enrich_table())
      b <- enrich_bundle()
      shiny::req(b)
      b$results$enrich_result_df
    })

    output$header <- shiny::renderUI({
      b <- enrich_bundle()
      demo <- is_demo()
      omics <- if (is.null(b)) "Proteomics"
               else enrich_omics_display(b$input_info$omics_type)
      method <- if (is.null(b)) "MSigDB Hallmark"
                else sprintf("%s \u00B7 %s",
                             toupper(b$params$type %||% "ora"),
                             b$params$database %||% "hallmark")
      view_header(
        title    = "Pathway enrichment",
        subtitle = htmltools::tagList(
          omics,
          htmltools::HTML(" &middot; "),
          method,
          htmltools::HTML(" &middot; "),
          htmltools::tags$span(
            class = "muted",
            if (demo) "demo fixture (built-in)"
            else "live result"
          )
        )
      )
    })

    # The size of the gene list this enrichment is over.
    #
    # An empty result reads identically whether no feature met the
    # threshold or three thousand did and none of the sets were
    # enriched. Reported as "100+ differential genes but every database
    # comes back with nothing", which is a question this answers before
    # anyone has to ask it.
    selected_features <- shiny::reactive({
      b <- diff_bundle()
      if (is.null(b)) return(NULL)
      df <- b$results$diff_result_df
      if (is.null(df)) return(NULL)
      thr <- diff_thresholds()
      # GSEA ranks every feature it can place, so there is no selection
      # to report -- only how many made it into the ranking. Saying "N
      # of M selected" here would describe a step the method does not
      # take.
      if (!identical(input$type %||% "ora", "ora")) {
        rankable <- !is.na(df$effect) & !is.na(df$feature_symbol)
        return(list(gsea = TRUE, n = sum(rankable), total = nrow(df),
                    mapped = sum(rankable & df$feature_symbol != df$feature_id)))
      }
      pcol <- if (identical(thr$p_preference, "raw")) "p_value" else "adj_p_value"
      pv <- df[[pcol]]
      keep <- !is.na(pv) & pv < (thr$p_cutoff %||% 0.05)
      if (!is.null(thr$effect_cutoff) && is.finite(thr$effect_cutoff)) {
        keep <- keep & !is.na(df$effect) & abs(df$effect) > thr$effect_cutoff
      }
      list(gsea = FALSE, n = sum(keep), total = nrow(df), p_col = pcol,
           p_cutoff = thr$p_cutoff %||% 0.05,
           effect_cutoff = thr$effect_cutoff,
           mapped = sum(keep & !is.na(df$feature_symbol) &
                          df$feature_symbol != df$feature_id))
    })

    output$input_summary <- shiny::renderUI({
      s <- selected_features()
      if (is.null(s)) return(NULL)
      if (isTRUE(s$gsea)) {
        return(htmltools::tags$div(
          class = "muted", style = "font-size:12.5px;margin:2px 0 10px",
          htmltools::tags$strong(sprintf("%d of %d features ranked",
                                         s$n, s$total)),
          " \u2014 GSEA reads the whole list, so the Differential ",
          "thresholds do not apply here.",
          if (s$mapped == 0L) {
            htmltools::tags$span(
              style = "color:var(--warn)",
              " None carry a gene symbol, so no set will match."
            )
          }
        ))
      }
      bits <- sprintf("%s < %s", s$p_col, format(s$p_cutoff))
      if (!is.null(s$effect_cutoff) && is.finite(s$effect_cutoff)) {
        bits <- paste0(bits, sprintf(", |effect| > %.3f", s$effect_cutoff))
      }
      htmltools::tags$div(
        class = "muted", style = "font-size:12.5px;margin:2px 0 10px",
        htmltools::tags$strong(sprintf("%d of %d features", s$n, s$total)),
        sprintf(" selected at %s", bits),
        if (s$n == 0L) {
          htmltools::tags$span(
            style = "color:var(--warn)",
            " \u2014 nothing to enrich. Loosen the thresholds in the ",
            "Differential view."
          )
        } else if (s$mapped == 0L) {
          htmltools::tags$span(
            style = "color:var(--warn)",
            " \u2014 none carry a gene symbol, so no set will match. ",
            "Check the Gene symbol line on the Import page."
          )
        }
      )
    })

    output$notices <- shiny::renderUI({
      tagged <- htmltools::tagList()
      err <- enrich_error()
      if (!is.null(err)) {
        tagged <- htmltools::tagAppendChild(
          tagged,
          notice(title  = if (!have_cp) "clusterProfiler unavailable"
                          else "run_enrichment failed",
                 detail = err,
                 kind   = "warn")
        )
      }
      if (is.null(diff_bundle()) && have_cp) {
        tagged <- htmltools::tagAppendChild(
          tagged,
          notice(
            title  = "Run a differential analysis first",
            detail = paste0("This view enriches the top hits from the ",
                            "Differential view. Showing the demo fixture ",
                            "until a real bundle arrives."),
            kind   = "info"
          )
        )
      }
      tagged
    })

    # Both paths hand a bundle to the same dispatcher. When the demo
    # had its own drawing code it was free to drift from the live view,
    # and a demo that no longer resembles the product is worse than no
    # demo.
    plot_bundle <- shiny::reactive({
      if (isTRUE(is_demo())) example_enrich_bundle() else enrich_bundle()
    })

    # The threshold the panel is read at. numericInput reports NA while
    # the box is mid-edit, and NA would empty the panel with nothing to
    # say why.
    show_p <- shiny::reactive(input$show_p %||% "adjusted")
    show_cutoff <- shiny::reactive({
      v <- input$show_cutoff
      if (is.null(v) || !is.finite(v) || v <= 0 || v > 1) 0.05 else v
    })

    output$dot <- shiny::renderPlot({
      b <- plot_bundle()
      shiny::req(b)
      omicsCore::plot_enrichment(b, view = "dot", top_n = 12L,
                                 p_preference = show_p(),
                                 p_cutoff = show_cutoff())
    })

    output$hits <- DT::renderDT({
      df <- table_data()
      shiny::req(nrow(df) > 0L)
      df <- omicsCore::filter_enrich_results(
        df, p_cutoff = show_cutoff(), p_preference = show_p())
      shiny::req(nrow(df) > 0L)
      pcol <- if (identical(show_p(), "raw")) "p_value" else "adj_p_value"
      df <- df[order(df[[pcol]]), , drop = FALSE]
      out <- data.frame(
        Pathway   = df$pathway_name,
        NES       = sprintf("%+.2f", df$effect),
        # Named for the column it holds, so the table cannot say adj.P
        # over raw values.
        P         = signif(df[[pcol]], 3),
        Direction = df$direction,
        Overlap   = sprintf("%d/%d",
                            df$overlap_size %||% NA_integer_,
                            df$gene_set_size %||% NA_integer_),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      names(out)[names(out) == "P"] <-
        if (identical(show_p(), "raw")) "p" else "adj.P"
      DT::datatable(
        out,
        rownames  = FALSE,
        selection = "single",
        options   = list(
          pageLength = 10,
          dom        = "ftip",
          scrollX    = TRUE,
          columnDefs = list(list(className = "dt-right",
                                 targets = c(1, 2)))
        )
      )
    }, server = TRUE)

    # Expose the bundle for slice 3F (report).
    shiny::reactive(enrich_bundle())
  })
}

# ---- internal helpers ------------------------------------------------

utils::globalVariables(".data")

`%||%` <- function(a, b) if (is.null(a)) b else a

enrich_omics_display <- function(t) {
  switch(t %||% "",
         proteomics = "Proteomics",
         rnaseq     = "RNA-seq",
         "\u2014")
}

enrich_params_card <- function(ns) {
  bslib::card(
    # A dropdown opened from this card renders *inside* it, and a card
    # clips its overflow -- so the Database menu was cut off at the card
    # edge. Only an issue since the controls moved to a top row: down a
    # tall rail there was always card below the menu.
    class = "param-row-card",
    bslib::card_header(
      htmltools::tags$h3(class = "card-title", "Parameters"),
      htmltools::tags$span(class = "card-sub",
                           "test type, database, direction")
    ),
    bslib::card_body(
      htmltools::tags$div(
        class = "param-row",
        param_group(
          "Test",
          shiny::radioButtons(
            ns("type"), label = NULL,
            choices  = c("ORA" = "ora", "GSEA" = "gsea"),
            selected = "ora", inline = TRUE
          )
        ),
        param_group(
          "Database",
          shiny::selectInput(
            ns("database"), label = NULL,
            choices  = c("hallmark", "kegg", "reactome",
                         "go_bp", "go_mf", "go_cc",
                         "wikipathways"),
            selected = "hallmark"
          )
        ),
        param_group(
          "Direction",
          shiny::radioButtons(
            ns("direction"), label = NULL,
            choices  = c("both", "up", "down"),
            selected = "both", inline = TRUE
          )
        ),
        param_group(
          "Significance",
          # A display threshold, like the Differential view's sliders:
          # applied to the result that is already computed, so switching
          # it costs nothing and re-runs nothing. GSEA over 2695 GO sets
          # keeps 9 at adjusted p and 131 at raw; which of those is the
          # right answer is the reader's call, not ours.
          shiny::radioButtons(
            ns("show_p"), label = NULL,
            choices  = c("adjusted p" = "adjusted", "raw p" = "raw"),
            selected = "adjusted", inline = TRUE
          ),
          shiny::numericInput(
            ns("show_cutoff"), label = NULL,
            value = 0.05, min = 0, max = 1, step = 0.01
          )
        ),
        htmltools::tags$div(
          class = "param-action",
          shiny::actionButton(
            ns("rerun"), "Re-run",
            class = "btn btn-primary"
          )
        )
      )
    )
  )
}

enrich_dot_card <- function(ns) {
  bslib::card(
    bslib::card_header(
      htmltools::tags$h3(class = "card-title", "Pathway dotplot"),
      htmltools::tags$span(class = "card-sub",
                           "top 12 by adj.P \u00B7 size = overlap")
    ),
    bslib::card_body(
      shiny::plotOutput(ns("dot"), height = "420px")
    )
  )
}

enrich_hits_card <- function(ns) {
  bslib::card(
    bslib::card_header(
      htmltools::tags$h3(class = "card-title", "Enriched sets"),
      htmltools::tags$span(class = "card-sub",
                           "ranked by adj.P")
    ),
    bslib::card_body(
      DT::DTOutput(ns("hits"))
    )
  )
}
