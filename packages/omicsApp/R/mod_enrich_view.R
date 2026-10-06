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
    ),
    # Only when the differential run holds several comparisons.
    shiny::uiOutput(ns("compare_card"))
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
                               diff_layer = shiny::reactive(NULL),
                               invalidate = shiny::reactiveVal(0L),
                               navigate = NULL,
                               diff_all = shiny::reactive(NULL),
                               current_project = shiny::reactive(NULL)) {
  shiny::moduleServer(id, function(input, output, session) {

    # ---- every comparison side by side --------------------------------
    # The main result is the comparison on screen in the Differential
    # view. With several comparisons the other question is which pathways
    # they share: each is enriched with the same settings, and the dot
    # plot sets them next to each other.
    compare_bundle <- shiny::reactiveVal(NULL)
    compare_error <- shiny::reactiveVal(NULL)
    all_comparisons <- shiny::reactive({
      b <- diff_all()
      if (!omicsCore::is_analysis_bundle(b)) return(character(0))
      omicsCore::diff_comparisons(b)
    })
    # A result that arrives after the comparisons it was asked about
    # were replaced is dropped (see run_epoch()).
    compare_epoch <- run_epoch()
    compare_running <- shiny::reactiveVal(FALSE)
    shiny::observeEvent(diff_all(), {
      compare_epoch$bump()
      compare_bundle(NULL)
      compare_error(NULL)
    }, ignoreNULL = FALSE, ignoreInit = TRUE)

    output$compare_card <- shiny::renderUI({
      if (length(all_comparisons()) < 2L || !have_cp) return(NULL)
      b <- compare_bundle()
      bslib::card(
        bslib::card_header(
          htmltools::tags$h3(class = "card-title", "Across comparisons"),
          htmltools::tags$span(
            class = "card-sub",
            sprintf("%d comparisons \u00B7 same test, database and thresholds",
                    length(all_comparisons())))
        ),
        bslib::card_body(
          htmltools::tags$div(
            style = "display:flex;gap:12px;align-items:center;flex-wrap:wrap",
            disabled_if(shiny::actionButton(session$ns("run_compare"),
                                            if (is.null(b)) "Enrich every comparison"
                                            else "Re-run for every comparison",
                                            class = "btn btn-sm btn-outline-primary"),
                        shiny::isolate(compare_running())),
            htmltools::tags$span(class = "muted", style = "font-size:12px",
                                 paste("Shared pathways sort to the top; an",
                                       "empty cell means the pathway was not",
                                       "found in that comparison."))
          ),
          if (!is.null(compare_error())) {
            notice("The comparison could not be computed", kind = "error",
                   technical = compare_error())
          },
          if (!is.null(b) && length(b$warnings)) {
            notice("Some comparisons found nothing",
                   detail = paste(b$warnings, collapse = " "), kind = "info")
          },
          if (!is.null(b)) {
            shiny::plotOutput(session$ns("compare_plot"), height = "480px")
          }
        )
      )
    })

    shiny::observeEvent(input$run_compare, {
      b <- diff_all()
      shiny::req(omicsCore::is_analysis_bundle(b), have_cp)
      thr <- diff_thresholds()
      type <- input$type %||% "ora"
      args <- list(type = type, database = input$database %||% "hallmark",
                   direction = run_direction(b),
                   organism = organism())
      if (identical(type, "ora")) {
        args <- c(args, list(p_cutoff = thr$p_cutoff,
                             p_preference = thr$p_preference,
                             effect_cutoff = thr$effect_cutoff))
      }
      my_run <- compare_epoch$start()
      set_button_busy("run_compare", TRUE, compare_running)
      run_async(
        detached_call(
          function() do.call(omicsCore::compare_enrichment,
                             c(list(diff_bundle = bundle), args)),
          bundle = b, args = args
        ),
        on_success = function(res) {
          if (compare_epoch$is_last_started(my_run)) set_button_busy("run_compare", FALSE, compare_running)
          if (!compare_epoch$is_current(my_run)) return(invisible())
          compare_error(NULL)
          compare_bundle(res)
        },
        on_error = function(msg) {
          if (compare_epoch$is_last_started(my_run)) set_button_busy("run_compare", FALSE, compare_running)
          if (compare_epoch$is_current(my_run)) compare_error(msg)
        },
        message = "Enriching every comparison..."
      )
    })

    output$compare_plot <- shiny::renderPlot(res = PLOT_RES, alt = "Pathways enriched in each comparison, side by side", fit_to_width("compare_plot", {
      b <- compare_bundle()
      shiny::req(b)
      omicsCore::plot_enrichment_comparison(
        b, p_preference = input$show_p %||% "adjusted")
    }))

    have_cp <- has_pkg("clusterProfiler")

    # The species the gene sets come from. Suggested from the symbols --
    # mouse genes are written Trp53, human TP53, zebrafish tp53 -- because
    # the default human sets matched no mouse gene, and the result was an
    # empty table with no word of why. Every species msigdbr carries the
    # MSigDB sets to (through orthologs) is offered; mouse and rat write
    # their genes alike, so a rat table is suggested mouse and switched by
    # hand.
    suggested_organism <- shiny::reactive({
      b <- diff_bundle()
      if (is.null(b)) return("Hs")
      guess_organism(b$results$diff_result_df$feature_symbol)
    })
    organism <- shiny::reactive(input$organism %||% suggested_organism())
    output$ui_organism <- shiny::renderUI({
      sel <- shiny::isolate(input$organism) %||% suggested_organism()
      shiny::selectInput(session$ns("organism"), label = "Species",
                         choices = species_choices(), selected = sel)
    })
    shiny::observeEvent(suggested_organism(), {
      shiny::updateSelectInput(session, "organism", selected = suggested_organism())
    }, ignoreInit = TRUE)

    enrich_bundle <- shiny::reactiveVal(NULL)
    enrich_error  <- shiny::reactiveVal(NULL)
    is_demo       <- shiny::reactiveVal(TRUE)
    # No demo pathways once the user has a project of their own.
    shiny::observeEvent(current_project(), {
      if (is.null(enrich_bundle())) is_demo(is.null(current_project()) || !have_cp)
    }, ignoreNULL = FALSE)

    # The layer the upstream diff was computed on has been replaced, so
    # this enrichment is no longer about anything in the project. Back
    # to the module's own start-up state.
    enrich_epoch <- run_epoch()
    shiny::observeEvent(invalidate(), {
      enrich_epoch$bump()
      compare_epoch$bump()
      enrich_bundle(NULL)
      enrich_error(NULL)
      is_demo(is.null(current_project()))
      compare_bundle(NULL)
    }, ignoreInit = TRUE)

    do_run <- function() {
      my_run <- enrich_epoch$start()
      db_arg <- input$database %||% "hallmark"
      type   <- input$type %||% "ora"
      thr    <- diff_thresholds()
      bundle <- diff_bundle()
      dir_   <- run_direction(bundle)
      if (is.null(bundle) || !have_cp) {
        # Demo fallback: synthetic table re-shaped to match the
        # standardized enrich schema.
        enrich_error(if (!have_cp) {
          paste0("clusterProfiler is not installed; showing the demo ",
                 "fixture. Install with `omicsCore::install_optional",
                 "('enrichment')`.")
        } else NULL)
        # The demo only when there is no project of the user's: its
        # pathways under a restored project read as the user's result.
        is_demo(is.null(current_project()) || !have_cp)
        enrich_bundle(NULL)
        return(invisible())
      }
      set_button_busy("rerun", TRUE)
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
                organism      = org,
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
                organism    = org,
                direction   = dir_
              )
            }
          },
          bundle = bundle, type = type, db_arg = db_arg, dir_ = dir_,
          thr = thr, org = organism()
        ),
        on_success = function(result) {
          if (enrich_epoch$is_last_started(my_run)) set_button_busy("rerun", FALSE)
          if (!enrich_epoch$is_current(my_run)) return(invisible())
          enrich_error(NULL)
          is_demo(FALSE)
          enrich_bundle(result)
        },
        on_error = function(msg) {
          if (enrich_epoch$is_last_started(my_run)) set_button_busy("rerun", FALSE)
          if (!enrich_epoch$is_current(my_run)) return(invisible())
          # The error, and nothing under it: the demo's pathways beneath
          # a failure notice read as the user's result.
          enrich_error(msg)
          is_demo(FALSE)
          enrich_bundle(NULL)
        },
        message = "Running pathway enrichment..."
      )
    }

    # Auto-run once on first diff_bundle change so the view lands
    # populated; subsequent updates require Re-run.
    # Also when the upstream result goes away (a layer change clears it):
    # with ignoreNULL the old enrichment stayed on screen, now labelled
    # with the new layer, and went into the project without its diff.
    shiny::observeEvent(diff_bundle(), {
      # A project opened or restored brings its enrichment with it: shown
      # as it was saved, not recomputed at whatever the controls say.
      saved <- current_project()$bundles$enrich
      b <- diff_bundle()
      if (!is.null(saved) && !is.null(b) && is.null(enrich_bundle()) &&
          identical(saved$params$comparison, b$params$comparison) &&
          identical(saved$input_info$omics_type, b$input_info$omics_type)) {
        enrich_epoch$bump()
        restore_controls(saved$params)
        enrich_error(NULL)
        is_demo(FALSE)
        enrich_bundle(saved)
        sc <- current_project()$bundles$enrich_compare
        if (!is.null(sc) && is.null(compare_bundle())) compare_bundle(sc)
        return(invisible())
      }
      do_run()
    }, ignoreInit = TRUE, ignoreNULL = FALSE)

    restore_controls <- function(params) {
      if (!is.null(params$type)) shiny::updateRadioButtons(session, "type", selected = params$type)
      if (!is.null(params$database)) shiny::updateSelectInput(session, "database",
                                                              selected = params$database[[1L]])
      if (!is.null(params$direction)) shiny::updateRadioButtons(session, "direction",
                                                                selected = params$direction)
      org <- species_code(params$organism)
      if (!is.null(org)) shiny::updateSelectInput(session, "organism", selected = org)
    }

    # The gene list(s) to test, as the Direction control says -- except
    # for a result with no direction (a spline fit), which has one list
    # only: it is enriched pooled, whatever the control is set to.
    run_direction <- function(bundle) {
      if (diff_undirected(bundle)) return("both")
      input$direction %||% "separate"
    }
    shiny::observeEvent(input$rerun, do_run())
    if (is.function(navigate)) {
      shiny::observeEvent(input$go_next, navigate("integration"))
    }

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
      # Which layer, stated rather than offered. Enrichment runs on a
      # differential result, so its layer is whichever one Differential
      # was run on -- a picker here could be set to disagree with that,
      # and the panel would then be labelled rnaseq while showing
      # proteomics pathways. To change it, change it there.
      layer <- diff_layer_tag()
      view_header(
        title    = "Pathway enrichment",
        actions  = if (is.function(navigate) && !demo) {
          shiny::actionButton(session$ns("go_next"), "Next: Integration \u2192",
                              class = "btn btn-ghost")
        },
        subtitle = htmltools::tagList(
          omics,
          htmltools::HTML(" &middot; "),
          method,
          if (!is.null(layer)) {
            htmltools::tagList(
              htmltools::HTML(" &middot; "),
              htmltools::tags$span(
                class = "muted",
                sprintf("from the %s differential", layer))
            )
          },
          htmltools::HTML(" &middot; "),
          htmltools::tags$span(
            class = "muted",
            if (demo) "demo data (built-in)"
            else "live result"
          )
        )
      )
    })

    # The tag of the layer the upstream differential ran on, handed over
    # by that view. NULL for the demo.
    diff_layer_tag <- shiny::reactive({
      if (isTRUE(is_demo())) return(NULL)
      tag <- diff_layer()
      if (is.null(tag) || !nzchar(tag)) NULL else tag
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
      # Once an ORA has run, the thresholds it ran at: the live ones may
      # have moved since, and the line describes the result below it.
      eb <- enrich_bundle()
      if (!is.null(eb) && identical(eb$params$type, "ora") &&
          identical(input$type %||% "ora", "ora")) {
        thr <- list(p_cutoff = eb$params$p_cutoff,
                    p_preference = eb$params$p_preference,
                    effect_cutoff = eb$params$effect_cutoff)
      }
      # GSEA ranks every feature it can place, so there is no selection
      # to report -- only how many made it into the ranking. Saying "N
      # of M selected" here would describe a step the method does not
      # take.
      if (!identical(input$type %||% "ora", "ora")) {
        rankable <- !is.na(df$effect) & !is.na(df$feature_symbol)
        return(list(gsea = TRUE, n = sum(rankable), total = nrow(df),
                    mapped = sum(rankable & is_gene_symbol(df$feature_symbol))))
      }
      pcol <- if (identical(thr$p_preference, "raw")) "p_value" else "adj_p_value"
      pv <- df[[pcol]]
      keep <- !is.na(pv) & pv < (thr$p_cutoff %||% 0.05)
      if (!is.null(thr$effect_cutoff) && is.finite(thr$effect_cutoff)) {
        keep <- keep & !is.na(df$effect) & abs(df$effect) >= thr$effect_cutoff
      }
      list(gsea = FALSE, n = sum(keep), total = nrow(df), p_col = pcol,
           p_cutoff = thr$p_cutoff %||% 0.05,
           effect_cutoff = thr$effect_cutoff,
           # A symbol is a symbol whether or not it is also the id: a
           # matrix keyed by gene symbols was told "none carry a gene
           # symbol" beside the pathways it had just found.
           mapped = sum(keep & is_gene_symbol(df$feature_symbol)))
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
      p_name <- switch(s$p_col, adj_p_value = "adjusted p", p_value = "p", s$p_col)
      bits <- sprintf("%s < %s", p_name, format(s$p_cutoff))
      if (!is.null(s$effect_cutoff) && is.finite(s$effect_cutoff)) {
        bits <- paste0(bits, sprintf(", |%s| \u2265 %.3f",
                                     omicsCore::effect_label(diff_bundle()), s$effect_cutoff))
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
          if (identical(err, CANCELLED_MESSAGE)) {
            notice("Cancelled", "Press Re-run to start again.", kind = "info")
          } else if (!have_cp) {
            notice(title = "clusterProfiler unavailable", detail = err,
                   kind = "warn")
          } else {
            notice(title  = "The enrichment could not be computed",
                   detail = if (grepl("No significant|no features|no genes", err, ignore.case = TRUE))
                     "No differential hits pass the current thresholds; loosen them in the Differential view or try GSEA."
                   else "See the technical details below.",
                   kind   = "error",
                   technical = err)
          }
        )
      }
      if (is.null(diff_bundle()) && have_cp) {
        tagged <- htmltools::tagAppendChild(
          tagged,
          notice(
            title  = "Run a differential analysis first",
            detail = if (is.null(current_project()))
              paste0("This view enriches the top hits from the ",
                     "Differential view. Showing the demo fixture ",
                     "until a real result arrives.")
            else "This view enriches the hits of the Differential view; run it there and the pathways appear here.",
            kind   = "info"
          )
        )
      }
      # Nothing found: say the likely reason instead of an empty table.
      eb <- enrich_bundle()
      if (!is.null(eb) && !isTRUE(is_demo()) &&
          nrow(eb$results$enrich_result_df %||% data.frame()) == 0L) {
        guess <- guess_organism(diff_bundle()$results$diff_result_df$feature_symbol)
        used <- species_code(eb$params$organism) %||% "Hs"
        tagged <- htmltools::tagAppendChild(tagged, notice(
          title = "No pathway was found",
          detail = if (!identical(guess, used))
            sprintf("The gene names look like %s genes (e.g. %s), but %s gene sets were used. Switch Species and re-run.",
                    tolower(species_label(guess)), species_example(guess),
                    tolower(species_label(used)))
          else "None of the sets passed the threshold. Try GSEA, another database, or looser Differential thresholds.",
          kind = "warn"))
      }
      # What the run itself had to say -- the gene names matched to the
      # gene sets only once upper and lower case were ignored, say.
      if (!is.null(eb) && !isTRUE(is_demo()) && length(eb$warnings)) {
        for (w in eb$warnings) {
          tagged <- htmltools::tagAppendChild(tagged, notice(
            "Note on this result", detail = w, kind = "info"))
        }
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
    show_p <- shiny::reactive({
      v <- input$show_p
      if (is.null(v) || !is.character(v) || length(v) != 1L ||
          !v %in% c("adjusted", "raw", "qvalue")) "adjusted" else v
    })
    show_cutoff <- shiny::reactive({
      v <- input$show_cutoff
      if (is.null(v) || !is.finite(v) || v <= 0 || v > 1) 0.05 else v
    })

    output$dot <- shiny::renderPlot(res = PLOT_RES, alt = "Dot plot of the most enriched pathways", fit_to_width("dot", {
      b <- plot_bundle()
      shiny::req(b)
      omicsCore::plot_enrichment(b, view = "dot", top_n = 12L,
                                 p_preference = show_p(),
                                 p_cutoff = show_cutoff())
    }))

    output$hits <- DT::renderDT({
      df <- table_data()
      shiny::req(nrow(df) > 0L)
      out <- enrich_hits_table(df, show_p(), show_cutoff())
      shiny::req(nrow(out) > 0L)
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

    output$download_table <- shiny::downloadHandler(
      filename = function() {
        b <- if (isTRUE(is_demo())) NULL else enrich_bundle()
        type <- toupper(b$params$type %||% "enrichment")
        sprintf("%s_%s.csv", tolower(type), format(Sys.Date(), "%Y%m%d"))
      },
      content = function(file) {
        df <- table_data()
        shiny::req(nrow(df) > 0L)

        # Unfiltered and every column, including the gene lists -- those
        # are what a follow-up actually needs, and they are the one thing
        # the on-screen table cannot show. Sorted so the file opens on
        # the strongest result rather than in database order.
        df <- df[order(df$database,
                       df$adj_p_value %||% df$p_value,
                       df$p_value), , drop = FALSE]
        utils::write.csv(df, file, row.names = FALSE, na = "")
      }
    )

    # Expose the bundle for slice 3F (report).
    list(
      bundle = shiny::reactive(enrich_bundle()),
      compare = shiny::reactive(compare_bundle())
    )
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

# The rows of the Enriched sets table: the pathways that pass the
# display threshold, strongest first.
enrich_hits_table <- function(df, show_p, show_cutoff) {
  df <- omicsCore::filter_enrich_results(
    df, p_cutoff = show_cutoff, p_preference = show_p)
  pcol <- if (identical(show_p, "raw")) "p_value" else "adj_p_value"
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
    if (identical(show_p, "raw")) "p" else "adjusted p"
  if (!nrow(df)) return(out)
  # ORA's direction is the gene list a pathway was found in, which is
  # worth saying in words: "up" next to an ORA pathway read as the
  # pathway's activity going up, which ORA does not measure.
  if (all(is.na(df$effect))) {
    out$Direction <- ifelse(df$direction %in% "up", "up-regulated genes",
                     ifelse(df$direction %in% "down", "down-regulated genes",
                            NA_character_))
    names(out)[names(out) == "Direction"] <- "Found among"
  }
  # ORA has no enrichment score, and a direction only when it was run on
  # one direction or on both separately; columns of NA said nothing.
  if (all(is.na(df$effect))) out$NES <- NULL
  if (all(is.na(df$direction))) {
    out <- out[setdiff(names(out), c("Direction", "Found among"))]
  }
  out
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
          help = paste("ORA tests whether the differential hits are over-represented",
                       "in a pathway. GSEA ranks every feature by effect and needs",
                       "no hit threshold."),
          shiny::radioButtons(
            ns("type"), label = NULL,
            choices  = c("ORA" = "ora", "GSEA" = "gsea"),
            selected = "ora", inline = TRUE
          )
        ),
        param_group(
          "Species",
          help = "Which species' gene sets to use. Suggested from how the gene names are written.",
          shiny::uiOutput(ns("ui_organism"))
        ),
        param_group(
          "Database",
          shiny::selectInput(
            ns("database"),
            label = htmltools::tags$span(class = "visually-hidden-label", "Database"),
            choices  = c("MSigDB Hallmark" = "hallmark",
                         "KEGG" = "kegg",
                         "Reactome" = "reactome",
                         "GO Biological Process" = "go_bp",
                         "GO Molecular Function" = "go_mf",
                         "GO Cellular Component" = "go_cc",
                         "WikiPathways" = "wikipathways"),
            selected = "hallmark"
          )
        ),
        param_group(
          "Direction",
          help = paste("Which hits ORA tests. Separately (the default) tests the",
                       "up-regulated and the down-regulated hits as two lists, so",
                       "opposite changes are not mixed in one pathway; pooled tests",
                       "them as one list. GSEA always keeps up and down apart."),
          shiny::radioButtons(
            ns("direction"), label = NULL,
            choices  = ENRICH_DIRECTION_CHOICES,
            selected = "separate", inline = TRUE
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
            ns("show_cutoff"),
            label = htmltools::tags$span(class = "visually-hidden-label", "Display p cutoff"),
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
                           "top 12 by adjusted p \u00B7 size = overlap")
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
                           "ranked by adjusted p"),
      # Everything, not what the table happens to be showing: the CSV is
      # what gets opened in Excel a week later, and a file silently
      # truncated to one significance threshold is the kind of thing
      # nobody notices until a number cannot be reproduced.
      shiny::downloadButton(
        ns("download_table"), "Download full table",
        class = "btn btn-ghost btn-sm"
      )
    ),
    bslib::card_body(
      DT::DTOutput(ns("hits"))
    )
  )
}

# Supplied to the worker function's environment by detached_call(), which
# codetools cannot see.
utils::globalVariables("bundle")

# The Direction control, in plain words. "both" keeps its old meaning --
# up and down pooled -- so a project saved before "separate" existed
# restores to what it ran.
ENRICH_DIRECTION_CHOICES <- c("Up and down separately" = "separate",
                              "Up only" = "up",
                              "Down only" = "down",
                              "Up and down pooled" = "both")

# A differential result that says whether features change, not which way
# (a global test, a spline fit): it can only be enriched as one list.
diff_undirected <- function(bundle) {
  at <- bundle$results$diff_result_df$analysis_type
  any(at %in% c("anova", "continuous_spline"))
}

# The species offered, as label = code, from omicsCore's own list so the
# menu and the engine cannot drift apart.
species_choices <- function() {
  sp <- omicsCore::enrichment_species()
  stats::setNames(sp$code, sp$label)
}

species_code <- function(organism) {
  if (is.null(organism) || !length(organism) || is.na(organism[[1L]])) return(NULL)
  sp <- omicsCore::enrichment_species()
  hit <- match(organism[[1L]], c(sp$species, sp$code))
  if (is.na(hit)) return(NULL)
  c(sp$code, sp$code)[[hit]]
}

species_label <- function(code) {
  sp <- omicsCore::enrichment_species()
  lab <- sp$label[match(code, sp$code)]
  if (is.na(lab)) code else sub(" \\(.*\\)$", "", lab)
}

species_example <- function(code) {
  switch(code, Hs = "TP53", Mm = "Trp53", Rn = "Tp53", Dr = "tp53",
         Dm = "p53", Sc = "RAD9", Ce = "cep-1", "TP53")
}

# Which species a set of gene symbols most likely comes from, by how the
# names are written: human upper case (TP53), mouse and rat capitalised
# (Trp53 -- the two cannot be told apart this way, so mouse is
# suggested), zebrafish lower case (tp53), worm lower case with a dash and
# number (unc-54), fly with many CG-numbered genes (CG1824), yeast with
# systematic ORF names (YAL001C).
guess_organism <- function(symbols) {
  x <- symbols[is_gene_symbol(symbols)]
  if (length(x) < 10L) return("Hs")
  orf_like <- mean(grepl("^Y[A-P][LR][0-9]{3}[WC](-[A-Z])?$", x))
  if (orf_like > 0.3) return("Sc")
  cg_like <- mean(grepl("^C[GR][0-9]+$", x))
  if (cg_like > 0.2) return("Dm")
  worm_like <- mean(grepl("^[a-z]{2,5}-[0-9]+(\\.[0-9]+)?$", x))
  if (worm_like > 0.5) return("Ce")
  lower_like <- mean(x == tolower(x) & grepl("[a-z]", x))
  mouse_like <- mean(grepl("^[A-Z][a-z0-9]+[a-z0-9-]*$", x) & grepl("[a-z]", x))
  human_like <- mean(x == toupper(x))
  if (lower_like > 0.5 && lower_like > human_like) return("Dr")
  if (mouse_like > 0.5 && mouse_like > human_like) "Mm" else "Hs"
}

# Values that look like gene symbols rather than accessions or blanks.
is_gene_symbol <- function(x) {
  x <- as.character(x)
  !is.na(x) & nzchar(x) & grepl("^[A-Za-z][A-Za-z0-9.-]{0,19}$", x) &
    !grepl("^(ENS[A-Z]*[0-9]{6,}|[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]|N[MRP]_|X[MR]_)", x)
}
