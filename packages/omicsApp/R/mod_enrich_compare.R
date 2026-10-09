# The "Across comparisons" card of the enrichment view.
#
# Split out of mod_enrich_view.R. The server half is a plain function
# called from inside enrich_view_server()'s moduleServer(), not a module,
# so the ids stay "enrich-run_compare", "enrich-compare_plot". It returns
# the comparison result and its state, which the module resets when the
# layer is replaced and refills when a project is restored.
enrich_compare_server <- function(input, output, session, diff_all, diff_thresholds,
                                  have_cp, organism, run_direction) {
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
          shiny::plotOutput(session$ns("compare_plot"), height = "auto")
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

  # Built by its own reactive so the figure can be reused (a download,
  # say) without drawing it twice.
  compare_plot <- shiny::reactive({
    b <- compare_bundle()
    shiny::req(b)
    omicsCore::plot_enrichment_comparison(
      b, p_preference = input$show_p %||% "adjusted")
  })

  # As tall as its rows need. At a fixed 480 px, five comparisons'
  # pathways (17 rows) left each row 12 px: the names overlapped and the
  # dots of neighbouring rows ran into each other. On a phone the legends
  # sit under the panel and take their own height.
  compare_plot_height <- function() {
    p <- tryCatch(compare_plot(), error = function(e) NULL)
    n <- if (inherits(p, "ggplot") && ".row" %in% names(p$data)) {
      length(unique(p$data$.row))
    } else 0L
    w <- session$clientData[[paste0("output_", session$ns("compare_plot"), "_width")]]
    narrow <- is.numeric(w) && length(w) && w < NARROW_PLOT_PX
    compare_plot_px(n, narrow)
  }

  output$compare_plot <- shiny::renderPlot(
    res = PLOT_RES, height = compare_plot_height,
    alt = "Pathways enriched in each comparison, side by side",
    fit_to_width("compare_plot", compare_plot()))

  list(compare_bundle = compare_bundle, compare_error = compare_error,
       all_comparisons = all_comparisons, compare_epoch = compare_epoch,
       compare_running = compare_running, compare_plot = compare_plot)
}

# Height of the comparison plot for `n_rows` pathways: room for the
# title, subtitle and slanted comparison names, then about 24 px a row
# (16 on a phone, where the text is smaller) plus, on a phone, the
# legends under the panel. Never below the 480 px the card had.
compare_plot_px <- function(n_rows, narrow = FALSE) {
  px <- if (narrow) 330 + 16 * n_rows else 190 + 24 * n_rows
  max(480, round(px))
}
