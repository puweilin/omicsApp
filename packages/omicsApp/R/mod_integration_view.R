#' Multi-omics integration view module
#'
#' Pairs the layer the Differential view ran on (the *primary*) with a
#' second layer of the project (the *partner*) and integrates them by one
#' of three methods:
#'
#' * **Fold-change concordance** -- re-runs the primary's contrast on the
#'   partner (same group column, control, case, covariates and, for a
#'   multi-contrast run, the same set of groups in the model) and compares
#'   the two differential results feature by feature.
#' * **Sample-level correlation** -- correlates each shared gene across
#'   the samples the pairing card says are the same person. Needs no
#'   differential result.
#' * **ActivePathways** -- pathway-level evidence merged across the two
#'   differential results (only offered when the package is installed).
#'
#' With no project the view shows the built-in demo fixture, labelled as
#' such. With a project it never does: a failed or impossible run shows
#' why, not somebody else's numbers.
#'
#' @param id Module namespace id.
#' @keywords internal
#' @noRd
integration_view_ui <- function(id) {
  ns <- shiny::NS(id)

  htmltools::tagList(
    shiny::uiOutput(ns("header")),
    shiny::uiOutput(ns("notices")),
    integration_setup_card(ns),
    # Before the results, not after: which sample was treated as which
    # person is an assumption the reader should see before reading
    # anything computed on it.
    integration_pairing_card(ns),
    # Which protein is which gene, for the same reason: every method
    # compares matched features, and isoforms of one gene each count.
    integration_feature_card(ns),
    shiny::uiOutput(ns("stats")),
    shiny::uiOutput(ns("results"))
  )
}

#' @rdname integration_view_ui
#' @param current_project Reactive yielding the live `omics_project`.
#' @param diff_bundle Reactive yielding the primary differential bundle
#'   (one contrast).
#' @param diff_layer Reactive yielding the project tag of the layer the
#'   differential bundle was computed on.
#' @param diff_thresholds Reactive yielding the Differential view's
#'   `list(p_cutoff, p_preference, effect_cutoff)`, so a "hit" means the
#'   same thing in both views.
#' @keywords internal
#' @noRd
integration_view_server <- function(id,
                                    current_project = shiny::reactiveVal(NULL),
                                    diff_bundle = shiny::reactiveVal(NULL),
                                    invalidate = shiny::reactiveVal(0L),
                                    diff_layer = shiny::reactiveVal(NULL),
                                    diff_thresholds = shiny::reactiveVal(NULL),
                                    navigate = NULL) {
  shiny::moduleServer(id, function(input, output, session) {

    method <- shiny::reactive(input$method %||% "concordance")

    # ---- which two layers --------------------------------------------
    # The primary is the layer the diff ran on -- the Differential view
    # says which. It used to be guessed from the bundle's omics type,
    # which picked the wrong layer whenever two layers shared a type.
    layers <- shiny::reactive({
      proj <- current_project()
      if (is.null(proj) || length(proj$experiments) < 2L) return(NULL)
      tags <- names(proj$experiments)
      primary <- diff_layer()
      if (is.null(primary) || !primary %in% tags) {
        hit <- omicsCore::bundle_layer(proj, diff_bundle())
        primary <- if (is.na(hit)) tags[[1L]] else hit
      }
      others <- setdiff(tags, primary)
      want <- input$partner
      partner <- if (!is.null(want) && want %in% others) want else others[[1L]]
      list(primary = primary, partner = partner, others = others)
    })

    output$ui_partner <- shiny::renderUI({
      l <- layers()
      if (is.null(l)) return(NULL)
      if (length(l$others) < 2L) {
        return(htmltools::tags$div(
          class = "muted", style = "font-size:12.5px;padding-top:6px",
          sprintf("%s \u00D7 %s", l$primary, l$partner)))
      }
      shiny::selectInput(session$ns("partner"),
                         label = sprintf("Integrate %s with", l$primary),
                         choices = l$others,
                         selected = shiny::isolate(l$partner))
    })

    # ---- prerequisites ------------------------------------------------
    # Everything a run needs, or the one reason it cannot happen
    # (mod_integration_run.R).
    can_run <- integration_can_run(current_project, layers, method, diff_bundle)

    integration_bundle <- shiny::reactiveVal(NULL)
    integration_error  <- shiny::reactiveVal(NULL)
    is_demo            <- shiny::reactive(is.null(current_project()))
    running            <- shiny::reactiveVal(FALSE)
    # The partner's own diff is the expensive half of a concordance run
    # and depends only on the partner and the design, so it is kept and
    # reused when only the thresholds (or the method) change.
    sec_cache <- shiny::reactiveVal(NULL)

    # ---- sample pairing and feature matching ----------------------------
    # The two cards between the setup and the results
    # (mod_integration_pairing.R).
    pairing <- integration_pairing_server(input, output, session, current_project, layers)
    features <- integration_feature_server(input, output, session, current_project, layers)
    feature_pairing <- features$feature_pairing
    link_upload     <- features$link_upload
    pending_link    <- features$pending_link

    # One of the layers this integration spanned has been replaced, so
    # the pairing it reports no longer exists. Back to the module's own
    # start-up state.
    shiny::observeEvent(invalidate(), {
      integration_bundle(NULL)
      integration_error(NULL)
      sec_cache(NULL)
      runs$last_key <- NULL
    }, ignoreInit = TRUE)

    # ---- running ------------------------------------------------------
    # Plain state, not reactive: which inputs the last run saw, and a
    # counter that lets a newer run win over an older one still in flight.
    runs <- new.env(parent = emptyenv())
    runs$last_key <- NULL
    runs$token <- 0L

    set_busy <- function(busy) {
      running(busy)
      tryCatch(
        if (busy) shinyjs::disable("rerun") else shinyjs::enable("rerun"),
        error = function(e) NULL)
    }

    # What a run depends on, and the run itself (mod_integration_run.R).
    run <- integration_run_server(input, output, session, current_project, diff_bundle,
                                  diff_thresholds, method, can_run, runs, sec_cache,
                                  set_busy, integration_bundle, integration_error)
    run_key <- run$run_key
    do_run  <- run$do_run

    shiny::observeEvent(input$rerun, {
      sec_cache(NULL)
      do_run()
    })
    shiny::observeEvent(run_key(), {
      key <- run_key()
      if (identical(key, runs$last_key)) return()
      do_run()
    })
    # Prerequisites lost (layer removed, diff cleared): the old result is
    # about something that is no longer on screen.
    shiny::observeEvent(can_run(), {
      if (!isTRUE(can_run()$ok)) {
        integration_bundle(NULL)
        runs$last_key <- NULL
      }
    }, ignoreInit = TRUE)

    # ---- what is shown ------------------------------------------------
    # The header, the notices, the stat cards, the plots and the tables
    # (mod_integration_results.R).
    figures <- integration_results_server(input, output, session, navigate, method, can_run,
                                          is_demo, integration_bundle, integration_error,
                                          running)

    # Each figure as a file (plot_download.R), named for the project, the
    # figure and the two layers. None for the demo, which is no one's data.
    integration_file <- function(what) function() {
      c(plot_download_project(current_project()), what,
        paste(integration_bundle()$params$experiments, collapse = "_"))
    }
    not_demo <- function(plot) shiny::reactive(!isTRUE(is_demo()) && plot_ready(plot))
    plot_download_server("scatter_download", figures$scatter,
                         integration_file("fold_change_concordance"),
                         width_in = 7, height_in = 5.5, available = not_demo(figures$scatter))
    plot_download_server("top_hits_download", figures$top_hits,
                         integration_file("top_hits_both_layers"),
                         width_in = 7, height_in = 5, available = not_demo(figures$top_hits))
    plot_download_server("cor_scatter_download", figures$cor_scatter,
                         integration_file("correlation_per_gene"),
                         width_in = 7, height_in = 5, available = not_demo(figures$cor_scatter))
    plot_download_server("ap_dot_download", figures$ap_dot,
                         integration_file("activepathways"),
                         width_in = 8, height_in = 5.5, available = not_demo(figures$ap_dot))

    # Exposed for the report and the project.
    shiny::reactive(integration_bundle())
  })
}

# ---- internal helpers ------------------------------------------------

utils::globalVariables(".data")

integration_setup_card <- function(ns) {
  choices <- c("Fold-change concordance" = "concordance",
               "Sample-level correlation" = "correlation")
  if (has_pkg("ActivePathways")) {
    choices <- c(choices, "ActivePathways (pathways)" = "active_pathways")
  }
  bslib::card(
    bslib::card_header(
      htmltools::tags$h3(class = "card-title", "Setup"),
      htmltools::tags$span(class = "card-sub", "which layers, and how")
    ),
    bslib::card_body(
      htmltools::tags$div(
        class = "row-grid r-6-6",
        htmltools::tags$div(
          htmltools::tags$label(class = "control-label", "Layers"),
          shiny::uiOutput(ns("ui_partner"))
        ),
        shiny::radioButtons(
          ns("method"), label = "Method", choices = choices,
          selected = "concordance", inline = TRUE
        )
      ),
      htmltools::tags$div(
        class = "muted", style = "font-size:12px",
        paste("Concordance repeats the Differential view's contrast on the",
              "second layer and compares the two results gene by gene, at the",
              "same thresholds. Correlation compares the two layers sample by",
              "sample and needs the pairing below.",
              if (has_pkg("ActivePathways"))
                paste("ActivePathways merges both layers' evidence pathway by",
                      "pathway, expecting the layers to change the same way:",
                      "genes that go up in one and down in the other count",
                      "against a pathway, and each pathway is reported as up,",
                      "down, or mixed."))
      )
    )
  )
}
