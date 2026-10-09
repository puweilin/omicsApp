#' QC view module
#'
#' Slice 3C: the view now reacts to the shared `current_project`
#' reactiveVal *and* a small controls panel (missing-rate slider
#' + outlier-method radio). When the project is `NULL` it falls
#' back to `example_qc_bundle()`; when it has at least one
#' experiment, we pick the active layer and re-run
#' `omicsCore::run_qc()` live against the user's inputs.
#'
#' Per slice-3 convention QC runs on every input change (no Run
#' button) -- it's cheap and the feedback loop is more useful that
#' way. A project that arrives with a saved QC result shows that result
#' instead, with its settings in the controls, until one is changed.
#'
#' Reference markup: `omicsApp/mockup/index.html:686-750`.
#'
#' @param id Module namespace id.
#' @keywords internal
#' @noRd
qc_view_ui <- function(id) {
  ns <- shiny::NS(id)
  htmltools::tagList(
    shiny::uiOutput(ns("header")),
    qc_controls_card(ns),
    shiny::uiOutput(ns("notices")),
    shiny::uiOutput(ns("stats")),
    htmltools::tags$div(
      class = "row-grid r-7-5",
      bslib::card(
        bslib::card_header(
          htmltools::tags$h3(class = "card-title", "PCA"),
          htmltools::tags$span(
            class = "card-sub",
            "samples projected on PC1 \u00D7 PC2"
          )
        ),
        bslib::card_body(
          shiny::uiOutput(ns("pca_color_picker")),
          hover_plot_output(ns("pca"), height = "360px")
        )
      ),
      bslib::card(
        bslib::card_header(
          shiny::uiOutput(ns("quality_title"), inline = TRUE),
          shiny::uiOutput(ns("quality_picker"), inline = TRUE)
        ),
        bslib::card_body(
          shiny::plotOutput(ns("missing"), height = "360px"),
          shiny::uiOutput(ns("missing_caption"))
        )
      )
    )
  )
}

#' @rdname qc_view_ui
#' @param current_project Reactive (or reactiveVal) yielding the
#'   live `omics_project` or `NULL`.
#' @keywords internal
#' @noRd
qc_view_server <- function(id, current_project = shiny::reactiveVal(NULL),
                           invalidate = shiny::reactiveVal(0L),
                           requested_layer = shiny::reactiveVal(NULL),
                           navigate = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    if (is.function(navigate)) {
      shiny::observeEvent(input$go_next, navigate("diff"))
    }

    # ---- a saved result -------------------------------------------------
    # A project opened or restored brings its QC result, computed at the
    # settings the user chose then. Recomputing it on arrival, at whatever
    # the controls happened to say, replaced that result -- in the project
    # and in the next autosave -- with one nobody asked for. So the saved
    # result is shown, its settings are put back into the controls, and QC
    # runs again only once the user changes one of them or moves to a
    # layer it was not computed on.
    #
    # restore() holds the saved bundle, the layer it belongs to, its
    # settings as a request (see qc_request below), and, per control, the
    # value saved and the value the control held when the project arrived.
    # Until a control is seen to take its saved value (or the user moves
    # it somewhere else), its old value is read as the saved one: the
    # browser applies the update a moment later, and in that moment the
    # old value is not a request to recompute.
    restore <- shiny::reactiveVal(NULL)
    pick <- function(name, value) {
      r <- restore()
      if (is.null(r) || !name %in% names(r$saved)) return(value)
      if (identical(value, r$old[[name]])) r$saved[[name]] else value
    }

    # Active experiment selection. A tag asked for from elsewhere (the
    # Project view's "View" link) wins, provided it still names a layer
    # in the current project -- otherwise the request is stale and
    # silently ignored rather than emptying the view. With no request,
    # pick the first proteomics layer (PCA + missingness panels are
    # designed for that); fall back to the first experiment of any
    # kind; fall back to the built-in proteomics fixture.
    active <- shiny::reactive({
      proj <- current_project()
      if (is.null(proj)) {
        return(list(input = NULL, tag = NULL, is_demo = TRUE))
      }
      exps <- proj$experiments
      if (length(exps) == 0L) {
        return(list(input = NULL, tag = NULL, is_demo = TRUE))
      }
      # The picker wins when it names a layer that exists: it is the one
      # the user is looking at. requested_layer() is how another view
      # hands over ("show me this one"), and it seeds the picker rather
      # than fighting it.
      want <- pick("layer", input$layer)
      if (is.null(want) || !want %in% names(exps)) want <- requested_layer()
      if (!is.null(want) && length(want) == 1L && want %in% names(exps)) {
        return(list(input = exps[[want]], tag = want, is_demo = FALSE))
      }
      idx <- qc_default_layer_idx(exps)
      list(
        input   = exps[[idx]],
        tag     = names(exps)[idx],
        is_demo = FALSE
      )
    })

    output$ui_layer <- shiny::renderUI({
      proj <- current_project()
      if (is.null(proj) || length(proj$experiments) < 2L) return(NULL)
      tags_avail <- names(proj$experiments)
      # isolate(), for the reason the Differential view documents: this
      # output writes input$layer and active() reads it, so reading it
      # here would close the loop -- re-rendering the control re-sends
      # its value, which invalidates active(), which re-renders it.
      sel <- shiny::isolate(pick("layer", input$layer))
      if (is.null(sel) || !sel %in% tags_avail) {
        sel <- shiny::isolate(requested_layer())
      }
      if (is.null(sel) || !sel %in% tags_avail) {
        sel <- tags_avail[[qc_default_layer_idx(proj$experiments)]]
      }
      shiny::selectInput(session$ns("layer"), label = "Omics layer",
                         choices = tags_avail, selected = sel)
    })

    # The QC bundle. tryCatch keeps a parameter mishap (e.g. a
    # slider that drops all features) from blowing up the view --
    # the notices panel surfaces the message and the previous
    # plots stay visible.
    last_bundle <- shiny::reactiveVal(NULL)
    last_error  <- shiny::reactiveVal(NULL)
    # What the result on screen was computed from (see qc_request()). A
    # change to the project that leaves it alone -- another view's result
    # being attached -- is not a reason to run QC again.
    shown_request <- NULL

    # The layer this result was computed on has been replaced, so the
    # result is no longer about anything in the project. NULL is the
    # module's own start-up state, so this only rewinds it.
    shiny::observeEvent(invalidate(), {
      last_bundle(NULL)
      last_error(NULL)
      restore(NULL)
      shown_request <<- NULL
    }, ignoreInit = TRUE)

    # The Filters card: imputation, the missing-value filter, and the
    # request they add up to (mod_qc_controls.R).
    controls <- qc_controls_server(input, output, session, active, pick)
    impute_choices        <- controls$impute_choices
    impute_applies        <- controls$impute_applies
    thr_r                 <- controls$thr_r
    missing_group_choices <- controls$missing_group_choices
    qc_request            <- controls$qc_request

    # Priority -6: after the restore below (-5), which itself runs after
    # a generation bump has cleared the view, so a project arriving with
    # a saved result is matched against it before anything is computed.
    shiny::observe({
      req <- qc_request()
      r <- restore()
      if (!is.null(r)) {
        # Controls that have taken their saved value, or been moved
        # elsewhere, no longer stand in for it.
        now <- list(layer = input$layer, missing_threshold = thr_r(),
                    outlier_method = input$outlier_method,
                    impute_method = input$impute_method,
                    missing_filter = input$missing_filter,
                    missing_group_col = input$missing_group_col)
        settled <- vapply(names(r$saved), function(nm) {
          identical(now[[nm]], r$saved[[nm]]) || !identical(now[[nm]], r$old[[nm]])
        }, logical(1))
        if (any(settled)) {
          r$saved <- r$saved[!settled]
          r$old <- r$old[!settled]
          restore(r)
        }
        # Asked for exactly what was saved, on the layer it was saved
        # from: that result, not a new run of it.
        if (identical(req$tag, r$tag) && identical(req$input, r$input) &&
            identical(req$params, r$params)) {
          shown_request <<- req
          last_error(NULL)
          last_bundle(r$bundle)
          return(invisible())
        }
      }
      if (identical(req, shown_request)) return(invisible())
      shown_request <<- req
      qc_input <- req$input
      p <- req$params

      run <- function() omicsCore::run_qc(
        qc_input,
        missing_threshold = p$missing_threshold,
        outlier_method    = p$outlier_method,
        impute_method     = p$impute_method,
        missing_filter    = p$missing_filter,
        group_col         = p$group_col
      )
      # A big layer takes several seconds the first time (up to 10 s for
      # 60 samples), and all it showed was the busy pulse. It now says
      # which step it is on. A small one finishes before a panel could
      # be read, so it gets none rather than a flash.
      big <- length(qc_input$expr_mat) >= QC_PROGRESS_CELLS
      bundle <- tryCatch(
        if (big) with_step_progress("Running quality control", run()) else run(),
        error = function(e) e)

      if (inherits(bundle, "error")) {
        last_error(conditionMessage(bundle))
      } else {
        # The layer's name travels with the result (a project can hold
        # two layers of one omics type); the demo's layers are not the
        # user's.
        if (!is.null(current_project())) bundle$input_info$layer <- req$tag
        last_error(NULL)
        last_bundle(bundle)
      }
    }, priority = -6)

    # A project opened or restored brings its QC result: shown as it was
    # saved, on the layer it was computed on, with its settings back in
    # the controls. After the generation bump and the clearing it causes
    # (priority -5, as in the Differential view), so that does not undo
    # this. A result this view published itself comes back through the
    # project too, and is recognised as the one already on screen.
    shiny::observeEvent(current_project(), {
      proj <- current_project()
      saved <- proj$bundles$qc
      if (!omicsCore::is_analysis_bundle(saved) || identical(saved, last_bundle())) return()
      tag <- qc_bundle_layer(proj$experiments, saved)
      if (is.null(tag)) return()
      p <- saved$params
      saved_vals <- list(
        layer = tag,
        missing_threshold = p$missing_threshold,
        outlier_method = qc_outlier_choice(p$outlier_method),
        impute_method = p$impute_method,
        missing_filter = p$missing_filter %||% "global",
        missing_group_col = p$group_col
      )
      saved_vals <- saved_vals[!vapply(saved_vals, is.null, logical(1))]
      old <- list(layer = input$layer, missing_threshold = thr_r(),
                  outlier_method = input$outlier_method,
                  impute_method = input$impute_method,
                  missing_filter = input$missing_filter,
                  missing_group_col = input$missing_group_col)[names(saved_vals)]
      restore(list(
        bundle = saved, tag = tag, input = proj$experiments[[tag]],
        params = qc_request_params(p$missing_threshold, p$outlier_method,
                                   p$impute_method, p$missing_filter %||% "global",
                                   p$group_col),
        saved = saved_vals, old = old))
      shown_request <<- NULL
      last_error(NULL)
      last_bundle(saved)
      qc_restore_controls(session, saved_vals)
    }, priority = -5, ignoreNULL = TRUE)

    # The header, the notices (with outlier exclusion) and the stat cards
    # (mod_qc_results.R).
    qc_results_server(input, output, session, navigate, current_project, active,
                      last_bundle, last_error)

    # The PCA and the missingness / depth panel (mod_qc_plots.R).
    plots <- qc_plots_server(input, output, session, active, last_bundle)
    pca_color_choices <- plots$pca_color_choices
    quality_view      <- plots$quality_view
    pca_plot          <- plots$pca_plot

    # Expose the QC bundle for slice 3F (report).
    shiny::reactive(last_bundle())
  })
}

# ---- internal helpers ------------------------------------------------

# The radio value that runs these outlier methods: "all", a single
# method, or the methods themselves when the radio has no button for
# them (a project saved before the leave-one-out test, or a script).
qc_outlier_choice <- function(methods) {
  if (is.null(methods)) return(NULL)
  if (setequal(methods, QC_ALL_OUTLIER_METHODS)) return("all")
  methods
}

# The layer a saved QC bundle was computed on: same kind, same shape, and
# its samples among the layer's. NULL when no layer fits -- the result is
# then about data the project no longer holds.
qc_bundle_layer <- function(exps, bundle) {
  info <- bundle$input_info
  kept <- qc_kept_samples(bundle)
  fits <- function(e) {
    identical(e$omics_type %||% "", info$omics_type %||% "") &&
      isTRUE(ncol(e$expr_mat) == info$n_samples_in) &&
      isTRUE(nrow(e$expr_mat) == info$n_features_in) &&
      all(kept %in% colnames(e$expr_mat))
  }
  # The layer it records, when that layer still holds the same data;
  # otherwise (older results) the first layer that fits.
  rec <- info$layer
  if (length(rec) == 1L && !is.na(rec)) {
    return(if (rec %in% names(exps) && fits(exps[[rec]])) rec else NULL)
  }
  for (tag in names(exps)) {
    e <- exps[[tag]]
    if (identical(e$omics_type %||% "", info$omics_type %||% "") &&
        isTRUE(ncol(e$expr_mat) == info$n_samples_in) &&
        isTRUE(nrow(e$expr_mat) == info$n_features_in) &&
        all(kept %in% colnames(e$expr_mat))) {
      return(tag)
    }
  }
  NULL
}

# What a QC result says about the cleaned input, read from the record
# run_qc() keeps -- or, for results saved before it stopped storing the
# cleaned input, from that copy.
qc_kept_samples <- function(bundle) {
  bundle$results$cleaning$kept_samples %||%
    colnames(bundle$results$cleaned_input$expr_mat)
}

qc_bundle_meta <- function(bundle) {
  bundle$results$plot_data$meta_df %||% bundle$results$cleaned_input$meta_df
}

# Put a saved result's settings back into the controls. The controls a
# result cannot be expressed in (several outlier methods the radio has no
# button for) are left as they are.
qc_restore_controls <- function(session, vals) {
  if (!is.null(vals$layer)) shiny::updateSelectInput(session, "layer", selected = vals$layer)
  if (!is.null(vals$missing_threshold)) {
    shiny::updateSliderInput(session, "missing_threshold", value = vals$missing_threshold)
  }
  if (length(vals$outlier_method) == 1L) {
    shiny::updateRadioButtons(session, "outlier_method", selected = vals$outlier_method)
  }
  if (!is.null(vals$impute_method)) {
    shiny::updateSelectInput(session, "impute_method", selected = vals$impute_method)
  }
  if (!is.null(vals$missing_filter)) {
    shiny::updateSelectInput(session, "missing_filter", selected = vals$missing_filter)
  }
  if (!is.null(vals$missing_group_col)) {
    shiny::updateSelectInput(session, "missing_group_col", selected = vals$missing_group_col)
  }
  invisible(vals)
}

# Proteomics first when nothing else has been asked for: the missingness
# panels are the ones this view was built around, and they say nothing
# about a counts matrix.
qc_default_layer_idx <- function(exps) {
  types <- vapply(exps, function(e) e$omics_type %||% "", character(1L))
  idx <- which(types == "proteomics")
  if (length(idx) == 0L) 1L else idx[[1L]]
}

# Matrix cells above which a QC run shows its steps: about 5,000
# features x 40 samples, which takes a second or more.
QC_PROGRESS_CELLS <- 2e5
