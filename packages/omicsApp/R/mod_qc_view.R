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
          shiny::plotOutput(ns("pca"), height = "360px")
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

    # Imputation is offered for proteomics and withheld everywhere else.
    #
    # A missing intensity in DIA usually means "below the detection
    # limit". Leaving it NA is not the neutral option it looks like:
    # limma drops a feature it cannot fit, so "none" is complete-case
    # analysis chosen silently. run_qc() therefore defaults proteomics to
    # MinProb and this control opens on the same value.
    #
    # A missing count is a different thing. A zero is an observation, and
    # imputing counts feeds DESeq2 numbers its model never saw -- so the
    # control is not offered there at all.
    #
    # DEP's method set and DEP's spelling, so a choice made here means
    # what it means in the proteomics literature and in every paper the
    # analyst has read. Grouped by assumption, because that is the choice
    # actually being made: MNAR says a value is missing *because* it was
    # low, MAR says it is missing for reasons unrelated to its size.
    impute_choices <- shiny::reactive({
      grouped <- list(
        "Left-censored (MNAR)" = c(
          "MinProb \u2014 draw near the minimum" = "MinProb",
          "MinDet \u2014 low quantile"           = "MinDet",
          "QRILC \u2014 quantile regression"     = "QRILC",
          "min \u2014 feature minimum"           = "min",
          "zero"                                 = "zero"),
        "Random (MAR)" = c(
          "knn \u2014 k-nearest neighbours" = "knn",
          "MLE \u2014 maximum likelihood"   = "MLE",
          "bpca \u2014 Bayesian PCA"        = "bpca"),
        "Other" = c(
          "mixed \u2014 MAR/MNAR per feature" = "mixed",
          "man \u2014 manual shift/scale"     = "man",
          "none \u2014 leave NA"              = "none")
      )
      # An option that errors on selection is worse than one not offered:
      # the backends stop with an install hint, which arrives as a red
      # notice over a view that was working a moment ago.
      needs <- c(MinProb = "imputeLCMD", MinDet = "imputeLCMD",
                 QRILC = "imputeLCMD", knn = "imputeLCMD",
                 MLE = "imputeLCMD", mixed = "imputeLCMD",
                 bpca = "pcaMethods")
      out <- lapply(grouped, function(g) {
        keep <- vapply(g, function(v) {
          # Single bracket: `needs[["min"]]` errors on a name that is not
          # there, where `needs["min"]` gives NA and lets the method
          # through as needing nothing.
          pkg <- unname(needs[v])
          is.na(pkg) || has_pkg(pkg)
        }, logical(1L))
        g[keep]
      })
      out[lengths(out) > 0L]
    })

    impute_applies <- shiny::reactive({
      a <- active()
      inp <- if (a$is_demo) example_qc_input() else a$input
      identical(inp$omics_type %||% "", "proteomics")
    })

    output$ui_impute <- shiny::renderUI({
      if (!impute_applies()) return(NULL)
      choices <- impute_choices()
      offered <- unlist(choices, use.names = FALSE)
      sel <- shiny::isolate(pick("impute_method", input$impute_method))
      # Defaults to what run_qc() would resolve on its own, so the
      # control opens showing what is actually running rather than
      # imposing a different answer the moment it renders.
      if (is.null(sel) || !sel %in% offered) {
        sel <- omicsCore::resolve_impute_method("proteomics")
      }
      if (!sel %in% offered) sel <- "none"
      shiny::selectInput(session$ns("impute_method"),
                         label = "Imputation (proteomics)",
                         choices = choices, selected = sel)
    })

    # The slider, debounced: QC runs on the main process, and each tick of
    # a drag used to queue a full run_qc() -- 62 s of frozen session for
    # a five-tick drag on 8,000 x 300. (Tests set the delay to 0.)
    qc_delay <- getOption("omicsApp.qc_debounce_ms", 400)
    thr_in <- shiny::reactive(input$missing_threshold)
    thr_r <- if (isTRUE(qc_delay > 0)) shiny::debounce(thr_in, qc_delay) else thr_in

    # The columns the group-wise missing filter can group by: the layer's
    # recorded design first, then the Differential view's guesses. Empty
    # when nothing splits the samples, and the control is then not shown.
    missing_group_choices <- shiny::reactive({
      a <- active()
      inp <- if (a$is_demo) example_qc_input() else a$input
      meta <- inp$meta_df
      if (is.null(meta) || !ncol(meta)) return(character(0))
      cands <- grouping_candidates(meta)
      if (!length(cands)) cands <- grouping_candidates(meta, min_per_level = 1L,
                                                      replicated = TRUE)
      design <- tryCatch(omicsCore::study_design(inp), error = function(e) NULL)
      if (!is.null(design)) cands <- c(design$group_col, setdiff(cands, design$group_col))
      cands
    })

    output$ui_missing_filter <- shiny::renderUI({
      cols <- missing_group_choices()
      if (!length(cols)) return(NULL)
      mode <- shiny::isolate(pick("missing_filter", input$missing_filter))
      if (is.null(mode) || !mode %in% QC_MISSING_FILTERS) mode <- "global"
      gc <- shiny::isolate(pick("missing_group_col", input$missing_group_col))
      if (is.null(gc) || !gc %in% cols) gc <- cols[[1L]]
      htmltools::tags$div(
        class = "row-grid r-6-6",
        shiny::selectInput(
          session$ns("missing_filter"),
          label = htmltools::tagList(
            "Apply the cutoff",
            info_tip(paste(
              "Across all samples, or within each group. \"In at least one group\"",
              "keeps a protein seen in enough samples of one condition, even if it is",
              "absent from the other -- often the most interesting kind."))),
          choices = c("Across all samples" = "global",
                      "In at least one group" = "any_group",
                      "In every group" = "all_groups"),
          selected = mode),
        shiny::conditionalPanel(
          condition = "input.missing_filter != 'global'", ns = session$ns,
          shiny::selectInput(session$ns("missing_group_col"), label = "Groups from",
                             choices = cols, selected = gc))
      )
    })

    # Everything a QC result depends on, as the controls currently say it
    # -- with the values a saved result is waiting for read through
    # pick(), and every default resolved, so that two requests for the
    # same thing compare equal however they were arrived at.
    qc_request <- shiny::reactive({
      a <- active()
      # The demo runs through run_qc() like a real project rather than
      # returning a fixed bundle. Both controls above are enabled, and
      # a control that is enabled and does nothing reads as a broken
      # app; the demo input is 50 x 12, so a re-run is milliseconds.
      qc_input <- if (a$is_demo) example_qc_input() else a$input
      thr <- pick("missing_threshold", thr_r()) %||% 0.5
      out_m <- pick("outlier_method", input$outlier_method) %||% "all"
      # All of them by default: after vsn the per-sample means are equal,
      # and the IQR test on them -- the old default -- could not see a
      # sample that PCA and connectivity both flagged; and with ten or
      # fewer samples only the leave-one-out test can flag anything.
      if (identical(out_m, "all")) out_m <- QC_ALL_OUTLIER_METHODS
      # Read here rather than trusted from the input: the control is
      # hidden when the layer is not proteomics, but Shiny keeps an
      # input's last value, so switching from a proteomics layer with
      # `knn` selected to a counts layer would otherwise impute counts
      # with a control the user can no longer see.
      # Unset resolves the way run_qc() resolves it per modality, so the
      # control and a plain run_qc(input) agree instead of quietly
      # differing.
      imp <- if (!impute_applies()) "none" else
        pick("impute_method", input$impute_method) %||%
          omicsCore::resolve_impute_method(qc_input$omics_type)
      # The same reasoning as imputation: hidden when nothing splits the
      # samples, so a value left over from another layer is not used.
      cols <- missing_group_choices()
      mf <- if (length(cols)) pick("missing_filter", input$missing_filter) %||% "global"
            else "global"
      gc <- NULL
      if (!identical(mf, "global")) {
        gc <- pick("missing_group_col", input$missing_group_col)
        if (is.null(gc) || !gc %in% cols) gc <- cols[[1L]]
      }
      list(tag = a$tag, input = qc_input,
           params = qc_request_params(thr, out_m, imp, mf, gc))
    })

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

    output$header <- shiny::renderUI({
      a <- active()
      bundle <- last_bundle()
      n_in  <- if (is.null(bundle)) NA_integer_ else bundle$input_info$n_samples_in
      n_feat <- if (is.null(bundle)) NA_integer_ else bundle$input_info$n_features_in
      omics <- if (is.null(bundle)) "\u2014" else omics_display(bundle$input_info$omics_type)
      view_header(
        title    = "Quality control",
        actions  = if (is.function(navigate) && !a$is_demo) {
          shiny::actionButton(session$ns("go_next"), "Next: Differential \u2192",
                              class = "btn btn-ghost")
        },
        subtitle = htmltools::tagList(
          omics,
          htmltools::HTML(" &middot; "),
          sprintf("%s features \u00B7 %s samples",
                  format(n_feat, big.mark = ","),
                  format(n_in,  big.mark = ",")),
          htmltools::HTML(" &middot; "),
          htmltools::tags$span(
            class = "muted",
            if (a$is_demo) "demo data (built-in)"
            else sprintf("layer = %s", a$tag)
          )
        )
      )
    })

    output$notices <- shiny::renderUI({
      err <- last_error()
      if (!is.null(err)) {
        return(notice(
          title  = "QC could not run",
          detail = err,
          kind   = "warn"
        ))
      }
      # What run_qc() did to the data on the way (a log scale for the
      # outlier tests, imputation on log2, samples flagged and kept).
      b <- last_bundle()
      notes <- b$warnings
      flagged <- b$results$qc_summary$outliers$flagged_samples
      # In the app's words: the engine's note names R arguments.
      notes <- notes[!grepl("^Flagged as outlier", notes)]
      notes <- sub("`raw_count` values were put on a log2 scale \\(log2-CPM\\) for outlier detection.",
                   "Counts were converted to log2 counts-per-million for the outlier checks.", notes)
      notes <- qc_plain_outlier_notes(notes)
      out <- htmltools::tagList()
      if (length(flagged) && !active()$is_demo) {
        out <- htmltools::tagAppendChild(out, notice(
          title = sprintf("Possible outlier%s: %s", if (length(flagged) > 1L) "s" else "",
                          paste(flagged, collapse = ", ")),
          detail = htmltools::tagList(
            qc_loo_explanation(b$results$qc_summary$outliers),
            "Kept in the analysis. Look at the PCA: if the sample is broken (a failed run, a swap) rather than biologically different, exclude it. ",
            shiny::actionButton(session$ns("exclude_flagged"),
                                "Exclude from this layer\u2026",
                                class = "btn btn-sm btn-outline-danger")),
          kind = "warn"))
      }
      if (length(notes)) {
        out <- htmltools::tagAppendChild(out, notice(
          title  = "About these QC results",
          detail = htmltools::tags$ul(lapply(notes, htmltools::tags$li)),
          kind   = "info"
        ))
      }
      out
    })

    # Excluding a sample replaces the layer by its subset. Results computed
    # with the sample are cleared, and the script repeats the exclusion.
    shiny::observeEvent(input$exclude_flagged, {
      flagged <- last_bundle()$results$qc_summary$outliers$flagged_samples
      shiny::req(length(flagged))
      shiny::showModal(shiny::modalDialog(
        title = "Exclude these samples?",
        htmltools::tags$p(sprintf("%s will be removed from layer '%s'.",
                                  paste(flagged, collapse = ", "), active()$tag)),
        htmltools::tags$p("Results already computed on this layer are cleared and must be re-run. The original file is unchanged."),
        footer = htmltools::tagList(
          shiny::modalButton("Cancel"),
          shiny::actionButton(session$ns("confirm_exclude"), "Exclude",
                              class = "btn btn-danger")),
        easyClose = TRUE))
    })
    shiny::observeEvent(input$confirm_exclude, {
      shiny::removeModal()
      proj <- current_project()
      tag <- active()$tag
      flagged <- last_bundle()$results$qc_summary$outliers$flagged_samples
      shiny::req(proj, tag, length(flagged))
      inp <- proj$experiments[[tag]]
      keep <- setdiff(colnames(inp$expr_mat), flagged)
      new <- tryCatch(omicsCore::subset_omics(inp, samples = keep), error = function(e) e)
      if (inherits(new, "error")) {
        shiny::showNotification(conditionMessage(new), type = "error")
        return()
      }
      new$excluded_samples <- unique(c(inp$excluded_samples, intersect(flagged, colnames(inp$expr_mat))))
      # A new identity, so every view lets go of results computed with
      # the samples still in.
      new$source_fingerprint <- paste0(inp$source_fingerprint %||% "", ":excluded=",
                                       paste(sort(new$excluded_samples), collapse = ","))
      proj$experiments[[tag]] <- new
      proj$bundles <- drop_layer_bundles(proj$bundles, tag, current_project())
      current_project(proj)
      shiny::showNotification(sprintf("Excluded %s from '%s'.", paste(flagged, collapse = ", "), tag),
                              type = "message")
    })

    output$stats <- shiny::renderUI({
      bundle <- last_bundle()
      if (is.null(bundle)) return(NULL)
      info <- bundle$input_info
      summary <- bundle$results$qc_summary
      n_flagged_samp <- length(summary$recommended_filters$remove_samples)
      n_flagged_feat <- length(summary$recommended_filters$remove_features)
      impute <- bundle$params$impute_method %||% "none"
      htmltools::tags$div(
        class = "stat-grid",
        stat_card(
          label  = "Samples kept",
          value  = sprintf("%d / %d", info$n_samples_out, info$n_samples_in),
          # Flagged samples are kept (run_qc(remove_outliers = FALSE)):
          # dropping one changes the design, so the card asks for a look
          # rather than reporting a removal that did not happen.
          trend  = if (n_flagged_samp == 0L) "no outliers flagged"
                   else sprintf("%d flagged (%s) \u2014 kept; check the PCA",
                                n_flagged_samp,
                                qc_method_labels(summary$outliers$method)),
          accent = if (n_flagged_samp == 0L) "ok" else "warn",
          mono   = TRUE
        ),
        stat_card(
          label = "Features kept",
          value = format(info$n_features_out, big.mark = ","),
          trend = sprintf("%d filtered at %.0f%% missing%s",
                          n_flagged_feat,
                          100 * (bundle$params$missing_threshold %||% 0.5),
                          switch(bundle$params$missing_filter %||% "global",
                                 any_group = " in every group",
                                 all_groups = " in any group",
                                 "")),
          mono  = TRUE
        ),
        stat_card(
          label = "Imputation",
          value = impute,
          trend = if (impute == "none") "NAs left visible"
                  else "expression matrix imputed"
        ),
        stat_card(
          label  = "Outlier method",
          value  = qc_method_labels(summary$outliers$method),
          trend  = sprintf("threshold = %g",
                           bundle$params$outlier_sd_threshold %||% 3),
          accent = "ok"
        )
      )
    })

    # Which metadata column colours the samples. It used to be a column
    # literally called `group` or nothing -- so the app's own template,
    # whose column is `condition`, drew an uncoloured PCA -- and the
    # legend under the plot was four fixed CSS colours that stopped
    # matching the points from the fifth group on. Now the column is the
    # user's choice (the Differential view's best guess by default) and
    # the legend is the plot's own.
    pca_color_choices <- shiny::reactive({
      bundle <- last_bundle()
      shiny::req(bundle)
      meta <- bundle$results$cleaned_input$meta_df
      if (is.null(meta) || !ncol(meta)) return(character(0))
      cands <- grouping_candidates(meta)
      design <- tryCatch(omicsCore::study_design(bundle$results$cleaned_input),
                         error = function(e) NULL)
      if (is.null(design)) {
        design <- tryCatch(omicsCore::study_design(active()$input),
                           error = function(e) NULL)
      }
      if (!is.null(design) && design$group_col %in% names(meta)) {
        cands <- c(design$group_col, setdiff(cands, design$group_col))
      }
      extra <- setdiff(names(meta)[vapply(meta, function(x) {
        n <- length(unique(stats::na.omit(x)))
        n >= 2L && n < nrow(meta)
      }, logical(1))], cands)
      c(cands, extra)
    })

    output$pca_color_picker <- shiny::renderUI({
      ch <- pca_color_choices()
      if (!length(ch)) return(NULL)
      sel <- shiny::isolate(input$pca_color_by)
      if (is.null(sel) || !sel %in% c(ch, "(none)")) sel <- ch[[1L]]
      htmltools::tags$div(
        class = "inline-control",
        shiny::selectInput(session$ns("pca_color_by"), label = "Colour by",
                           choices = c(ch, "(none)"), selected = sel,
                           width = "220px"))
    })

    output$pca <- shiny::renderPlot(res = PLOT_RES, alt = "Principal component plot of the samples", fit_to_width("pca", {
      bundle <- last_bundle()
      shiny::req(bundle)
      ch <- pca_color_choices()
      color_by <- input$pca_color_by
      if (is.null(color_by) || !color_by %in% ch) {
        color_by <- if (length(ch)) ch[[1L]] else NULL
      }
      if (identical(input$pca_color_by, "(none)")) color_by <- NULL
      p <- omicsCore::plot_qc(bundle, view = "pca", color_by = color_by)
      p + ggplot2::theme(legend.position = "bottom")
    }))

    # Which quality panel this modality is actually asking about.
    #
    # Missingness is the proteomics question: a peptide that was not
    # detected is a hole in the matrix. A counts matrix has no holes --
    # every gene has a number for every sample, most of them zero -- so
    # the panel reported "63,241 features, all at 0%", which is true and
    # says nothing, in the space that should have been showing whether a
    # library was under-sequenced.
    #
    # A default per modality, not a lock: an intensity matrix has a
    # meaningful total too, and someone with an imputed counts matrix may
    # well want the missingness view.
    default_quality_view <- shiny::reactive({
      if (identical(active()$input$omics_type %||% "", "rnaseq")) "depth"
      else "missing"
    })
    quality_view <- shiny::reactive({
      v <- input$quality_view
      if (is.null(v) || !v %in% c("missing", "depth")) default_quality_view()
      else v
    })

    output$quality_title <- shiny::renderUI({
      depth <- identical(quality_view(), "depth")
      htmltools::tagList(
        htmltools::tags$h3(class = "card-title",
                           if (depth) "Depth" else "Missingness"),
        htmltools::tags$span(
          class = "card-sub",
          if (depth) "library size and features detected"
          else "per-sample and per-feature missing rate")
      )
    })

    output$quality_picker <- shiny::renderUI({
      # isolate(), for the reason the layer picker documents: this output
      # writes input$quality_view and quality_view() reads it.
      sel <- shiny::isolate(input$quality_view)
      if (is.null(sel) || !sel %in% c("missing", "depth")) {
        sel <- default_quality_view()
      }
      shiny::radioButtons(
        session$ns("quality_view"), label = NULL,
        choices = c("Depth" = "depth", "Missingness" = "missing"),
        selected = sel, inline = TRUE)
    })

    output$missing <- shiny::renderPlot(res = PLOT_RES, alt = "Missing values per sample and per feature", fit_to_width("missing", {
      bundle <- last_bundle()
      shiny::req(bundle)
      omicsCore::plot_qc(bundle, view = quality_view())
    }))

    output$missing_caption <- shiny::renderUI({
      a <- active()
      bundle <- last_bundle()
      shiny::req(bundle)
      caption <- if (identical(quality_view(), "depth")) {
        d <- bundle$results$qc_summary$depth
        if (is.null(d) || nrow(d) == 0L) {
          "No depth summary for this layer."
        } else {
          low <- omicsCore::qc_depth_outliers(d)
          sprintf("%d samples \u00b7 median library %s \u00b7 %s",
                  nrow(d),
                  format(round(stats::median(d$library_size)), big.mark = ","),
                  if (length(low) == 0L) "none shallow"
                  else sprintf("shallow: %s", paste(low, collapse = ", ")))
        }
      } else if (a$is_demo) {
        "Demo fixture: ~5% of cells set to NA at random."
      } else {
        raw <- a$input$expr_mat
        raw_pct <- 100 * mean(is.na(raw))
        n_na <- sum(is.na(bundle$results$cleaned_input$expr_mat))
        n_cells <- length(bundle$results$cleaned_input$expr_mat)
        if (!is.null(bundle$results$qc_summary$imputation)) {
          sprintf("Imported layer: %.1f%% of cells missing; after imputation for this view: %.1f%%.",
                  raw_pct, 100 * n_na / max(n_cells, 1L))
        } else {
          sprintf("Imported layer: %d / %d cells missing (%.1f%%).",
                  sum(is.na(raw)), length(raw), raw_pct)
        }
      }
      htmltools::tags$div(
        class = "muted",
        style = "font-size:12px;margin-top:6px",
        caption
      )
    })

    # Expose the QC bundle for slice 3F (report).
    shiny::reactive(last_bundle())
  })
}

# ---- internal helpers ------------------------------------------------

qc_controls_card <- function(ns) {
  bslib::card(
    bslib::card_header(
      htmltools::tags$h3(class = "card-title", "Filters"),
      htmltools::tags$span(
        class = "card-sub",
        "which layer, and how it is filtered"
      )
    ),
    bslib::card_body(
      # A project holds several layers and QC describes exactly one of
      # them. Which one was decided elsewhere -- by the arrow in the
      # Projects table, or by falling back to the first proteomics layer
      # -- so the answer to "what am I looking at" was not on this page.
      shiny::uiOutput(ns("ui_layer")),
      htmltools::tags$div(
        class = "row-grid r-6-6",
        shiny::sliderInput(
          ns("missing_threshold"),
          label = htmltools::tagList(
            "Feature missing-rate cutoff",
            info_tip("Features missing in more than this fraction of samples (or of each group's samples, when the cutoff is applied by group) are filtered out of the QC view.")),
          min   = 0,
          max   = 1,
          value = 0.5,
          step  = 0.05
        ),
        shiny::radioButtons(
          ns("outlier_method"),
          label   = htmltools::tagList(
            "Outlier detection",
            info_tip(paste(
              "How samples are flagged: IQR of per-sample summaries, distance in PCA space,",
              "low connectivity (correlation) to the other samples, or leave-one-out:",
              "each sample set aside in turn and compared with how closely the remaining",
              "samples resemble each other. Only leave-one-out can flag a sample in a",
              "study of ten or fewer."))),
          choices = c("All four" = "all",
                      "IQR" = "iqr",
                      "PCA" = "pca",
                      "Connectivity" = "connectivity",
                      "Leave-one-out" = "loo"),
          selected = "all",
          inline   = TRUE
        )
      ),
      # Shown only when the layer has a column that splits its samples
      # into groups; rendered from the server for that reason.
      shiny::uiOutput(ns("ui_missing_filter")),
      # Proteomics only, and rendered from the server because the choices
      # depend on the layer and on which optional packages are installed.
      shiny::uiOutput(ns("ui_impute")),
      # Said here because the controls look as though they reach the
      # analysis: they do not. Differential and Integration read the
      # imported matrix, and a user who filtered here and then read a
      # volcano would otherwise assume the filter was applied.
      htmltools::tags$div(
        class = "muted", style = "font-size:11.5px;margin-top:6px",
        bsicons::bs_icon("info-circle"),
        " These settings shape the QC view only. The Differential and",
        " Integration views analyse the imported matrix."
      )
    )
  )
}

# What "All" runs in the outlier control.
QC_ALL_OUTLIER_METHODS <- c("pca", "connectivity", "iqr", "loo")
QC_MISSING_FILTERS <- c("global", "any_group", "all_groups")

# A QC request's settings in one canonical form, whether they came from
# the controls or from a saved bundle's params, so that the two compare
# equal when they ask for the same run.
qc_request_params <- function(thr, outlier, impute, missing_filter, group_col) {
  list(missing_threshold = as.numeric(thr),
       outlier_method = as.character(outlier),
       impute_method = impute,
       missing_filter = missing_filter %||% "global",
       group_col = if (!identical(missing_filter %||% "global", "global")) group_col)
}

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
  kept <- colnames(bundle$results$cleaned_input$expr_mat)
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

# The outlier methods as the control names them.
qc_method_labels <- function(methods) {
  labels <- c(pca = "PCA", connectivity = "connectivity", iqr = "IQR",
              loo = "leave-one-out", none = "none")
  out <- ifelse(methods %in% names(labels), labels[methods], methods)
  paste(out, collapse = " + ")
}

# The engine's notes about small studies, in the app's words: they speak
# of z-scores and thresholds, and the decision they inform -- whether a
# sample could have been flagged at all -- is simpler than that.
qc_plain_outlier_notes <- function(notes) {
  notes <- sub(paste0("^With (\\d+) samples a z-score cannot exceed [0-9.]+, so the PCA and ",
                      "connectivity tests flag nothing at a threshold of [0-9.e+-]+; ",
                      "the leave-one-out test still applies\\.$"),
               paste("With \\1 samples the PCA and connectivity checks cannot flag a sample.",
                     "The leave-one-out check can: it sets each sample aside in turn and asks",
                     "whether it resembles its closest sample much less than the other samples",
                     "resemble theirs."),
               notes)
  notes <- sub(paste0("^With (\\d+) samples a z-score cannot exceed [0-9.]+, so a threshold ",
                      "of [0-9.e+-]+ flags nothing; inspect the PCA plot instead\\.$"),
               paste("With \\1 samples this check cannot flag a sample. Look at the PCA plot,",
                     "or choose leave-one-out, which works with as few as four samples."),
               notes)
  sub("^The leave-one-out test needs at least (\\d+) samples; with (\\d+) it flags nothing\\.$",
      "The leave-one-out check needs at least \\1 samples, so with \\2 it flagged nothing.",
      notes)
}

# One sentence per sample the leave-one-out check flagged, saying what it
# saw in numbers a reader can check against the PCA.
qc_loo_explanation <- function(outliers) {
  st <- if (identical(outliers$method, "loo")) outliers$stats else outliers$by_method$loo$stats
  if (is.null(st) || !any(st$is_outlier)) return(NULL)
  st <- st[st$is_outlier, , drop = FALSE]
  htmltools::tags$p(lapply(seq_len(nrow(st)), function(i) {
    sprintf(paste("Leave-one-out: %s correlates %.3f with its closest sample (%s),",
                  "where the other samples typically reach %.3f with theirs. "),
            st$sample_id[i], st$nearest_correlation[i], st$nearest_sample[i],
            st$reference_correlation[i])
  }))
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

omics_display <- function(t) {
  switch(t %||% "",
         proteomics = "Proteomics",
         rnaseq     = "RNA-seq",
         "\u2014")
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# Matrix cells above which a QC run shows its steps: about 5,000
# features x 40 samples, which takes a second or more.
QC_PROGRESS_CELLS <- 2e5
