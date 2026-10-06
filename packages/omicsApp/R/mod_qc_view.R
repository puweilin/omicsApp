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
#' way.
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
      want <- input$layer
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
      sel <- shiny::isolate(input$layer)
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

    # The layer this result was computed on has been replaced, so the
    # result is no longer about anything in the project. NULL is the
    # module's own start-up state, so this only rewinds it.
    shiny::observeEvent(invalidate(), {
      last_bundle(NULL)
      last_error(NULL)
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
      sel <- shiny::isolate(input$impute_method)
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
    thr_in <- shiny::reactive(input$missing_threshold %||% 0.5)
    thr_r <- if (isTRUE(qc_delay > 0)) shiny::debounce(thr_in, qc_delay) else thr_in

    shiny::observe({
      a <- active()
      thr   <- thr_r()
      out_m <- input$outlier_method   %||% "all"
      # All three by default: after vsn the per-sample means are equal,
      # and the IQR test on them -- the old default -- could not see a
      # sample that PCA and connectivity both flagged.
      if (identical(out_m, "all")) out_m <- c("pca", "connectivity", "iqr")
      # Read here rather than trusted from the input: the control is
      # hidden when the layer is not proteomics, but Shiny keeps an
      # input's last value, so switching from a proteomics layer with
      # `knn` selected to a counts layer would otherwise impute counts
      # with a control the user can no longer see.
      # NULL lets run_qc() resolve it per modality, which is where that
      # decision belongs -- and means the control and a plain
      # run_qc(input) agree instead of quietly differing.
      imp <- if (!impute_applies()) "none" else input$impute_method

      # The demo runs through run_qc() like a real project rather than
      # returning a fixed bundle. Both controls above are enabled, and
      # a control that is enabled and does nothing reads as a broken
      # app; the demo input is 50 x 12, so a re-run is milliseconds.
      qc_input <- if (a$is_demo) example_qc_input() else a$input

      run <- function() omicsCore::run_qc(
        qc_input,
        missing_threshold = thr,
        outlier_method    = out_m,
        impute_method     = imp
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
        last_error(NULL)
        last_bundle(bundle)
      }
    })

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
      out <- htmltools::tagList()
      if (length(flagged) && !active()$is_demo) {
        out <- htmltools::tagAppendChild(out, notice(
          title = sprintf("Possible outlier%s: %s", if (length(flagged) > 1L) "s" else "",
                          paste(flagged, collapse = ", ")),
          detail = htmltools::tagList(
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
      proj$bundles <- drop_layer_bundles(proj$bundles, tag)
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
                                paste(summary$outliers$method, collapse = " + ")),
          accent = if (n_flagged_samp == 0L) "ok" else "warn",
          mono   = TRUE
        ),
        stat_card(
          label = "Features kept",
          value = format(info$n_features_out, big.mark = ","),
          trend = sprintf("%d filtered at %.0f%% missing",
                          n_flagged_feat,
                          100 * (bundle$params$missing_threshold %||% 0.5)),
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
          value  = paste(summary$outliers$method, collapse = " + "),
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
            info_tip("Features missing in more than this fraction of samples are filtered out of the QC view.")),
          min   = 0,
          max   = 1,
          value = 0.5,
          step  = 0.05
        ),
        shiny::radioButtons(
          ns("outlier_method"),
          label   = htmltools::tagList(
            "Outlier detection",
            info_tip("How samples are flagged: IQR of per-sample summaries, distance in PCA space, or low connectivity (correlation) to the other samples.")),
          choices = c("All three" = "all",
                      "IQR" = "iqr",
                      "PCA" = "pca",
                      "Connectivity" = "connectivity"),
          selected = "all",
          inline   = TRUE
        )
      ),
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
