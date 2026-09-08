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
          shiny::plotOutput(ns("pca"), height = "360px"),
          shiny::uiOutput(ns("pca_legend"))
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
                           activation = shiny::reactiveVal(0L)) {
  shiny::moduleServer(id, function(input, output, session) {

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
      tag <- if (!is.null(want) && length(want) == 1L && want %in% names(exps)) {
        want
      } else {
        names(exps)[qc_default_layer_idx(exps)]
      }
      # baseline_input(), not effective_input(): QC is the one consumer that
      # must not read its own output. Re-running over already-filtered values
      # cannot return the features the previous run removed, so loosening the
      # cutoff would appear to do nothing, and switching the imputation method
      # would impute cells that are no longer missing. Each apply replaces the
      # layer's preprocessing rather than stacking onto it.
      base <- omicsCore::baseline_input(proj, tag)
      list(input = base$input, tag = tag, base_version = base$version,
           is_demo = FALSE)
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
      shiny::selectInput(session$ns("layer"), label = "Experiment layer",
                         choices = tags_avail, selected = sel)
    })

    # The QC bundle. tryCatch keeps a parameter mishap (e.g. a
    # slider that drops all features) from blowing up the view --
    # the notices panel surfaces the message and the previous
    # plots stay visible.
    last_bundle <- shiny::reactiveVal(NULL)
    last_error  <- shiny::reactiveVal(NULL)
    # The preview, with the identity it was computed under. `applied` is what
    # leaves this module; nothing downstream sees `last_bundle`.
    preview     <- shiny::reactiveVal(NULL)
    applied     <- shiny::reactiveVal(NULL)

    # The layer this result was computed on has been replaced, so the
    # result is no longer about anything in the project. NULL is the
    # module's own start-up state, so this only rewinds it.
    # A generation bump means the layer this view was describing has been
    # replaced. Applying QC bumps it too -- that is how the downstream
    # views are told to let go -- so the result of the apply must survive
    # the signal the apply itself raised. Without this the module cleared
    # `applied()` in the same flush, and the observer that files results
    # into the project found nothing: the cleaned data arrived, the QC
    # bundle and its diagnostics did not.
    just_applied <- shiny::reactiveVal(FALSE)
    shiny::observeEvent(invalidate(), {
      last_bundle(NULL)
      last_error(NULL)
      preview(NULL)
      if (isTRUE(just_applied())) {
        just_applied(FALSE)
      } else {
        applied(NULL)
      }
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

    # The mirror image of imputation: offered on raw counts and withheld
    # everywhere else.
    #
    # `min_count = 10` means ten reads, which is a low value on
    # `raw_count` and a high one on `tpm` or `logcpm` -- the same number
    # would filter almost nothing on one layer and almost everything on
    # another. So the control is not offered where it cannot be read the
    # same way, and run_qc() resolves those layers to off.
    low_expression_applies <- shiny::reactive({
      a <- active()
      inp <- if (a$is_demo) example_qc_input() else a$input
      identical(inp$omics_type %||% "", "rnaseq") &&
        identical(inp$assay_type %||% "", "raw_count")
    })

    output$ui_low_expression <- shiny::renderUI({
      if (!low_expression_applies()) {
        # Silence here is the wrong answer. On an RNA-seq layer the filter
        # is the control a reader is looking for, and its absence is
        # indistinguishable from the feature not existing -- there is
        # nothing on the page to check against. Say which layer this is and
        # why the filter does not apply to it.
        a <- active()
        inp <- if (a$is_demo) example_qc_input() else a$input
        if (!identical(inp$omics_type %||% "", "rnaseq")) return(NULL)
        return(htmltools::tags$div(
          class = "muted", style = "font-size:12px",
          sprintf(paste0("Low-expression filtering is offered on raw counts; ",
                         "this layer is %s. A count threshold means nothing ",
                         "on a normalised scale — re-import as raw_count ",
                         "to filter here."),
                  inp$assay_type %||% "unlabelled")
        ))
      }
      base <- shiny::isolate(applied_pp()$params) %||% list()
      htmltools::tagList(
        htmltools::tags$div(
          class = "row-grid r-6-6",
          shiny::numericInput(
            session$ns("min_count"),
            label = "Minimum count",
            value = base[["min_count"]] %||% 10,
            min = 0, step = 1),
          shiny::numericInput(
            session$ns("min_samples"),
            label = "…in at least this many samples",
            value = base[["min_samples"]] %||% 3,
            min = 0, step = 1)
        ),
        # Said out loud because the effect is large and invisible: on a
        # full human index this drops most of the rows, and someone who
        # reads "63,241 features" on import and "15,593" here should be
        # able to find out why without asking.
        htmltools::tags$div(
          class = "muted", style = "font-size:12px;margin-top:-8px",
          "Features detected in too few samples are dropped before ",
          "testing. Both DESeq2 and edgeR ask for this: an undetected ",
          "feature cannot support a contrast, and testing it anyway adds ",
          "to the multiple-testing correction. Set either to 0 to test ",
          "everything."
        )
      )
    })

    # The settings as run_qc() would resolve them, in one place, so the
    # preview's signature and the record kept on apply cannot disagree about
    # what was asked for.
    collect_qc_params <- function() {
      # The applied settings are the base; the three visible controls
      # override their own fields and nothing else.
      #
      # Built the other way round -- three controls, full stop -- re-applying
      # a restored project silently reset every parameter the interface does
      # not show. A layer saved with outlier_sd_threshold = 7 came back at
      # the default 3, and sample_missing_threshold vanished entirely, which
      # is the opposite of the decision that these keep their saved values
      # precisely because there is no control to set them with.
      base <- applied_pp()$params %||% list()
      base[["missing_threshold"]] <- input$missing_threshold %||%
        base[["missing_threshold"]] %||% 0.5
      base[["outlier_method"]] <- input$outlier_method %||%
        base[["outlier_method"]] %||% "iqr"
      # Read rather than trusted from the control: it is hidden when the
      # layer is not proteomics, and Shiny keeps an input's last value.
      base[["impute_method"]] <- if (!impute_applies()) {
        "none"
      } else {
        input$impute_method %||% base[["impute_method"]]
      }
      # Same reasoning, and the same trap: Shiny keeps a hidden input's
      # last value, so a count threshold set on an RNA-seq layer would
      # otherwise follow the user to a proteomics one and filter it on a
      # scale where 10 is not a low number.
      if (!low_expression_applies()) {
        base[["min_count"]] <- 0
        base[["min_samples"]] <- 0
      } else {
        # numericInput yields NA while the box is empty or mid-edit, and
        # NA reaches assert_number() as a number -- the preview failed
        # with a validation error on every keystroke that cleared the
        # field.
        num <- function(v, fallback) {
          if (is.null(v) || length(v) != 1L || is.na(v)) fallback else v
        }
        base[["min_count"]] <- num(input$min_count,
                                   base[["min_count"]] %||% 10)
        base[["min_samples"]] <- num(input$min_samples,
                                     base[["min_samples"]] %||% 3)
      }
      base
    }

    # What this preview is about: which project, which layer, which version
    # of it, and which settings. Captured when the preview runs, so applying
    # can only confirm it -- never mint a fresh identity for a stale result.
    preview_signature <- function(a, params) {
      list(project_id = current_project()$metadata$project_id,
           tag = a$tag, base_version = a$base_version, params = params)
    }

    shiny::observe({
      a <- active()
      # One list, used for the run and for the signature. Computing them
      # separately is how the preview came to be signed with the merged
      # parameters while being run with only the three visible ones --
      # every unexposed setting reset to its default on the next apply.
      params <- collect_qc_params()

      # The demo runs through run_qc() like a real project rather than
      # returning a fixed bundle. The controls are enabled, and a control
      # that is enabled and does nothing reads as a broken app; the demo
      # input is 50 x 12, so a re-run is milliseconds.
      qc_input <- if (a$is_demo) example_qc_input() else a$input

      bundle <- tryCatch(
        do.call(omicsCore::run_qc, c(list(input = qc_input), params)),
        error = function(e) e)

      # A failed preview keeps the previous plots on screen -- useful -- but
      # must not leave the previous *bundle* looking like the answer to the
      # current settings. It did, and the apply button then stamped the
      # current layer's identity onto a result computed from another one.
      sig <- if (a$is_demo) NULL else preview_signature(a, params)
      if (inherits(bundle, "error")) {
        last_error(conditionMessage(bundle))
        preview(list(bundle = NULL, sig = sig, ok = FALSE,
                     error = conditionMessage(bundle)))
      } else {
        last_error(NULL)
        last_bundle(bundle)
        preview(list(bundle = bundle, sig = sig, ok = TRUE, error = NULL))
      }
    })

    output$header <- shiny::renderUI({
      a <- active()
      bundle <- last_bundle()
      n_in  <- if (is.null(bundle)) NA_integer_ else bundle$input_info$n_samples_in
      n_feat <- if (is.null(bundle)) NA_integer_ else bundle$input_info$n_features_in
      omics <- if (is.null(bundle)) "\u2014" else omics_display(bundle$input_info$omics_type)
      # The assay type, said out loud. Two of the controls on this page are
      # offered or withheld on the strength of it -- imputation on
      # proteomics, the count filter on raw counts -- and it was the one
      # thing about the layer the page did not show, so a withheld control
      # looked like a missing feature with nothing to check.
      assay <- if (is.null(bundle)) NULL else bundle$input_info$assay_type
      view_header(
        title    = "Quality control",
        subtitle = htmltools::tagList(
          omics,
          if (!is.null(assay) && nzchar(assay)) htmltools::HTML(" &middot; "),
          if (!is.null(assay) && nzchar(assay)) assay,
          htmltools::HTML(" &middot; "),
          sprintf("%s features \u00B7 %s samples",
                  format(n_feat, big.mark = ","),
                  format(n_in,  big.mark = ",")),
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
      err <- last_error()
      if (is.null(err)) return(NULL)
      notice(
        title  = "QC could not run",
        detail = err,
        kind   = "warn"
      )
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
          label  = "Samples passing",
          value  = sprintf("%d / %d", info$n_samples_out, info$n_samples_in),
          trend  = if (n_flagged_samp == 0L) "no outliers flagged"
                   else sprintf("%d flagged (%s)",
                                n_flagged_samp,
                                summary$outliers$method),
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
          value  = summary$outliers$method,
          trend  = sprintf("threshold = %g",
                           bundle$params$outlier_sd_threshold %||% 3),
          accent = "ok"
        )
      )
    })

    output$pca <- shiny::renderPlot({
      bundle <- last_bundle()
      shiny::req(bundle)
      color_by <- if ("group" %in% names(bundle$results$cleaned_input$meta_df))
        "group" else NULL
      omicsCore::plot_qc(bundle, view = "pca", color_by = color_by)
    })

    output$pca_legend <- shiny::renderUI({
      bundle <- last_bundle()
      shiny::req(bundle)
      meta <- bundle$results$cleaned_input$meta_df
      if (!"group" %in% names(meta)) return(NULL)
      groups <- sort(unique(as.character(meta$group)))
      colors <- c("var(--brand-600)", "var(--omics-down)",
                  "var(--accent-500)", "var(--ok)")
      htmltools::tags$div(
        class = "legend",
        lapply(seq_along(groups), function(i) {
          legend_swatch(groups[i], colors[((i - 1L) %% length(colors)) + 1L])
        })
      )
    })

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

    output$missing <- shiny::renderPlot({
      bundle <- last_bundle()
      shiny::req(bundle)
      omicsCore::plot_qc(bundle, view = quality_view())
    })

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
        n_na <- sum(is.na(bundle$results$cleaned_input$expr_mat))
        n_cells <- length(bundle$results$cleaned_input$expr_mat)
        sprintf("Live layer: %d / %d cells missing (%.1f%%).",
                n_na, n_cells, 100 * n_na / max(n_cells, 1L))
      }
      htmltools::tags$div(
        class = "muted",
        style = "font-size:12px;margin-top:6px",
        caption
      )
    })

    # ---- apply ---------------------------------------------------------

    # The preprocessing already applied to the layer on screen, if any.
    applied_pp <- shiny::reactive({
      proj <- current_project()
      a <- active()
      if (is.null(proj) || isTRUE(a$is_demo)) return(NULL)
      pp <- proj$preprocessing[[a$tag]]
      # Preprocessing left behind by a layer replacement describes data the
      # project no longer holds; effective_input() already refuses to serve
      # it, and the controls must not offer it back either.
      if (is.null(pp) || !identical(pp$base_input_version, a$base_version)) {
        return(NULL)
      }
      pp
    })

    shiny::observeEvent(input$apply_to_downstream, {
      p <- preview()
      a <- active()
      if (isTRUE(a$is_demo)) {
        shiny::showNotification(
          "This is the built-in demo. Import a project to apply QC.",
          type = "message")
        return()
      }
      if (is.null(p) || !isTRUE(p$ok)) {
        shiny::showNotification(
          paste0("These settings did not run, so there is nothing to apply",
                 if (!is.null(p$error)) paste0(": ", p$error) else "."),
          type = "error", duration = NULL)
        return()
      }
      # The preview is only the answer to the settings that produced it. A
      # signature check here, rather than re-deriving the identity from the
      # current controls, is what stops a result computed before a layer
      # switch from being applied to the layer after it.
      if (!identical(p$sig, preview_signature(a, collect_qc_params()))) {
        shiny::showNotification(
          "The settings changed after this preview ran. Wait for it to finish, then apply.",
          type = "warning", duration = NULL)
        return()
      }
      b <- p$bundle
      b$input_info$layer_tag <- a$tag
      b$input_info$base_input_version <- a$base_version
      # QC reads the import, so its record refers to the baseline. Recorded
      # against the effective input it would be retired by its own result:
      # applying it is what makes the effective version differ.
      b <- bind_result(b, current_project(), "qc", tags = a$tag,
                       versions = a$base_version, stage = "baseline")
      just_applied(TRUE)
      applied(list(tag = a$tag, bundle = b, params = b$params,
                   base_version = a$base_version,
                   at = Sys.time()))
      shiny::showNotification(
        sprintf("QC applied to '%s': %d \u2192 %d features, %d \u2192 %d samples. Downstream analyses will use this.",
                a$tag,
                b$input_info$n_features_in, b$input_info$n_features_out,
                b$input_info$n_samples_in,  b$input_info$n_samples_out),
        type = "message", duration = 8)
    })

    output$applied_status <- shiny::renderUI({
      a <- active()
      if (isTRUE(a$is_demo)) return(NULL)
      pp <- applied_pp()
      txt <- if (is.null(pp)) {
        "Not applied \u2014 analyses are using the raw import."
      } else {
        sprintf("Applied: %d features, %d samples \u00B7 %s imputation \u00B7 cutoff %s",
                nrow(pp$cleaned$expr_mat), ncol(pp$cleaned$expr_mat),
                pp$params$impute_method %||% "none",
                format(pp$params$missing_threshold %||% NA))
      }
      htmltools::tags$div(class = "muted",
                          style = "font-size:12px;margin-top:6px", txt)
    })

    # ---- restoring what was applied -------------------------------------

    # On arriving at a layer -- a new project, a restored one, a switch of
    # tag -- the controls show what is applied to it, not whatever they were
    # last left on. Without this, opening a project saved with `min` while
    # `knn` happened to be selected showed `knn`, and the page then described
    # settings the analyses were not using.
    #
    # It runs on arrival only, so the user's edits afterwards are their own.
    # Keyed on the activation batch and the layer, not on the project
    # object. Every analysis write-back updates the project, so an observer
    # watching that fired on unrelated events -- finishing a differential,
    # or renaming the project -- and pushed the applied settings back over
    # whatever the user was in the middle of typing.
    shiny::observeEvent(list(activation(), active()$tag), {
      pp <- applied_pp()
      params <- pp$params
      if (is.null(params)) return()
      if (!is.null(params$missing_threshold)) {
        shiny::updateSliderInput(session, "missing_threshold",
                                 value = params$missing_threshold)
      }
      if (!is.null(params$outlier_method)) {
        shiny::updateRadioButtons(session, "outlier_method",
                                  selected = params$outlier_method)
      }
      if (!is.null(params$impute_method) && impute_applies()) {
        shiny::updateSelectInput(session, "impute_method",
                                 selected = params$impute_method)
      }
      if (low_expression_applies()) {
        if (!is.null(params$min_count)) {
          shiny::updateNumericInput(session, "min_count",
                                    value = params$min_count)
        }
        if (!is.null(params$min_samples)) {
          shiny::updateNumericInput(session, "min_samples",
                                    value = params$min_samples)
        }
      }
    }, priority = 25L, ignoreNULL = FALSE)

    # What leaves the module is the applied preprocessing, never the preview.
    shiny::reactive(applied())
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
          label = "Feature missing-rate cutoff",
          min   = 0,
          max   = 1,
          value = 0.5,
          step  = 0.05
        ),
        shiny::radioButtons(
          ns("outlier_method"),
          label   = "Outlier detection",
          choices = c("IQR" = "iqr",
                      "PCA" = "pca",
                      "Connectivity" = "connectivity"),
          selected = "iqr",
          inline   = TRUE
        )
      ),
      # Proteomics only, and rendered from the server because the choices
      # depend on the layer and on which optional packages are installed.
      shiny::uiOutput(ns("ui_impute")),
      # Raw counts only, for the mirror-image reason.
      shiny::uiOutput(ns("ui_low_expression")),
      # Everything above is a preview. Until this is pressed the panels
      # describe what these settings would do, and the analyses keep using
      # whatever was applied last -- which for a fresh layer is the import.
      #
      # The separation exists because the preview runs on every control
      # change, from the moment the view is constructed, whether or not the
      # user has opened it. Feeding that straight to the analyses would mean
      # the default cutoff silently filtered every project.
      htmltools::tags$div(
        class = "qc-apply",
        shiny::actionButton(ns("apply_to_downstream"),
                            "Apply to downstream analyses",
                            class = "btn btn-primary"),
        shiny::uiOutput(ns("applied_status"))
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
