# The Filters card of the QC view.
#
# Split out of mod_qc_view.R. The server half is a plain function called
# from inside qc_view_server()'s moduleServer(), not a module, so the ids
# stay "qc-impute_method", "qc-missing_filter". It returns the reactives
# the module keeps, the request QC runs on among them.
qc_controls_server <- function(input, output, session, active, pick) {
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

  list(impute_choices = impute_choices, impute_applies = impute_applies,
       thr_r = thr_r, missing_group_choices = missing_group_choices,
       qc_request = qc_request)
}

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
