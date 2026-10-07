# The PCA and quality panels of the QC view.
#
# Split out of mod_qc_view.R and called from inside its module server,
# so the ids stay "qc-pca", "qc-quality_view". Returns the colouring
# choices and the quality panel shown.
qc_plots_server <- function(input, output, session, active, last_bundle) {
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
    meta <- qc_bundle_meta(bundle)
    if (is.null(meta) || !ncol(meta)) return(character(0))
    cands <- grouping_candidates(meta)
    # Results saved before run_qc() stopped storing the cleaned input
    # carry its design; otherwise the layer's is the one.
    design <- if (!is.null(bundle$results$cleaned_input)) {
      tryCatch(omicsCore::study_design(bundle$results$cleaned_input),
               error = function(e) NULL)
    }
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
      n_na <- qc_missing_after(bundle)
      n_cells <- bundle$input_info$n_features_out * bundle$input_info$n_samples_out
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

  list(pca_color_choices = pca_color_choices, quality_view = quality_view)
}

# Cells still missing after imputation: the ones the method could not fill.
qc_missing_after <- function(bundle) {
  if (!is.null(bundle$results$cleaning)) {
    sum(is.na(bundle$results$cleaning$imputed_values))
  } else {
    sum(is.na(bundle$results$cleaned_input$expr_mat))
  }
}
