# The header, the notices and the stat cards of the QC view, and the
# exclusion of flagged samples.
#
# Split out of mod_qc_view.R and called from inside its module server,
# so the ids stay "qc-exclude_flagged", "qc-stats".
qc_results_server <- function(input, output, session, navigate, current_project, active,
                              last_bundle, last_error) {
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
    # The method that ran: MinProb falls back to MinDet when the data
    # cannot support it, and the notices say why.
    impute_used <- summary$imputation$method %||% impute
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
        value = impute_used,
        trend = if (impute == "none") "NAs left visible"
                else if (!identical(impute_used, impute))
                  sprintf("in place of %s (see the note above)", impute)
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
  invisible()
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
