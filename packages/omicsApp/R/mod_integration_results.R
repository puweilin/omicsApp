# What the integration view shows: header, notices, stat cards, plots and
# tables, for a live result or the demo.
#
# Split out of mod_integration_view.R and called from inside its module
# server, so the outputs keep their ids ("integration-scatter",
# "integration-top_table").
integration_results_server <- function(input, output, session, navigate, method, can_run,
                                       is_demo, integration_bundle, integration_error,
                                       running) {
  # Demo and live share one dispatcher.
  plot_bundle <- shiny::reactive({
    if (isTRUE(is_demo())) example_integration_bundle()
    else integration_bundle()
  })

  shown_method <- shiny::reactive({
    b <- plot_bundle()
    if (is.null(b)) method() else b$params$method %||% method()
  })

  conc_df <- shiny::reactive({
    b <- plot_bundle()
    shiny::req(b)
    b$results$integration_df
  })

  output$header <- shiny::renderUI({
    info <- can_run()
    subtitle <- if (isTRUE(is_demo())) {
      "Proteomics \u00D7 RNA-seq \u00B7 demo data \u00B7 60 paired features"
    } else if (isTRUE(info$ok)) {
      what <- switch(method(),
                     correlation = "sample-level correlation",
                     active_pathways = "ActivePathways",
                     gsub("_vs_", " vs ", info$comparison %||%
                            paste0(info$case, "_vs_", info$control)))
      sprintf("%s \u00D7 %s \u00B7 %s", info$primary_tag,
              info$secondary_tag, what)
    } else {
      "not run yet"
    }
    view_header(
      title    = "Multi-omics integration",
      subtitle = subtitle,
      actions  = shiny::actionButton(
        session$ns("rerun"),
        if (is.null(integration_bundle())) "Run integration" else "Re-run integration",
        icon = shiny::icon("play"),
        class = "btn btn-primary"
      )
    )
  })

  output$notices <- shiny::renderUI({
    tagged <- htmltools::tagList()
    err <- integration_error()
    if (identical(err, CANCELLED_MESSAGE)) {
      tagged <- htmltools::tagAppendChild(tagged, notice(
        "Cancelled", "Press Run integration to start again.", kind = "info"))
    } else if (!is.null(err)) {
      tagged <- htmltools::tagAppendChild(
        tagged,
        notice(title  = "The integration could not be computed",
               detail = integration_error_hint(err),
               kind   = "error",
               technical = err)
      )
    }
    if (isTRUE(running())) {
      tagged <- htmltools::tagAppendChild(
        tagged, notice("Running\u2026", kind = "info"))
    }
    info <- can_run()
    if (isTRUE(is_demo())) {
      tagged <- htmltools::tagAppendChild(tagged, notice(
        title  = "Showing demo data",
        detail = paste("Import two omics layers (or load the example",
                       "project from the Project view) to integrate your own."),
        kind   = "info"))
    } else if (!isTRUE(info$ok)) {
      reason <- switch(info$reason %||% "",
        layers = "Integration needs two layers in the project. Import a second omics layer.",
        diff   = paste("Fold-change concordance compares two differential results.",
                       "Run a differential analysis first, or switch to sample-level correlation."),
        info$detail %||% "The prerequisites are not met.")
      tagged <- htmltools::tagAppendChild(tagged, notice(
        title = "Nothing to integrate yet", detail = reason, kind = "info"))
      if (identical(info$reason, "diff") && is.function(navigate)) {
        tagged <- htmltools::tagAppendChild(tagged, htmltools::tags$div(
          style = "margin:6px 0 10px",
          shiny::actionButton(session$ns("go_diff"), "Go to Differential \u2192",
                              class = "btn btn-sm btn-ghost")))
      }
    } else if (length(info$dropped_covs)) {
      tagged <- htmltools::tagAppendChild(tagged, notice(
        title = "Covariates not repeated on the partner layer",
        detail = sprintf("'%s' has no %s column, so its differential run is not adjusted for it.",
                         info$secondary_tag,
                         paste(sprintf("'%s'", info$dropped_covs), collapse = ", ")),
        kind = "warn"))
    }
    b <- integration_bundle()
    if (!is.null(b) && length(b$warnings)) {
      for (w in b$warnings) {
        tagged <- htmltools::tagAppendChild(tagged, notice(w, kind = "warn"))
      }
    }
    tagged
  })

  if (is.function(navigate)) {
    shiny::observeEvent(input$go_diff, navigate("diff"))
  }

  output$stats <- shiny::renderUI({
    b <- plot_bundle()
    shiny::req(b)
    df <- b$results$integration_df
    shiny::req(nrow(df) > 0L)
    integration_stat_cards(df, b)
  })

  output$results <- shiny::renderUI({
    m <- shown_method()
    if (!isTRUE(is_demo()) && is.null(integration_bundle())) return(NULL)
    ns <- session$ns
    if (identical(m, "correlation")) {
      htmltools::tagList(
        htmltools::tags$div(class = "row-grid r-6-6",
          integration_plot_card(ns("cor_scatter"),
                                "Correlation per gene",
                                "r across paired samples \u00B7 y = -log10 adjusted p",
                                plot_download_ui(ns("cor_scatter_download"))),
          integration_table_card(ns("top_table"), "Top features",
                                 "ranked by adjusted p")))
    } else if (identical(m, "active_pathways")) {
      htmltools::tagList(
        htmltools::tags$div(class = "row-grid r-6-6",
          integration_plot_card(ns("ap_dot"), "ActivePathways",
                                "colour = direction \u00B7 shape = which layers found it",
                                plot_download_ui(ns("ap_dot_download"))),
          integration_table_card(ns("top_table"), "Pathways",
                                 "ranked by adjusted p")))
    } else {
      htmltools::tagList(
        htmltools::tags$div(class = "row-grid r-6-6",
          integration_plot_card(ns("scatter"), "Fold-change concordance",
                                "A vs B \u00B7 dashed = agreement \u00B7 coloured = hit in both",
                                plot_download_ui(ns("scatter_download"))),
          integration_plot_card(ns("top_hits"), "Top hits in both layers",
                                "one row per feature \u00B7 one dot per layer",
                                plot_download_ui(ns("top_hits_download")))),
        htmltools::tags$div(class = "row-grid r-6-6",
          integration_table_card(ns("top_table"), "Features",
                                 "hits in both layers first"),
          if (isTRUE(is_demo())) integration_ap_card(ns)))
    }
  })

  # Each figure's ggplot, kept apart from its renderPlot() so the hover
  # read-out (and anything else that needs the figure) uses the same one.
  scatter_plot <- shiny::reactive({
    b <- plot_bundle()
    shiny::req(b, identical(b$params$method, "concordance"))
    # Six names fit a half-width card without the labels piling up.
    omicsCore::plot_integration(b, view = "effect_pair", top_n = 6L)
  })

  # Twelve rows is what a 320 px card holds at a readable size.
  top_hits_plot <- shiny::reactive({
    b <- plot_bundle()
    shiny::req(b, identical(b$params$method, "concordance"))
    omicsCore::plot_integration(b, view = "top_hits", top_n = 12L)
  })

  cor_scatter_plot <- shiny::reactive({
    b <- plot_bundle()
    shiny::req(b, identical(b$params$method, "correlation"))
    omicsCore::plot_integration(b, view = "scatter")
  })

  ap_dot_plot <- shiny::reactive({
    b <- plot_bundle()
    shiny::req(b, identical(b$params$method, "active_pathways"))
    omicsCore::plot_integration(b, view = "dotplot")
  })

  output$scatter <- shiny::renderPlot(
    fit_to_width("scatter", scatter_plot()), res = PLOT_RES,
    alt = paste("Each feature's effect in one layer against its effect in the other,",
                "the hits in both layers coloured and the top ones named"))
  output$top_hits <- shiny::renderPlot(
    fit_to_width("top_hits", top_hits_plot()), res = PLOT_RES,
    alt = "The top hits in both layers, one row each, with a dot for each layer's effect")
  output$cor_scatter <- shiny::renderPlot(
    fit_to_width("cor_scatter", cor_scatter_plot()), res = PLOT_RES,
    alt = "Per-feature correlation between the layers across paired samples")
  output$ap_dot <- shiny::renderPlot(
    fit_to_width("ap_dot", ap_dot_plot()), res = PLOT_RES,
    alt = "Pathways found by combining the two layers")

  # Only the top few points are named on the plots; hovering (or
  # tapping) any other names it and gives its numbers.
  plot_hover_server("scatter", scatter_plot, effect_pair_hover_text, input, output, session)
  plot_hover_server("top_hits", top_hits_plot, top_hits_hover_text, input, output, session)
  plot_hover_server("cor_scatter", cor_scatter_plot, correlation_hover_text, input, output, session)
  plot_hover_server("ap_dot", ap_dot_plot, active_pathways_hover_text, input, output, session)

  output$top_table <- DT::renderDT({
    b <- plot_bundle()
    shiny::req(b)
    out <- integration_result_table(b$results$integration_df,
                                    b$params$method %||% "concordance",
                                    b$params$experiments)
    DT::datatable(out, rownames = FALSE, selection = "single",
                  options = list(pageLength = 10, dom = "ftip",
                                 scrollX = TRUE))
  }, server = TRUE)

  # The pathway fixture is only ever part of the demo; a live
  # ActivePathways run draws from its own result.
  ap_df <- shiny::reactive(example_integration_tables()$active_pathways_df)

  output$ap_table <- DT::renderDT({
    ap <- ap_df()
    out <- data.frame(
      Pathway          = ap$pathway_name,
      `p (A)`          = signif(ap$p_a, 3),
      `p (B)`          = signif(ap$p_b, 3),
      `p (combined)`   = signif(ap$p_combined, 3),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    DT::datatable(
      out,
      rownames  = FALSE,
      selection = "single",
      options   = list(
        pageLength = 10,
        dom        = "tip",
        scrollX    = TRUE,
        columnDefs = list(list(className = "dt-right",
                               targets = c(1, 2, 3)))
      )
    )
  }, server = TRUE)
  # The figures, for whatever else needs the ggplot a card draws.
  invisible(list(scatter = scatter_plot, top_hits = top_hits_plot,
                 cor_scatter = cor_scatter_plot, ap_dot = ap_dot_plot))
}

# A plain-language reading of the errors a run most often ends in.
integration_error_hint <- function(msg) {
  msg <- msg %||% ""
  if (grepl("No shared .* features", msg)) {
    return(paste("The two layers have no gene symbols in common. Check that both",
                 "carry a symbol column (feature_symbol) from the same organism,",
                 "or match them with your own table under Feature matching."))
  }
  if (grepl("paired by the feature link|feature link has no row", msg)) {
    return(paste("The mapping table matches no features of these two layers.",
                 "Check which column holds which layer's identifiers, under",
                 "Feature matching."))
  }
  if (grepl("guess", msg, fixed = TRUE)) {
    return("Accept the suggested sample pairing below, or add a donor column to both layers.")
  }
  if (grepl("No sample pairing|paired samples", msg)) {
    return("Sample-level correlation needs samples paired across the layers; see the pairing card below.")
  }
  if (grepl("ActivePathways", msg, fixed = TRUE)) {
    return("The ActivePathways package is not installed on this server.")
  }
  "See the technical details below."
}

integration_stat_cards <- function(df, bundle) {
  method <- bundle$params$method %||% "concordance"
  exps <- bundle$params$experiments %||% c("A", "B")
  if (identical(method, "correlation")) {
    sig <- df$is_significant %in% TRUE
    return(htmltools::tags$div(
      class = "stat-grid",
      stat_card("Genes correlated", format(nrow(df), big.mark = ","),
                trend = sprintf("%s \u2194 %s", exps[1], exps[2]),
                accent = "brand", mono = TRUE),
      stat_card("Paired samples", bundle$params$method_info$n_samples %||% "\u2014",
                trend = sprintf("pairing: %s",
                                bundle$params$method_info$pairing_source %||% "\u2014"),
                mono = TRUE),
      stat_card("Positive (sig)", sum(sig & df$effect > 0, na.rm = TRUE),
                trend = "RNA and protein rise together", accent = "up", mono = TRUE),
      stat_card("Median r", sprintf("%.2f", stats::median(df$effect, na.rm = TRUE)),
                trend = "all genes", mono = TRUE)
    ))
  }
  if (identical(method, "active_pathways") && "evidence" %in% names(df)) {
    # Which way the significant pathways went. The test expected the two
    # layers to change the same way, so "up in both" and "down in both"
    # are the findings, and a pathway whose layers disagree is the
    # exception worth a look.
    sig <- df$is_significant %in% TRUE
    return(htmltools::tags$div(
      class = "stat-grid",
      stat_card("Significant pathways", sum(sig),
                trend = sprintf("of %d tested \u00B7 %d only when the layers are merged",
                                nrow(df), sum(sig & df$evidence %in% "combined")),
                accent = "brand", mono = TRUE),
      stat_card("Up in both layers", sum(sig & df$direction %in% "up"),
                trend = "the driving genes rise in both", accent = "up", mono = TRUE),
      stat_card("Down in both layers", sum(sig & df$direction %in% "down"),
                trend = "the driving genes fall in both", accent = "down", mono = TRUE),
      stat_card("Layers disagree", sum(sig & df$direction %in% "mixed"),
                trend = "mixed, or opposite directions in the two layers",
                accent = if (any(sig & df$direction %in% "mixed")) "warn" else "ok",
                mono = TRUE)
    ))
  }
  if (identical(method, "active_pathways")) {
    sig <- df$is_significant %in% TRUE
    return(htmltools::tags$div(
      class = "stat-grid",
      stat_card("Pathways tested", nrow(df), accent = "brand", mono = TRUE),
      stat_card("Significant", sum(sig), mono = TRUE),
      stat_card("Shared by both", sum(sig & df$direction %in% "shared"),
                accent = "up", mono = TRUE),
      stat_card("Combined only", sum(sig & df$direction %in% "combined"),
                trend = "found only when the layers are merged", mono = TRUE)
    ))
  }
  # Concordance. The counts are of features that are hits in *both*
  # layers: every feature has some sign pair, so counting quadrants over
  # all of them reported half of an unrelated background as "concordant".
  both <- if (all(c("significant_a", "significant_b") %in% names(df))) {
    df$significant_a %in% TRUE & df$significant_b %in% TRUE
  } else {
    df$is_significant %in% TRUE
  }
  quad <- df$quadrant
  up_n <- sum(both & quad %in% "up_up")
  down_n <- sum(both & quad %in% "down_down")
  disc_n <- sum(both & quad %in% c("up_down", "down_up"))
  rho <- if (all(c("effect_a", "effect_b") %in% names(df)) && sum(both) >= 3L) {
    suppressWarnings(stats::cor(df$effect_a[both], df$effect_b[both],
                                method = "spearman", use = "pairwise.complete.obs"))
  } else NA_real_
  htmltools::tags$div(
    class = "stat-grid",
    stat_card(
      label  = "Paired features",
      value  = format(nrow(df), big.mark = ","),
      trend  = sprintf("%s \u2194 %s \u00B7 %d hits in both", exps[1], exps[2], sum(both)),
      accent = "brand", mono = TRUE
    ),
    stat_card(label = "Concordant \u2191", value = up_n,
              trend = "hit and up in both layers", accent = "up", mono = TRUE),
    stat_card(label = "Concordant \u2193", value = down_n,
              trend = "hit and down in both layers", accent = "down", mono = TRUE),
    stat_card(
      label = "Discordant",
      value = disc_n,
      trend = if (is.na(rho)) "hit in both, opposite signs"
              else sprintf("opposite signs \u00B7 Spearman \u03C1 of hits %.2f", rho),
      accent = if (disc_n > 0L) "warn" else "ok",
      mono = TRUE
    )
  )
}

integration_result_table <- function(df, method, experiments) {
  exps <- experiments %||% c("A", "B")
  if (identical(method, "concordance") &&
      all(c("effect_a", "effect_b") %in% names(df))) {
    both <- if (all(c("significant_a", "significant_b") %in% names(df))) {
      df$significant_a %in% TRUE & df$significant_b %in% TRUE
    } else df$is_significant %in% TRUE
    ord <- order(!both, df$p_value, na.last = TRUE)
    d <- df[ord, , drop = FALSE]
    out <- data.frame(
      Feature = feature_row_label(d),
      a = round(d$effect_a, 3),
      b = round(d$effect_b, 3),
      Quadrant = d$quadrant,
      `Hit in both` = ifelse(both[ord], "yes", ""),
      `Combined p` = signif(d$p_value, 3),
      check.names = FALSE, stringsAsFactors = FALSE)
    names(out)[2:3] <- paste0("Effect (", exps[1:2], ")")
    return(out)
  }
  if (identical(method, "active_pathways") && "evidence" %in% names(df)) {
    d <- df[order(df$adj_p_value, na.last = TRUE), , drop = FALSE]
    word <- function(x) unname(ifelse(is.na(x), "\u2014", x))
    out <- data.frame(
      Pathway = d$feature_symbol,
      `Adj. p` = signif(d$adj_p_value, 3),
      Direction = word(c(up = "up in both layers", down = "down in both layers",
                         mixed = "mixed / layers disagree")[d$direction]),
      a = word(d$direction_a),
      b = word(d$direction_b),
      `Found by` = word(c(shared = "both layers", unique = "one layer",
                          combined = "only combined")[d$evidence]),
      check.names = FALSE, stringsAsFactors = FALSE)
    names(out)[4:5] <- paste0("In ", exps[1:2])
    return(out)
  }
  d <- df[order(df$adj_p_value, na.last = TRUE), , drop = FALSE]
  data.frame(
    Feature = feature_row_label(d),
    Effect = round(d$effect, 3),
    `Adj. p` = signif(d$adj_p_value, 3),
    Direction = d$direction,
    check.names = FALSE, stringsAsFactors = FALSE)
}

# A result row's name: the gene, or -- where a gene has several pairs --
# the pair's id ("TP53 (P04637-2)"), so its rows can be told apart.
feature_row_label <- function(d) {
  sym <- d$feature_symbol
  shared <- !is.na(sym) & (duplicated(sym) | duplicated(sym, fromLast = TRUE))
  ifelse(shared, d$feature_id, sym)
}

integration_plot_card <- function(output_id, title, sub, download = NULL) {
  bslib::card(
    bslib::card_header(
      htmltools::tags$h3(class = "card-title", title),
      htmltools::tags$span(class = "card-sub", sub),
      download
    ),
    bslib::card_body(hover_plot_output(output_id, height = "320px"))
  )
}

integration_table_card <- function(output_id, title, sub) {
  bslib::card(
    bslib::card_header(
      htmltools::tags$h3(class = "card-title", title),
      htmltools::tags$span(class = "card-sub", sub)
    ),
    bslib::card_body(DT::DTOutput(output_id))
  )
}

integration_ap_card <- function(ns) {
  bslib::card(
    bslib::card_header(
      htmltools::tags$h3(class = "card-title",
                         "ActivePathways \u00B7 combined p"),
      htmltools::tags$span(class = "card-sub",
                           "Brown's method \u00B7 demo data")
    ),
    bslib::card_body(
      DT::DTOutput(ns("ap_table"))
    )
  )
}
