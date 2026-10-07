# What the enrichment view shows: header, the gene-list summary, notices,
# the dot plot, the table and its download.
#
# Split out of mod_enrich_view.R and called from inside its module
# server, so the outputs keep their ids ("enrich-dot", "enrich-hits").
# Returns the reactives the module keeps.
enrich_results_server <- function(input, output, session, navigate, diff_bundle, diff_layer,
                                  diff_thresholds, current_project, have_cp, enrich_bundle,
                                  enrich_error, is_demo, settings_changed) {
  # The table powering both the dot card and the hits card. In
  # demo mode this is the static fixture; in live mode it's
  # the enrichment bundle's standardized result data frame. We
  # guard the live path with `req()` so a transient NULL during
  # an async flush doesn't crash the reactive graph (downstream
  # renderers also `req()` on `nrow(df) > 0`).
  table_data <- shiny::reactive({
    if (isTRUE(is_demo())) return(example_enrich_table())
    b <- enrich_bundle()
    shiny::req(b)
    b$results$enrich_result_df
  })

  output$header <- shiny::renderUI({
    b <- enrich_bundle()
    demo <- is_demo()
    omics <- if (is.null(b)) "Proteomics"
             else omics_display(b$input_info$omics_type)
    method <- if (is.null(b)) "MSigDB Hallmark"
              else sprintf("%s \u00B7 %s",
                           toupper(b$params$type %||% "ora"),
                           b$params$database %||% "hallmark")
    # Which layer, stated rather than offered. Enrichment runs on a
    # differential result, so its layer is whichever one Differential
    # was run on -- a picker here could be set to disagree with that,
    # and the panel would then be labelled rnaseq while showing
    # proteomics pathways. To change it, change it there.
    layer <- diff_layer_tag()
    view_header(
      title    = "Pathway enrichment",
      actions  = if (is.function(navigate) && !demo) {
        shiny::actionButton(session$ns("go_next"), "Next: Integration \u2192",
                            class = "btn btn-ghost")
      },
      subtitle = htmltools::tagList(
        omics,
        htmltools::HTML(" &middot; "),
        method,
        if (!is.null(layer)) {
          htmltools::tagList(
            htmltools::HTML(" &middot; "),
            htmltools::tags$span(
              class = "muted",
              sprintf("from the %s differential", layer))
          )
        },
        htmltools::HTML(" &middot; "),
        htmltools::tags$span(
          class = "muted",
          if (demo) "demo data (built-in)"
          else "live result"
        )
      )
    )
  })

  # The tag of the layer the upstream differential ran on, handed over
  # by that view. NULL for the demo.
  diff_layer_tag <- shiny::reactive({
    if (isTRUE(is_demo())) return(NULL)
    tag <- diff_layer()
    if (is.null(tag) || !nzchar(tag)) NULL else tag
  })

  # The size of the gene list this enrichment is over.
  #
  # An empty result reads identically whether no feature met the
  # threshold or three thousand did and none of the sets were
  # enriched. Reported as "100+ differential genes but every database
  # comes back with nothing", which is a question this answers before
  # anyone has to ask it.
  selected_features <- shiny::reactive({
    b <- diff_bundle()
    if (is.null(b)) return(NULL)
    df <- b$results$diff_result_df
    if (is.null(df)) return(NULL)
    thr <- diff_thresholds()
    # Once an ORA has run, the thresholds it ran at: the live ones may
    # have moved since, and the line describes the result below it.
    eb <- enrich_bundle()
    if (!is.null(eb) && identical(eb$params$type, "ora") &&
        identical(input$type %||% "ora", "ora")) {
      thr <- list(p_cutoff = eb$params$p_cutoff,
                  p_preference = eb$params$p_preference,
                  effect_cutoff = eb$params$effect_cutoff)
    }
    # GSEA ranks every feature it can place, so there is no selection
    # to report -- only how many made it into the ranking. Saying "N
    # of M selected" here would describe a step the method does not
    # take.
    if (!identical(input$type %||% "ora", "ora")) {
      rankable <- !is.na(df$effect) & !is.na(df$feature_symbol)
      return(list(gsea = TRUE, n = sum(rankable), total = nrow(df),
                  mapped = sum(rankable & is_gene_symbol(df$feature_symbol))))
    }
    pcol <- if (identical(thr$p_preference, "raw")) "p_value" else "adj_p_value"
    pv <- df[[pcol]]
    keep <- !is.na(pv) & pv < (thr$p_cutoff %||% 0.05)
    if (!is.null(thr$effect_cutoff) && is.finite(thr$effect_cutoff)) {
      keep <- keep & !is.na(df$effect) & abs(df$effect) >= thr$effect_cutoff
    }
    list(gsea = FALSE, n = sum(keep), total = nrow(df), p_col = pcol,
         p_cutoff = thr$p_cutoff %||% 0.05,
         effect_cutoff = thr$effect_cutoff,
         # A symbol is a symbol whether or not it is also the id: a
         # matrix keyed by gene symbols was told "none carry a gene
         # symbol" beside the pathways it had just found.
         mapped = sum(keep & is_gene_symbol(df$feature_symbol)))
  })

  output$input_summary <- shiny::renderUI({
    s <- selected_features()
    if (is.null(s)) return(NULL)
    if (isTRUE(s$gsea)) {
      return(htmltools::tags$div(
        class = "muted", style = "font-size:12.5px;margin:2px 0 10px",
        htmltools::tags$strong(sprintf("%d of %d features ranked",
                                       s$n, s$total)),
        " \u2014 GSEA reads the whole list, so the Differential ",
        "thresholds do not apply here.",
        if (s$mapped == 0L) {
          htmltools::tags$span(
            style = "color:var(--warn)",
            " None carry a gene symbol, so no set will match."
          )
        }
      ))
    }
    p_name <- switch(s$p_col, adj_p_value = "adjusted p", p_value = "p", s$p_col)
    bits <- sprintf("%s < %s", p_name, format(s$p_cutoff))
    if (!is.null(s$effect_cutoff) && is.finite(s$effect_cutoff)) {
      bits <- paste0(bits, sprintf(", |%s| \u2265 %.3f",
                                   omicsCore::effect_label(diff_bundle()), s$effect_cutoff))
    }
    htmltools::tags$div(
      class = "muted", style = "font-size:12.5px;margin:2px 0 10px",
      htmltools::tags$strong(sprintf("%d of %d features", s$n, s$total)),
      sprintf(" selected at %s", bits),
      if (s$n == 0L) {
        htmltools::tags$span(
          style = "color:var(--warn)",
          " \u2014 nothing to enrich. Loosen the thresholds in the ",
          "Differential view."
        )
      } else if (s$mapped == 0L) {
        htmltools::tags$span(
          style = "color:var(--warn)",
          " \u2014 none carry a gene symbol, so no set will match. ",
          "Check the Gene symbol line on the Import page."
        )
      }
    )
  })

  output$notices <- shiny::renderUI({
    tagged <- htmltools::tagList()
    err <- enrich_error()
    if (is.null(err) && isTRUE(settings_changed())) {
      tagged <- htmltools::tagAppendChild(tagged, notice(
        "The settings have changed since this result was computed",
        paste("The pathways below are from the previous settings (test, database,",
              "species, gene lists or the Differential thresholds). Press Re-run to update them."),
        kind = "warn"))
    }
    if (!is.null(err)) {
      tagged <- htmltools::tagAppendChild(
        tagged,
        if (identical(err, CANCELLED_MESSAGE)) {
          notice("Cancelled", "Press Re-run to start again.", kind = "info")
        } else if (!have_cp) {
          notice(title = "clusterProfiler unavailable", detail = err,
                 kind = "warn")
        } else {
          notice(title  = "The enrichment could not be computed",
                 detail = if (grepl("No significant|no features|no genes", err, ignore.case = TRUE))
                   "No differential hits pass the current thresholds; loosen them in the Differential view or try GSEA."
                 else "See the technical details below.",
                 kind   = "error",
                 technical = err)
        }
      )
    }
    if (is.null(diff_bundle()) && have_cp) {
      tagged <- htmltools::tagAppendChild(
        tagged,
        notice(
          title  = "Run a differential analysis first",
          detail = if (is.null(current_project()))
            paste0("This view enriches the top hits from the ",
                   "Differential view. Showing the demo fixture ",
                   "until a real result arrives.")
          else "This view enriches the hits of the Differential view; run it there and the pathways appear here.",
          kind   = "info"
        )
      )
    }
    # Nothing found: say the likely reason instead of an empty table.
    eb <- enrich_bundle()
    if (!is.null(eb) && !isTRUE(is_demo()) &&
        nrow(eb$results$enrich_result_df %||% data.frame()) == 0L) {
      guess <- guess_organism(diff_bundle()$results$diff_result_df$feature_symbol)
      used <- species_code(eb$params$organism) %||% "Hs"
      tagged <- htmltools::tagAppendChild(tagged, notice(
        title = "No pathway was found",
        detail = if (!identical(guess, used))
          sprintf("The gene names look like %s genes (e.g. %s), but %s gene sets were used. Switch Species and re-run.",
                  tolower(species_label(guess)), species_example(guess),
                  tolower(species_label(used)))
        else "None of the sets passed the threshold. Try GSEA, another database, or looser Differential thresholds.",
        kind = "warn"))
    }
    # What the run itself had to say -- the gene names matched to the
    # gene sets only once upper and lower case were ignored, say.
    if (!is.null(eb) && !isTRUE(is_demo()) && length(eb$warnings)) {
      for (w in eb$warnings) {
        tagged <- htmltools::tagAppendChild(tagged, notice(
          "Note on this result", detail = w, kind = "info"))
      }
    }
    tagged
  })

  # Both paths hand a bundle to the same dispatcher. When the demo
  # had its own drawing code it was free to drift from the live view,
  # and a demo that no longer resembles the product is worse than no
  # demo.
  plot_bundle <- shiny::reactive({
    if (isTRUE(is_demo())) example_enrich_bundle() else enrich_bundle()
  })

  # The threshold the panel is read at. numericInput reports NA while
  # the box is mid-edit, and NA would empty the panel with nothing to
  # say why.
  show_p <- shiny::reactive({
    v <- input$show_p
    if (is.null(v) || !is.character(v) || length(v) != 1L ||
        !v %in% c("adjusted", "raw", "qvalue")) "adjusted" else v
  })
  show_cutoff <- shiny::reactive({
    v <- input$show_cutoff
    if (is.null(v) || !is.finite(v) || v <= 0 || v > 1) 0.05 else v
  })

  output$dot <- shiny::renderPlot(res = PLOT_RES, alt = "Dot plot of the most enriched pathways", fit_to_width("dot", {
    b <- plot_bundle()
    shiny::req(b)
    omicsCore::plot_enrichment(b, view = "dot", top_n = 12L,
                               p_preference = show_p(),
                               p_cutoff = show_cutoff())
  }))

  output$hits <- DT::renderDT({
    df <- table_data()
    shiny::req(nrow(df) > 0L)
    out <- enrich_hits_table(df, show_p(), show_cutoff())
    shiny::req(nrow(out) > 0L)
    DT::datatable(
      out,
      rownames  = FALSE,
      selection = "single",
      options   = list(
        pageLength = 10,
        dom        = "ftip",
        scrollX    = TRUE,
        columnDefs = list(list(className = "dt-right",
                               targets = c(1, 2)))
      )
    )
  }, server = TRUE)

  output$download_table <- shiny::downloadHandler(
    filename = function() {
      b <- if (isTRUE(is_demo())) NULL else enrich_bundle()
      type <- toupper(b$params$type %||% "enrichment")
      sprintf("%s_%s.csv", tolower(type), format(Sys.Date(), "%Y%m%d"))
    },
    content = function(file) {
      df <- table_data()
      shiny::req(nrow(df) > 0L)

      # Unfiltered and every column, including the gene lists -- those
      # are what a follow-up actually needs, and they are the one thing
      # the on-screen table cannot show. Sorted so the file opens on
      # the strongest result rather than in database order.
      df <- df[order(df$database,
                     df$adj_p_value %||% df$p_value,
                     df$p_value), , drop = FALSE]
      utils::write.csv(df, file, row.names = FALSE, na = "")
    }
  )

  list(table_data = table_data, diff_layer_tag = diff_layer_tag,
       selected_features = selected_features, plot_bundle = plot_bundle,
       show_p = show_p, show_cutoff = show_cutoff)
}

# The rows of the Enriched sets table: the pathways that pass the
# display threshold, strongest first.
enrich_hits_table <- function(df, show_p, show_cutoff) {
  df <- omicsCore::filter_enrich_results(
    df, p_cutoff = show_cutoff, p_preference = show_p)
  pcol <- if (identical(show_p, "raw")) "p_value" else "adj_p_value"
  df <- df[order(df[[pcol]]), , drop = FALSE]
  out <- data.frame(
    Pathway   = df$pathway_name,
    NES       = sprintf("%+.2f", df$effect),
    # Named for the column it holds, so the table cannot say adj.P
    # over raw values.
    P         = signif(df[[pcol]], 3),
    Direction = df$direction,
    Overlap   = sprintf("%d/%d",
                        df$overlap_size %||% NA_integer_,
                        df$gene_set_size %||% NA_integer_),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  names(out)[names(out) == "P"] <-
    if (identical(show_p, "raw")) "p" else "adjusted p"
  if (!nrow(df)) return(out)
  # ORA's direction is the gene list a pathway was found in, which is
  # worth saying in words: "up" next to an ORA pathway read as the
  # pathway's activity going up, which ORA does not measure.
  if (all(is.na(df$effect))) {
    out$Direction <- ifelse(df$direction %in% "up", "up-regulated genes",
                     ifelse(df$direction %in% "down", "down-regulated genes",
                            NA_character_))
    names(out)[names(out) == "Direction"] <- "Found among"
  }
  # ORA has no enrichment score, and a direction only when it was run on
  # one direction or on both separately; columns of NA said nothing.
  if (all(is.na(df$effect))) out$NES <- NULL
  if (all(is.na(df$direction))) {
    out <- out[setdiff(names(out), c("Direction", "Found among"))]
  }
  out
}

enrich_dot_card <- function(ns) {
  bslib::card(
    bslib::card_header(
      htmltools::tags$h3(class = "card-title", "Pathway dotplot"),
      htmltools::tags$span(class = "card-sub",
                           "top 12 by adjusted p \u00B7 size = overlap")
    ),
    bslib::card_body(
      shiny::plotOutput(ns("dot"), height = "420px")
    )
  )
}

enrich_hits_card <- function(ns) {
  bslib::card(
    bslib::card_header(
      htmltools::tags$h3(class = "card-title", "Enriched sets"),
      htmltools::tags$span(class = "card-sub",
                           "ranked by adjusted p"),
      # Everything, not what the table happens to be showing: the CSV is
      # what gets opened in Excel a week later, and a file silently
      # truncated to one significance threshold is the kind of thing
      # nobody notices until a number cannot be reproduced.
      shiny::downloadButton(
        ns("download_table"), "Download full table",
        class = "btn btn-ghost btn-sm"
      )
    ),
    bslib::card_body(
      DT::DTOutput(ns("hits"))
    )
  )
}
