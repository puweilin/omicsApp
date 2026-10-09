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

  # The pathways in the table, in its order: a selected row's index is
  # read against this.
  hits_rows <- shiny::reactive({
    df <- table_data()
    shiny::req(nrow(df) > 0L)
    enrich_hits_rows(df, show_p(), show_cutoff())
  })

  output$hits <- DT::renderDT({
    out <- enrich_hits_table(hits_rows(), show_p(), show_cutoff())
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

  # ---- the selected pathway -------------------------------------------
  # A row of the table opens the pathway under it: for GSEA its
  # running-score curve, for ORA -- which counts genes and has no curve
  # -- the pathway's genes that were in the list, with their change.
  selected_pathway <- shiny::reactive({
    i <- input$hits_rows_selected
    if (!length(i)) return(NULL)
    rows <- tryCatch(hits_rows(), error = function(e) NULL)
    if (is.null(rows) || i[[1L]] > nrow(rows)) return(NULL)
    rows[i[[1L]], , drop = FALSE]
  })

  # The differential result the enrichment was run on, when it is still
  # the one on screen; NULL once the Differential view has moved to
  # another comparison, whose genes and values would not be this
  # pathway's.
  matching_diff <- shiny::reactive({
    eb <- enrich_bundle()
    db <- diff_bundle()
    if (is.null(eb) || is.null(db)) return(NULL)
    if (!identical(eb$params$comparison, db$params$comparison)) return(NULL)
    db
  })

  # Built here, drawn by renderPlot below: kept as its own reactive so
  # the figure can be reused (a download, say) without drawing it twice.
  gsea_curve <- shiny::reactive({
    row <- selected_pathway()
    eb <- enrich_bundle()
    shiny::req(row, eb, identical(eb$params$type, "gsea"), !isTRUE(is_demo()))
    res <- tryCatch(
      omicsCore::plot_gsea(eb, row$pathway_id, database = row$database,
                           diff_bundle = matching_diff()),
      error = function(e) e)
    shiny::validate(shiny::need(
      !inherits(res, "error"),
      paste("The curve could not be drawn:",
            if (inherits(res, "error")) conditionMessage(res) else "")))
    res
  })

  output$gsea_curve <- shiny::renderPlot(
    res = PLOT_RES, alt = "Running enrichment score of the selected pathway",
    fit_to_width("gsea_curve", gsea_curve()))

  # A card of its own under the results row, only once a pathway is
  # picked. Inside the table's card it sat below a scrolling table, in a
  # card as tall as the dot plot beside it, where nobody would find it.
  output$pathway_detail <- shiny::renderUI({
    row <- selected_pathway()
    if (is.null(row)) return(NULL)
    gsea <- identical(plot_bundle()$params$type, "gsea")
    body <- if (isTRUE(is_demo())) {
      notice("Demo pathways only", kind = "info",
             detail = paste("These example pathways come without genes or a",
                            "ranking. Open a project and run the enrichment to",
                            "see a pathway's curve or genes."))
    } else if (gsea) {
      htmltools::tags$div(
        class = "pathway-detail",
        style = "display:flex;flex-wrap:wrap;gap:8px 24px;align-items:flex-start",
        htmltools::tags$div(
          style = "flex:1 1 420px;max-width:760px;min-width:0",
          shiny::plotOutput(session$ns("gsea_curve"), height = "320px")),
        htmltools::tags$p(
          class = "muted", style = "flex:1 1 220px;font-size:12.5px;margin:0",
          htmltools::tags$strong("How to read it. "),
          paste("Genes are ranked from most up-regulated (left) to most",
                "down-regulated (right). Each tick is a gene of this pathway; the",
                "curve climbs at every tick and falls between them, and its peak is",
                "the enrichment score (ES). A peak near the left means the",
                "pathway's genes are mostly up, a dip near the right mostly down.")))
    } else {
      ora_pathway_genes_ui(row, matching_diff(), diff_bundle())
    }
    bslib::card(
      bslib::card_header(
        htmltools::tags$h3(class = "card-title", "Selected pathway"),
        htmltools::tags$span(class = "card-sub",
                             if (gsea) "running enrichment score"
                             else "its genes in the list")),
      bslib::card_body(body))
  })

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
       show_p = show_p, show_cutoff = show_cutoff,
       selected_pathway = selected_pathway, gsea_curve = gsea_curve)
}

# The pathways of the Enriched sets table: those that pass the display
# threshold, strongest first.
enrich_hits_rows <- function(df, show_p, show_cutoff) {
  df <- omicsCore::filter_enrich_results(
    df, p_cutoff = show_cutoff, p_preference = show_p)
  pcol <- if (identical(show_p, "raw")) "p_value" else "adj_p_value"
  df[order(df[[pcol]]), , drop = FALSE]
}

# ORA's answer for one pathway: its genes that were in the tested list,
# with their change in the comparison. There is no running score to
# draw -- ORA counts genes, it does not walk a ranking -- and saying so
# is better than leaving the space empty.
ora_pathway_genes_ui <- function(row, diff_now, diff_any) {
  genes <- strsplit(row$overlap_features %||% "", "/", fixed = TRUE)[[1L]]
  genes <- unique(genes[nzchar(genes) & !is.na(genes)])
  list_name <- switch(row$direction %||% "",
                      up = "the up-regulated genes",
                      down = "the down-regulated genes",
                      "the gene list")
  head_line <- htmltools::tags$div(
    htmltools::tags$strong(row$pathway_name),
    htmltools::tags$div(
      class = "muted", style = "font-size:12.5px",
      sprintf("%s of the pathway's %s genes were among %s \u00B7 adjusted p %s",
              format(length(genes)), format(row$gene_set_size), list_name,
              format_p_plain(row$adj_p_value))))
  why <- htmltools::tags$p(
    class = "muted", style = "font-size:12.5px;margin:6px 0",
    paste("ORA counts how many of a pathway's genes are in the list, so there is no",
          "running score to draw. These are those genes."))
  if (!length(genes)) {
    return(htmltools::tags$div(head_line, notice(
      "No gene names were kept for this pathway", kind = "info",
      detail = "The result does not list which genes overlapped.")))
  }
  tab <- ora_gene_values(genes, diff_now, row$direction)
  if (is.null(tab)) {
    return(htmltools::tags$div(
      style = "max-width:640px", head_line, why,
      htmltools::tags$p(style = "font-size:12.5px", paste(genes, collapse = ", ")),
      htmltools::tags$p(
        class = "muted", style = "font-size:12px",
        if (is.null(diff_any))
          "Their changes are not shown: the differential result is not loaded."
        else paste("Their changes are not shown: the Differential view now shows",
                   "another comparison than the one this enrichment was run on."))))
  }
  eff <- omicsCore::effect_label(diff_now)
  rows <- lapply(seq_len(nrow(tab)), function(i) {
    v <- tab$effect[i]
    colour <- if (!is.finite(v)) "inherit" else if (v > 0) omics_colors$up else omics_colors$down
    htmltools::tags$tr(
      htmltools::tags$td(tab$gene[i]),
      htmltools::tags$td(style = sprintf("text-align:right;color:%s", colour),
                         if (is.finite(v)) sprintf("%+.2f", v) else "\u2013"),
      htmltools::tags$td(style = "text-align:right", format_p_plain(tab$adj_p[i])))
  })
  htmltools::tags$div(
    class = "pathway-detail", style = "max-width:640px",
    head_line, why,
    htmltools::tags$div(
      style = "max-height:260px;overflow-y:auto",
      htmltools::tags$table(
        class = "table table-sm ora-genes", style = "font-size:12.5px;margin:0",
        htmltools::tags$thead(htmltools::tags$tr(
          htmltools::tags$th("Gene"),
          htmltools::tags$th(style = "text-align:right", eff),
          htmltools::tags$th(style = "text-align:right", "adjusted p"))),
        htmltools::tags$tbody(rows))))
}

# The overlap genes' values in the differential result, most changed
# first in the direction of the list they were found in. A symbol
# measured more than once (protein isoforms) takes its most significant
# row. NULL when there is no matching result to read from.
ora_gene_values <- function(genes, diff_bundle, direction = NA_character_) {
  if (is.null(diff_bundle)) return(NULL)
  df <- diff_bundle$results$diff_result_df
  if (is.null(df) || !nrow(df)) return(NULL)
  df <- df[order(df$adj_p_value %||% df$p_value), , drop = FALSE]
  key <- if ("feature_symbol" %in% names(df)) df$feature_symbol else df$feature_id
  i <- match(genes, key)
  # The run may have matched the names ignoring case.
  miss <- is.na(i)
  if (any(miss)) i[miss] <- match(toupper(genes[miss]), toupper(key))
  out <- data.frame(gene = genes,
                    effect = df$effect[i],
                    adj_p = (df$adj_p_value %||% df$p_value)[i],
                    stringsAsFactors = FALSE)
  ord <- switch(direction %||% "",
                up = order(-out$effect),
                down = order(out$effect),
                order(-abs(out$effect)))
  out[ord, , drop = FALSE]
}

format_p_plain <- function(p) {
  if (length(p) != 1L || is.na(p)) return("\u2013")
  if (p < 0.001) formatC(p, digits = 1, format = "e") else formatC(p, digits = 3, format = "fg")
}

# The rows of the Enriched sets table, as shown.
enrich_hits_table <- function(df, show_p, show_cutoff) {
  df <- enrich_hits_rows(df, show_p, show_cutoff)
  pcol <- if (identical(show_p, "raw")) "p_value" else "adj_p_value"
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
                           "ranked by adjusted p \u00B7 click a pathway for more"),
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
