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

  # Built here, drawn below: the card's download saves this same figure.
  dot_plot <- shiny::reactive({
    b <- plot_bundle()
    shiny::req(b)
    omicsCore::plot_enrichment(b, view = "dot", top_n = 12L,
                               p_preference = show_p(),
                               p_cutoff = show_cutoff())
  })
  output$dot <- shiny::renderPlot(
    res = PLOT_RES, alt = "Dot plot of the most enriched pathways",
    height = function() {
      w <- session$clientData[[paste0("output_", session$ns("dot"), "_width")]]
      enrich_dot_height(tryCatch(dot_plot(), error = function(e) NULL),
                        narrow = is.numeric(w) && length(w) && w < NARROW_PLOT_PX)
    },
    fit_to_width("dot", dot_plot()))
  # The axis names are wrapped and long ones hard to read; hovering (or
  # tapping) a dot gives the pathway in full with its numbers.
  plot_hover_server("dot", dot_plot,
                    function(row, p) enrich_dot_hover_text(row, p, show_p()),
                    input, output, session)

  # The pathways in the table, in its order: a selected row's index is
  # read against this.
  hits_rows <- shiny::reactive({
    df <- table_data()
    shiny::req(nrow(df) > 0L)
    enrich_hits_rows(df, show_p(), show_cutoff())
  })

  output$hits <- DT::renderDT({
    tab <- enrich_hits_table(hits_rows(), show_p(), show_cutoff())
    shiny::req(nrow(tab$data) > 0L)
    enrich_hits_datatable(tab)
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
    # Its title is the pathway's name, which the card header does not say.
    fit_to_width("gsea_curve", gsea_curve(), keep_title = TRUE))

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
                             else "its genes in the list"),
        # Shown by its server only while there is a curve: ORA's genes
        # are a table, not a figure.
        plot_download_ui(session$ns("gsea_curve_download"))),
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
       selected_pathway = selected_pathway, gsea_curve = gsea_curve,
       dot_plot = dot_plot)
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
  if (length(p) != 1L) return("\u2013")
  format_p_value(p)
}

# The rows of the Enriched sets table, as shown: list(data =, column_defs =)
# from dt_sortable_text(), the p-values written out and sorted by number.
#
# The card is about 390 px wide on a 1280 px screen, beside the dot plot,
# and five columns did not fit in it. What a method does not need is
# left out: ORA has no enrichment score, GSEA's direction is the sign of
# its NES, and GSEA counts no overlap -- its pathway size is shown
# instead. The pathway name is the one column that can give up width;
# it wraps.
enrich_hits_table <- function(df, show_p, show_cutoff) {
  df <- enrich_hits_rows(df, show_p, show_cutoff)
  pcol <- if (identical(show_p, "raw")) "p_value" else "adj_p_value"
  # Named for the column it holds, so the table cannot say adjusted p
  # over raw values.
  p_name <- if (identical(show_p, "raw")) "p" else "adjusted p"
  na_col <- function(x) if (is.null(x)) rep(NA_real_, nrow(df)) else x
  gsea <- any(is.finite(df$effect))
  out <- data.frame(Pathway = df$pathway_name, check.names = FALSE,
                    stringsAsFactors = FALSE)
  keys <- list()
  if (gsea) {
    out$NES <- ifelse(is.finite(df$effect), sprintf("%+.2f", df$effect), "\u2013")
    keys$NES <- df$effect
  } else if (any(df$direction %in% c("up", "down"))) {
    # ORA's direction is the gene list a pathway was found in, and the
    # header says so: "up" under Direction read as the pathway's
    # activity going up, which ORA does not measure. A result pooled
    # from both lists has no such column.
    out$`Gene list` <- ifelse(df$direction %in% c("up", "down"), df$direction, "\u2013")
  }
  out[[p_name]] <- df[[pcol]]
  if (gsea) {
    out$`Set size` <- as.integer(na_col(df$gene_set_size))
  } else {
    ov <- na_col(df$overlap_size)
    size <- na_col(df$gene_set_size)
    out$Overlap <- ifelse(is.na(ov), "\u2013",
                          ifelse(is.na(size), format(ov),
                                 sprintf("%d/%d", as.integer(ov), as.integer(size))))
    keys$Overlap <- ov
  }
  dt_sortable_text(out, p_cols = p_name, sort_keys = keys)
}

# The Enriched sets widget for enrich_hits_table()'s rows.
#
# Numbers on the right and on one line; the pathway takes the rest of
# the width and wraps (enrich_hits_css() tightens the cells). Column
# widths are left to the browser: DataTables' own sizing set the table
# a couple of pixels wider than its card, which then scrolled sideways
# for nothing. On a phone the columns cannot all fit beside a readable
# name, so the table scrolls sideways inside its card with the pathway
# column held in place: a row of numbers without its name means nothing.
enrich_hits_datatable <- function(tab) {
  out <- tab$data
  shown <- names(out)[!startsWith(names(out), ".sort_")]
  col <- function(name) match(name, names(out)) - 1L
  numeric_cols <- intersect(c("NES", "adjusted p", "p", "Overlap", "Set size"), shown)
  defs <- c(list(list(className = "dt-right", targets = col(numeric_cols))),
            tab$column_defs)
  w <- DT::datatable(
    out,
    rownames   = FALSE,
    selection  = "single",
    extensions = "FixedColumns",
    options    = list(
      pageLength   = 10,
      dom          = "ftip",
      autoWidth    = FALSE,
      scrollX      = TRUE,
      fixedColumns = list(left = 1),
      columnDefs   = defs
    )
  )
  # A name with no space to break at (a long gene set id) would
  # otherwise set the column's width by itself.
  w <- DT::formatStyle(w, "Pathway", `overflow-wrap` = "break-word")
  if ("Gene list" %in% shown) {
    w <- DT::formatStyle(w, "Gene list", color = DT::styleEqual(
      c("up", "down"), c(omics_colors$up, omics_colors$down)))
  }
  w
}

# The dot plot's height for its rows, in pixels.
#
# At a fixed 420 px, two databases with up and down lists made three
# panels of long GO names, wrapped to two or three lines each, and the
# names were drawn over one another. The panels of a facet_wrap() are
# equally tall, so the tallest panel's need sets them all. Rows are
# evenly spaced, so a panel needs as many rows as it has times the
# space of its most crowded neighbours -- two adjacent three-line names
# need three lines each -- at the axis text's line height (smaller on a
# phone, where the legends also move under the plot). A result of
# short names in one panel keeps the old 420 px.
ENRICH_DOT_MIN_PX <- 420
enrich_dot_height <- function(p, narrow = FALSE) {
  df <- if (inherits(p, "ggplot")) p$data else NULL
  if (!is.data.frame(df) || !nrow(df) || !".label" %in% names(df)) return(ENRICH_DOT_MIN_PX)
  list_col <- if (".list" %in% names(df)) as.character(df$.list) else rep("", nrow(df))
  panel <- paste(df$database, list_col)
  # In the order the axis draws them.
  ord <- order(as.integer(factor(df$.label)))
  lines <- lengths(strsplit(as.character(df$.label), "\n", fixed = TRUE))[ord]
  need <- tapply(lines, panel[ord], function(l) {
    pair <- if (length(l) > 1L) max((l[-1L] + l[-length(l)]) / 2) else l
    # Rows of one line still need a little air between them.
    length(l) * max(1.25, pair)
  })
  line_px <- if (isTRUE(narrow)) 13 else 18
  h <- 110 + 40 * length(need) + max(need) * line_px * length(need) +
    # Room for the keys, which sit under the panel on every width.
    if (isTRUE(narrow)) 170 else 130
  round(max(ENRICH_DOT_MIN_PX, min(1800, h)))
}

enrich_dot_card <- function(ns) {
  bslib::card(
    bslib::card_header(
      htmltools::tags$h3(class = "card-title", "Pathway dotplot"),
      htmltools::tags$span(class = "card-sub",
                           # Size is genes in the list for ORA and set size for
                           # GSEA, and the p can be raw: the keys say which.
                           "top 12 per database \u00B7 hover a dot for details"),
      plot_download_ui(ns("dot_download"))
    ),
    bslib::card_body(
      hover_plot_output(ns("dot"), height = "auto")
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
      htmltools::tags$style(htmltools::HTML(enrich_hits_css(ns("hits")))),
      DT::DTOutput(ns("hits"))
    )
  )
}

# Tighter cells than DataTables' own, for a table that shares a row with
# the dot plot: its header kept 26 px beside every title for the sort
# arrows, which on four columns was a fifth of the card. Numbers stay on
# one line; the held pathway column takes the page colour, not white,
# so it does not stand out from the rows when nothing is scrolled.
enrich_hits_css <- function(id) {
  sel <- paste0("#", id, " table.dataTable")
  paste0(
    sel, "{font-size:13px}",
    # DataTables' own rules name the sort state and the table classes,
    # and outrank a rule this short.
    sel, ">thead>tr>th{padding:6px 14px 6px 5px!important;vertical-align:bottom}",
    sel, ">thead>tr>th:before,", sel, ">thead>tr>th:after{right:3px!important}",
    sel, ">tbody>tr>td{padding:6px 5px!important}",
    sel, ">tbody>tr>td.dt-right{white-space:nowrap}",
    sel, " .dtfc-fixed-left{background-color:var(--bg)}")
}
