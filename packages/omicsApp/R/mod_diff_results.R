# The results half of the differential view: header, notices, stat
# cards, volcano, hit table and the full-table downloads.
#
# Split out of mod_diff_view.R. These are plain functions called from
# inside diff_view_server()'s moduleServer(), not modules, so every
# output keeps its id ("diff-volcano", "diff-download_csv").

diff_results_server <- function(input, output, session, navigate, active, shown_bundle,
                                diff_bundle, diff_error, comparisons, settings_changed,
                                marked, p_col, p_label, fdr_cut_d, fc_cut_d) {
  output$header <- shiny::renderUI({
    a <- active()
    b <- shown_bundle()
    source_note <- htmltools::tags$span(
      class = "muted",
      if (a$is_demo) "demo project (built-in)"
      else sprintf("layer = %s", a$tag)
    )

    # Before the first result there is nothing to summarise. Three em
    # dashes separated by middots said that in a way that read as
    # damage rather than as an empty state, which is how it was
    # reported.
    if (is.null(b)) {
      return(view_header(
        title    = "Differential",
        subtitle = htmltools::tagList(
          htmltools::tags$span(class = "muted", "no result yet"),
          htmltools::HTML(" &middot; "),
          source_note
        )
      ))
    }

    stale <- !is.null(diff_error()) && !identical(diff_error(), CANCELLED_MESSAGE)
    comparison <- if (length(b$params$comparison)) gsub("_vs_", " vs ", b$params$comparison) else
      sprintf("%s vs %s",
              b$params$case_group %||% input$case %||% "case",
              b$params$control_group %||% input$control %||% "control")
    n_cmp <- length(comparisons())
    view_header(
      title    = "Differential",
      actions  = if (is.function(navigate)) {
        shiny::actionButton(session$ns("go_enrich"), "Next: Enrichment \u2192",
                            class = "btn btn-ghost")
      },
      subtitle = htmltools::tagList(
        omics_display(b$input_info$omics_type),
        if (n_cmp > 1L) htmltools::tagList(
          htmltools::HTML(" &middot; "),
          sprintf("%d comparisons", n_cmp)),
        htmltools::HTML(" &middot; "),
        comparison,
        htmltools::HTML(" &middot; "),
        b$params$method,
        htmltools::HTML(" &middot; "),
        source_note,
        if (stale) htmltools::tagList(
          htmltools::HTML(" &middot; "),
          htmltools::tags$span(style = "color:var(--warn);font-weight:600",
                               "previous result \u2014 the latest run failed"))
      )
    )
  })

  output$notices <- shiny::renderUI({
    err <- diff_error()
    missing_engines <- diff_missing_engines()
    tagged <- htmltools::tagList()
    if (identical(err, CANCELLED_MESSAGE)) {
      tagged <- htmltools::tagAppendChild(tagged, notice(
        "Cancelled", "Press Run analysis to start again.", kind = "info"))
    } else if (!is.null(err)) {
      # The engine's own sentence when it is one: hiding "'control'
      # is not a group. Groups: Control, TreatA" behind "See the
      # technical details below" made users open a fold to read it.
      plain <- !looks_internal_error(err)
      tagged <- htmltools::tagAppendChild(
        tagged,
        notice(title  = if (is.null(diff_bundle())) "The differential analysis could not run"
                        else "The latest run failed; the result below is from the previous run",
               detail = if (plain) err else diff_error_hint(err),
               kind   = "error",
               technical = if (!plain) err)
      )
    }
    if (isTRUE(settings_changed()) && is.null(err)) {
      tagged <- htmltools::tagAppendChild(tagged, notice(
        "The settings have changed since this result was computed",
        "The plots and tables below still show the previous result. Press Re-run to update them.",
        kind = "warn"))
    }
    # What the engine said about the result: an ignored covariate, a
    # scale conversion, genes too low to test.
    warns <- diff_bundle()$warnings
    if (length(warns)) {
      tagged <- htmltools::tagAppendChild(tagged, notice(
        title = "Notes on this result",
        detail = htmltools::tags$ul(lapply(unique(warns), htmltools::tags$li)),
        kind = "info"))
    }
    if (is.null(diff_bundle()) && is.null(err)) {
      tagged <- htmltools::tagAppendChild(
        tagged,
        notice(title  = "No result yet",
               detail = paste("Check the layer, the control group and the groups to",
                              "compare in the Parameters card, then press Run analysis."),
               kind   = "info")
      )
    }
    if (length(missing_engines)) {
      tagged <- htmltools::tagAppendChild(
        tagged,
        notice(
          title  = "Some engines are unavailable",
          detail = sprintf(
            "Not installed: %s. Install with `omicsCore::install_optional()`.",
            paste(missing_engines, collapse = ", ")
          ),
          kind = "info"
        )
      )
    }
    tagged
  })

  output$stats <- shiny::renderUI({
    b <- shown_bundle()
    if (is.null(b)) return(NULL)
    # A weighted contrast has no single case group; its "up" is the
    # direction of the contrast as written.
    case_lbl <- b$params$case_group %||% "contrast"
    if (length(case_lbl) > 1L) case_lbl <- "case"
    df <- marked()
    sig <- df[df$is_significant, , drop = FALSE]
    up_n   <- sum(sig$effect > 0, na.rm = TRUE)
    down_n <- sum(sig$effect < 0, na.rm = TRUE)
    top    <- if (nrow(sig) > 0L) sig[which.max(abs(sig$effect)), ] else NULL
    top_value <- if (is.null(top)) "\u2014" else as.character(top$feature_symbol[1L])
    top_trend <- if (is.null(top)) "no features pass thresholds"
                 else sprintf("%s %+.2f \u00B7 %s %.2g", omicsCore::effect_label(b),
                              top$effect[1L], p_label(), top[[p_col()]][1L])
    htmltools::tags$div(
      class = "stat-grid",
      {
        # Features with no p-value were not tested (too many missing
        # values, or set aside as too low to test); counting them as
        # tested overstated the screen.
        n_tested <- sum(!is.na(df$p_value))
        stat_card(
          label = "Tested features",
          value = if (n_tested < nrow(df))
                    sprintf("%s / %s", format(n_tested, big.mark = ","), format(nrow(df), big.mark = ","))
                  else format(nrow(df), big.mark = ","),
          trend = if (n_tested < nrow(df))
                    sprintf("%s not testable (missing values or too few counts)",
                            format(nrow(df) - n_tested, big.mark = ","))
                  else sprintf("%s, %s", b$params$method, b$params$comparison %||% "\u2014"),
          mono  = TRUE
        )
      },
      stat_card(
        label  = sprintf("Up in %s", case_lbl),
        value  = up_n,
        trend  = sprintf("%s > %.2f \u00B7 %s < %.3g", omicsCore::effect_label(b),
                         fc_cut_d(), p_label(), fdr_cut_d()),
        accent = "up"
      ),
      stat_card(
        label  = sprintf("Down in %s", case_lbl),
        value  = down_n,
        trend  = sprintf("%s < -%.2f \u00B7 %s < %.3g", omicsCore::effect_label(b),
                         fc_cut_d(), p_label(), fdr_cut_d()),
        accent = "down"
      ),
      stat_card(
        label = "Top hit",
        value = top_value,
        trend = top_trend,
        mono  = TRUE
      )
    )
  })

  # plotly drops a ggplot caption, so the card says what "significant"
  # means on the figure: the same cut as the table beside it.
  output$volcano_cut <- shiny::renderText({
    b <- shown_bundle()
    if (is.null(b)) return("")
    volcano_cut_text(p_label(), fdr_cut_d(), omicsCore::effect_label(b), fc_cut_d())
  })

  # The figure's legend, with the counts plot_volcano() puts in its own
  # (the same mask as the stat cards, so the numbers match them too).
  output$volcano_legend <- shiny::renderUI({
    if (is.null(shown_bundle())) return(volcano_legend())
    df <- marked()
    volcano_legend(
      up_n   = sum(df$is_significant & df$effect > 0, na.rm = TRUE),
      down_n = sum(df$is_significant & df$effect < 0, na.rm = TRUE),
      continuous = startsWith(as.character(df$analysis_type[1L] %||% ""), "continuous"))
  })

  output$volcano <- plotly::renderPlotly({
    b <- shown_bundle()
    shiny::validate(shiny::need(b, "Press Run analysis to draw the volcano."))
    # Drawn at the thresholds the controls set, so its colours agree
    # with the hit table and the stat cards beside it. The same values
    # are saved with the project (params$display_thresholds), so the
    # report and the exported script draw this same figure.
    # The labels are added as plotly annotations, not drawn by
    # plot_volcano(): ggplotly() cannot convert ggrepel's text layer and
    # dropped it with a warning, so "Label top hits" labelled nothing.
    p <- omicsCore::plot_volcano(b, top_n = 0L, p_basis = volcano_p_basis(p_col()),
                                 p_threshold = fdr_cut_d(),
                                 effect_threshold = volcano_effect_cut(fc_cut_d()))
    # plotly's own legend is hidden: the card's legend below the plot
    # says the same, with the counts, and does not take width from the
    # plot on a phone.
    # No plotly title: the card already says "Volcano", and at phone
    # width the title ran into the plotly toolbar.
    fig <- plotly::ggplotly(p + ggplotly_untitled(), tooltip = "text") |> drop_hoveron() |>
      plotly::layout(showlegend = FALSE)
    if (isTRUE(input$label_top)) {
      # Laid out for the width the browser reports, so a phone gets
      # fewer labels rather than overlapping ones.
      fig <- label_volcano(fig, b, 20L, p_col(),
                           width = session$clientData[[paste0("output_", session$ns("volcano"),
                                                              "_width")]])
    }
    # WebGL rather than one SVG node per point: 60,000 genes painted in
    # 0.5 s instead of 4 s, with the same points and hover text. Only
    # where the browser has WebGL: without it (remote desktops, some
    # locked-down or GPU-less machines) plotly drew a grey box saying
    # "WebGL is not supported" and no volcano at all.
    if (!identical(session$rootScope()$input$omics_webgl, FALSE)) {
      fig <- plotly::toWebGL(fig)
    }
    plotly::config(fig, displaylogo = FALSE,
                   modeBarButtonsToRemove = c("lasso2d", "select2d"))
  })

  # An empty card before the first run read as a broken one.
  output$hits_empty <- shiny::renderUI({
    if (!is.null(diff_bundle())) return(NULL)
    htmltools::tags$p(class = "muted", style = "font-size:13px;margin:4px 0",
                      "Run the analysis to see the features that pass the thresholds.")
  })

  output$hits <- DT::renderDT({
    df <- marked()
    sig <- df[df$is_significant, , drop = FALSE]
    sig <- sig[order(-abs(sig$effect)), , drop = FALSE]
    out <- data.frame(
      Feature   = sig$feature_symbol,
      Effect    = round(sig$effect, 3),
      p         = signif(sig[[p_col()]], 3),
      Direction = ifelse(sig$effect > 0, "up", "down"),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    # The columns are named for what they hold: the p-value the mask
    # was read from, and the effect in the words the cards use.
    names(out)[2] <- omicsCore::effect_label(shown_bundle())
    names(out)[3] <- p_label()
    DT::datatable(
      out,
      rownames  = FALSE,
      selection = "single",
      options   = list(
        pageLength = 10,
        dom        = "ftip",
        language   = list(emptyTable = "No feature passes the current thresholds."),
        scrollX    = TRUE,
        columnDefs = list(list(className = "dt-right", targets = 1:2))
      )
    )
  }, server = TRUE)
  invisible()
}

diff_downloads_server <- function(input, output, session, active, diff_bundle, p_col,
                                  fdr_cut_d, fc_cut_d) {
  # The whole result, every feature and every comparison, with the
  # significance at the current thresholds: the table a user takes to
  # a paper or a colleague. The hit list on screen is one page of one
  # comparison.
  full_table <- function() {
    b <- diff_bundle()
    shiny::req(b)
    df <- b$results$diff_result_df
    pv <- df[[p_col()]]
    df$significant <- !is.na(pv) & !is.na(df$effect) &
      pv < fdr_cut_d() & abs(df$effect) >= fc_cut_d()
    keep <- intersect(c("comparison", "feature_id", "feature_symbol", "effect",
                        "effect_type", "statistic", "p_value", "adj_p_value",
                        "base_mean", "significant", "method"), names(df))
    df[, keep, drop = FALSE]
  }
  table_name <- function(ext) {
    sprintf("differential_%s_%s.%s", active()$tag %||% "layer",
            format(Sys.Date(), "%Y%m%d"), ext)
  }
  output$download_csv <- shiny::downloadHandler(
    filename = function() table_name("csv"),
    content = function(file) utils::write.csv(full_table(), file, row.names = FALSE,
                                              fileEncoding = "UTF-8")
  )
  output$download_xlsx <- shiny::downloadHandler(
    filename = function() table_name("xlsx"),
    content = function(file) {
      df <- full_table()
      sheets <- split(df, df$comparison %||% "result")
      names(sheets) <- substr(gsub("[\\[\\]*?/:]", "_", names(sheets)), 1L, 31L)
      openxlsx::write.xlsx(sheets, file)
    }
  )
  invisible()
}

# Which Bioconductor diff backends aren't installed in this R
# session. Returned as a character vector for the notices strip.
# DESeq2 / edgeR / limma are the ones run_diff() can dispatch to.
diff_missing_engines <- function() {
  engines <- c(limma = "limma", DESeq2 = "DESeq2", edgeR = "edgeR")
  missing <- vapply(engines, function(pkg) !has_pkg(pkg), logical(1))
  names(engines)[missing]
}

diff_volcano_card <- function(ns) {
  bslib::card(
    bslib::card_header(
      htmltools::tags$h3(class = "card-title", "Volcano"),
      # The cut the figure is drawn at, in words: plotly does not show
      # the caption plot_volcano() writes it into.
      shiny::textOutput(ns("volcano_cut"), container = function(...)
        htmltools::tags$span(class = "card-sub", ...))
    ),
    bslib::card_body(
      plotly::plotlyOutput(ns("volcano"), height = "360px"),
      shiny::uiOutput(ns("volcano_legend"))
    )
  )
}

diff_hits_card <- function(ns) {
  bslib::card(
    bslib::card_header(
      htmltools::tags$h3(class = "card-title", "Top hits"),
      htmltools::tags$span(class = "card-sub",
                           "largest changes first, within current thresholds")
    ),
    bslib::card_body(
      shiny::uiOutput(ns("hits_empty")),
      DT::DTOutput(ns("hits"), fill = FALSE),
      htmltools::tags$div(
        style = "display:flex;gap:8px;margin-top:8px;flex-wrap:wrap",
        shiny::downloadButton(ns("download_csv"), "Download all results (CSV)",
                              class = "btn btn-sm btn-ghost"),
        shiny::downloadButton(ns("download_xlsx"), "Excel, one sheet per comparison",
                              class = "btn btn-sm btn-ghost")
      )
    )
  )
}

# An R-internal error (a subscript, a missing object) rather than a
# sentence written for the user.
looks_internal_error <- function(msg) {
  grepl(paste0("Error in|subscript|object '.*' not found|non-numeric argument|",
               "argument is of length zero|missing value where|unused argument|",
               "could not find function|values must be length|replacement has|",
               "non-conformable|NA/NaN/Inf|invalid 'type'|undefined columns"),
        msg %||% "")
}

# A plain-language reading of the errors a run most often ends in; the
# message itself stays available under "Technical details".
diff_error_hint <- function(msg) {
  msg <- msg %||% ""
  if (grepl("not a level of", msg, fixed = TRUE)) {
    return("One of the chosen groups is not in this layer's group column.")
  }
  if (grepl("confounded", msg, fixed = TRUE)) {
    return("A covariate cannot be separated from the groups being compared; remove it.")
  }
  if (grepl("missing values", msg, fixed = TRUE)) {
    return("A covariate or pairing column has empty cells for some samples.")
  }
  if (grepl("paired design", msg, fixed = TRUE)) {
    return("The pairing column does not pair the samples one-to-one across the groups.")
  }
  if (grepl("required|not installed", msg)) {
    return("The chosen method needs a package that is not installed; pick another method.")
  }
  "See the technical details below."
}

# The volcano's significance cut, from the view's controls, in the form
# plot_volcano() takes it. A |log2FC| cutoff of 0 is no cutoff: drawn,
# it was a pair of dashed lines on top of the zero line.
volcano_p_basis <- function(p_col) {
  if (identical(p_col, "p_value")) "raw" else "adjusted"
}
volcano_effect_cut <- function(fc_cut) {
  if (is.null(fc_cut) || !is.finite(fc_cut) || fc_cut <= 0) NULL else fc_cut
}

# "significant = adjusted p < 0.05 and |log2FC| >= 0.263", in the words
# the controls use.
volcano_cut_text <- function(p_label, p_cut, effect_name, fc_cut) {
  txt <- sprintf("significant = %s < %s", p_label, format(p_cut, digits = 3))
  if (!is.null(volcano_effect_cut(fc_cut))) {
    txt <- sprintf("%s and |%s| \u2265 %s", txt, effect_name, format(fc_cut, digits = 3))
  }
  txt
}

# The volcano's legend, in plot_volcano()'s classes and colours: up in
# red, down in blue, the rest grey. Counts when there is a result.
volcano_legend <- function(up_n = NULL, down_n = NULL, continuous = FALSE) {
  words <- if (isTRUE(continuous)) c("positive", "negative") else c("up", "down")
  count <- function(word, n) {
    if (is.null(n)) word else sprintf("%s (%s)", word, format(n, big.mark = ","))
  }
  htmltools::tags$div(
    class = "legend",
    legend_swatch(count(words[[1]], up_n), omics_colors$up),
    legend_swatch(count(words[[2]], down_n), omics_colors$down),
    legend_swatch("not significant", omics_colors$ns)
  )
}

# ---- labelling the volcano's top features ------------------------------
#
# plotly has no label repulsion, and the top hits of a real result sit
# together at the top of the cloud: twenty labels at one fixed offset
# from their points piled into an unreadable block. So the labels are
# laid out here, in pixels, before plotly sees them. Up-regulated
# features are labelled to the right of their points and down-regulated
# to the left; on each side the labels are stacked in a column, one
# line apart, in the order of their points' height. Two labels on one
# side never share a line, and the two sides never cross the middle,
# so no two labels overlap. A label that cannot be fitted is left off
# -- the least significant on its side first -- and its name is still
# in the point's hover.

VOLCANO_LABEL_FONT <- 11
VOLCANO_LABEL_LINE <- 16   # px between stacked labels (11 px text + air)
VOLCANO_LABEL_GAP <- 14    # px from a point to its label
VOLCANO_PLOT_HEIGHT <- 360 # the card's plotlyOutput height

# Approximate width of a label in px at the label font.
volcano_label_width <- function(text) nchar(text) * 6.5 + 4

# The plot area in px, from the figure ggplotly() built and the width
# the browser reports for the output (assumed when it reports none).
volcano_geometry <- function(fig, width = NULL) {
  m <- fig$x$layout$margin
  width <- if (is.numeric(width) && length(width) == 1L && is.finite(width) && width > 0)
    width else 600
  list(
    x_range = unlist(fig$x$layout$xaxis$range),
    y_range = unlist(fig$x$layout$yaxis$range),
    # A little under the margins ggplotly asked for: the axis titles
    # can take a few px more, and a plot area smaller than assumed only
    # spreads the labels further apart than needed.
    w = max(width - (m$l %||% 50) - (m$r %||% 10) - 8, 120),
    h = max(VOLCANO_PLOT_HEIGHT - (m$t %||% 40) - (m$b %||% 40) - 8, 120))
}

# Stack labels wanted at heights `want` (px from the top of the plot
# area) into lines at least `line` apart within [lo, hi], keeping their
# order. Drops the least significant (highest `rank`) until they fit.
# Returns the kept indices and their heights.
#
# A label pushed far from its point (more than `max_shift` px up or
# down) counts as not fitting too: a column of names reaching to the
# x axis, each on a long line across the cloud, was technically clear
# of overlaps and still hard to read.
stack_volcano_labels <- function(want, rank, lo, hi, line, max_shift = Inf) {
  keep <- seq_along(want)
  repeat {
    if (!length(keep)) return(list(keep = integer(0), y = numeric(0)))
    o <- keep[order(want[keep], rank[keep])]
    y <- pmin(pmax(want[o], lo), hi)
    for (i in seq_along(y)[-1L]) y[i] <- max(y[i], y[i - 1L] + line)
    if (y[length(y)] > hi) {
      y[length(y)] <- hi
      for (i in rev(seq_along(y))[-1L]) y[i] <- min(y[i], y[i + 1L] - line)
    }
    if (y[1L] >= lo - 1e-9 && all(abs(y - want[o]) <= max_shift)) {
      return(list(keep = o, y = y))
    }
    keep <- setdiff(keep, keep[which.max(rank[keep])])
  }
}

# The `n` most significant features of a result as plotly annotations,
# laid out as above for a plot area of `geom` (see volcano_geometry()).
# Their heads sit on the points, at the coordinates plot_volcano()
# draws them (the p column the figure was drawn from); the label text
# sits at the end of a thin leader line. Returns the annotations and
# the x range to draw at: wider than the data when that makes room for
# the labels beside the outermost points, but never so wide that the
# data is squeezed into less than ~40% of the plot.
volcano_annotations <- function(bundle, n, p_col = "adj_p_value", geom) {
  none <- list(annotations = list(), x_range = geom$x_range)
  df <- bundle$results$diff_result_df
  df <- df[!is.na(df[[p_col]]) & !is.na(df$effect), , drop = FALSE]
  top <- utils::head(df[order(df[[p_col]]), , drop = FALSE], n)
  if (!nrow(top)) return(none)
  lab <- ifelse(is.na(top$feature_symbol) | !nzchar(top$feature_symbol),
                top$feature_id, top$feature_symbol)
  x <- top$effect
  y <- -log10(pmax(top[[p_col]], .Machine$double.xmin))
  rank <- seq_len(nrow(top))
  right <- x >= 0
  lw <- volcano_label_width(lab)
  W <- geom$w
  H <- geom$h

  # Room for each side's column beside its outermost point: on a
  # symmetric axis of half-width R, a point at |x| leaves room for a
  # label of width w when |x|/(2R) * W + gap + w <= W/2.
  half <- max(abs(geom$x_range))
  need <- vapply(list(right, !right), function(s) {
    if (!any(s)) return(0)
    max(abs(x[s])) * W / max(W - 2 * (VOLCANO_LABEL_GAP + max(lw[s])), 1)
  }, numeric(1))
  # Capped so the data keeps at least ~40% of the width; a phone gets
  # the narrower cloud, and the toggle is the user asking for names.
  half <- min(max(half, need), half * 2.4)
  x_range <- c(-half, half)
  px <- (x + half) / (2 * half) * W
  y0 <- geom$y_range[[1L]]
  y1 <- geom$y_range[[2L]]
  py <- (1 - (y - y0) / (y1 - y0)) * H

  # Each side's labels line up in one column just beyond its outermost
  # labelled point (right side: their left edges; left side: their
  # right edges), so no label sits on a point it names. Where the plot's
  # edge leaves no room beyond a point, that point goes unlabelled
  # (outermost first) and the column is placed for the rest; a column
  # that would cross the middle is left off.
  edge <- numeric(length(x))
  fits <- logical(length(x))
  for (side in c(TRUE, FALSE)) {
    cand <- which(right == side)
    while (length(cand)) {
      e <- if (side) min(max(px[cand]) + VOLCANO_LABEL_GAP, W - max(lw[cand]))
           else max(min(px[cand]) - VOLCANO_LABEL_GAP, max(lw[cand]))
      clear <- if (side) px[cand] + VOLCANO_LABEL_GAP <= e + 1e-9
               else px[cand] - VOLCANO_LABEL_GAP >= e - 1e-9
      if (all(clear)) break
      # The outermost point that the column cannot get past.
      cand <- setdiff(cand, cand[!clear][which.max(abs(px[cand][!clear] - W / 2))])
    }
    if (length(cand)) {
      edge[cand] <- e
      fits[cand] <- if (side) e >= W / 2 + 2 else e <= W / 2 - 2
    }
  }

  out <- list()
  for (side in c(TRUE, FALSE)) {
    idx <- which(right == side & fits)
    if (!length(idx)) next
    st <- stack_volcano_labels(py[idx], rank[idx], lo = VOLCANO_LABEL_LINE / 2,
                               hi = H - VOLCANO_LABEL_LINE / 2, line = VOLCANO_LABEL_LINE,
                               max_shift = H * 0.4)
    for (j in seq_along(st$keep)) {
      i <- idx[[st$keep[[j]]]]
      out[[length(out) + 1L]] <- list(
        rank = rank[[i]],
        x = x[[i]], y = y[[i]], text = lab[[i]],
        # The label's anchor, in data coordinates (axref/ayref on the
        # axes), so it stays put relative to the points if the plot is
        # a few px off the size assumed here.
        ax = edge[[i]] / W * (2 * half) - half,
        ay = y1 - st$y[[j]] / H * (y1 - y0),
        axref = "x", ayref = "y",
        xanchor = if (side) "left" else "right", yanchor = "middle",
        showarrow = TRUE, arrowhead = 0, arrowwidth = 0.7, standoff = 3,
        arrowcolor = "#9AA3AE",
        font = list(size = VOLCANO_LABEL_FONT, color = "#1A2541"))
    }
  }
  # Most significant first, as they were chosen.
  out <- out[order(vapply(out, `[[`, numeric(1), "rank"))]
  out <- lapply(out, function(a) { a$rank <- NULL; a })
  list(annotations = out, x_range = x_range)
}

# Adds the labels to a ggplotly() volcano, widening its x range when
# that gives them room.
label_volcano <- function(fig, bundle, n, p_col, width = NULL) {
  geom <- volcano_geometry(fig, width)
  lab <- volcano_annotations(bundle, n, p_col, geom)
  if (!length(lab$annotations)) return(fig)
  plotly::layout(fig, annotations = lab$annotations,
                 xaxis = list(range = lab$x_range))
}

# ggplotly() sets `hoveron` on its traces; scattergl has no such
# attribute, so once toWebGL() converts them plotly warns about it on
# every build. Dropped first, it is never there to warn about.
ggplotly_untitled <- function() ggplot2::labs(title = NULL, subtitle = NULL)

drop_hoveron <- function(fig) {
  fig$x$data <- lapply(fig$x$data, function(tr) {
    tr$hoveron <- NULL
    tr
  })
  fig
}
