# Hover read-outs for the static figures.
#
# The PCA and the integration scatters are ggplots drawn as images, and
# a reader looking at an outlying sample or a gene off the diagonal had
# no way to ask which one it was: the names drawn on the plot are only
# the top few, and the rest of the points were anonymous. Hovering a
# point (or tapping it, on a phone) now shows a small card with its
# name and values.
#
# The figures stay what they are -- ggplot through renderPlot -- rather
# than turning into plotly widgets: Shiny already reports where the
# pointer is in the plot's own coordinates, and nearPoints() finds the
# point there, so the image, its labels and its legend are unchanged.
#
# Use:
#   UI      hover_plot_output(ns("pca"), height = "360px")
#   server  pca_plot <- shiny::reactive(<the ggplot>)
#           output$pca <- shiny::renderPlot(fit_to_width("pca", pca_plot()))
#           plot_hover_server("pca", pca_plot, pca_hover_text, input, output, session)
#
# `describe(row, p)` turns the nearest row of the plot's data into
# list(title =, rows = c(key = value, ...), note =), or NULL for no card.

# Debounced: the card follows the pointer once it rests, not on every
# pixel it crosses, which would be a server round trip per mouse move.
PLOT_HOVER_DELAY_MS <- 80L

# How near, in screen pixels, the pointer has to be to a point. A finger
# is less precise than a mouse, so a tap reaches further.
PLOT_HOVER_RADIUS_PX <- c(hover = 10, click = 18)

# The card's widest. Narrower where the plot leaves less room.
PLOT_HOVER_MAX_WIDTH_PX <- 260L

#' A plot output whose points can be hovered or tapped
#'
#' A `plotOutput()` that reports the pointer (`<id>_hover`) and taps or
#' clicks (`<id>_click`), with the read-out layer (`<id>_tip`) laid over
#' it. `output_id` is the namespaced id.
#' @keywords internal
#' @noRd
hover_plot_output <- function(output_id, height) {
  htmltools::tags$div(
    class = "plot-hover-wrap",
    shiny::plotOutput(
      output_id, height = height,
      hover = shiny::hoverOpts(paste0(output_id, "_hover"),
                               delay = PLOT_HOVER_DELAY_MS,
                               delayType = "debounce", nullOutside = TRUE),
      click = shiny::clickOpts(paste0(output_id, "_click"))
    ),
    # A pointer-only aid: the plot's alt text and the table beside it
    # carry the same names for a keyboard or screen-reader user, so the
    # card is hidden from assistive technology and never takes focus.
    shiny::uiOutput(paste0(output_id, "_tip"), class = "plot-hover-layer",
                    `aria-hidden` = "true")
  )
}

#' Serve the read-out for one hover_plot_output()
#'
#' @param output_id The un-namespaced id of the plot output.
#' @param plot A reactive returning the ggplot that output draws (before
#'   fit_to_width(), which changes only its theme).
#' @param describe function(row, p) -> list(title, rows, note) or NULL.
#' @keywords internal
#' @noRd
plot_hover_server <- function(output_id, plot, describe, input, output, session) {
  hover_id <- paste0(output_id, "_hover")
  click_id <- paste0(output_id, "_click")
  # The latest of the two: a mouse moving over the plot, or a tap. The
  # hover goes back to NULL when the pointer leaves the plot, and the
  # card with it.
  pointer <- shiny::reactiveVal(NULL)
  shiny::observeEvent(input[[hover_id]], {
    ev <- input[[hover_id]]
    pointer(if (!is.null(ev)) list(event = ev, kind = "hover"))
  }, ignoreNULL = FALSE, ignoreInit = TRUE)
  # Shiny sets both back to NULL when a new image replaces the old one
  # (another colouring, a re-run): a tapped card would otherwise stay up
  # over a plot whose points have moved.
  shiny::observeEvent(input[[click_id]], {
    ev <- input[[click_id]]
    if (!is.null(ev)) pointer(list(event = ev, kind = "click"))
    else if (identical(shiny::isolate(pointer())$kind, "click")) pointer(NULL)
  }, ignoreNULL = FALSE, ignoreInit = TRUE)

  output[[paste0(output_id, "_tip")]] <- shiny::renderUI({
    ptr <- pointer()
    if (is.null(ptr)) return(NULL)
    p <- tryCatch(plot(), error = function(e) NULL)
    row <- plot_hover_row(p, ptr$event, PLOT_HOVER_RADIUS_PX[[ptr$kind]])
    if (is.null(row)) return(NULL)
    text <- tryCatch(describe(row, p), error = function(e) NULL)
    if (is.null(text)) return(NULL)
    cd <- session$clientData
    size <- function(what) {
      v <- tryCatch(cd[[sprintf("output_%s_%s", session$ns(output_id), what)]],
                    error = function(e) NULL)
      if (is.numeric(v) && length(v) == 1L && is.finite(v)) v else NA_real_
    }
    plot_hover_card(text, ptr$event, width = size("width"), height = size("height"))
  })
  invisible(pointer)
}

#' The row of a ggplot's data nearest a hover or click event
#'
#' The plot's own data and its own x and y columns, so the point found is
#' the one drawn. NULL when nothing is within `radius` CSS pixels, or the
#' plot has no point data (an empty-result placeholder).
#' @keywords internal
#' @noRd
plot_hover_row <- function(p, event, radius = PLOT_HOVER_RADIUS_PX[["hover"]]) {
  if (is.null(event) || is.null(event$x) || !inherits(p, "ggplot")) return(NULL)
  df <- p$data
  if (!is.data.frame(df) || !nrow(df)) return(NULL)
  xvar <- plot_aes_column(p, "x")
  yvar <- plot_aes_column(p, "y")
  if (is.null(xvar) || is.null(yvar) || !all(c(xvar, yvar) %in% names(df))) return(NULL)
  hit <- tryCatch(
    shiny::nearPoints(df, event, xvar = xvar, yvar = yvar,
                      threshold = radius, maxpoints = 1L),
    error = function(e) NULL)
  if (is.null(hit) || !nrow(hit)) return(NULL)
  hit
}

#' The data column a ggplot maps to an aesthetic
#'
#' `.data$PC1`, `.data[["effect"]]` or a bare `PC1` all name a column;
#' anything computed (`-log10(p)`) does not, and gives NULL.
#' @keywords internal
#' @noRd
plot_aes_column <- function(p, aes) {
  q <- p$mapping[[aes]]
  if (is.null(q)) return(NULL)
  expr <- if (inherits(q, "formula")) unclass(q)[[2L]] else q
  env <- attr(q, ".Environment") %||% emptyenv()
  if (is.symbol(expr)) return(as.character(expr))
  if (is.call(expr) && length(expr) == 3L && identical(expr[[2L]], quote(.data))) {
    if (identical(expr[[1L]], as.name("$"))) return(as.character(expr[[3L]]))
    if (identical(expr[[1L]], as.name("[["))) {
      col <- tryCatch(eval(expr[[3L]], env), error = function(e) NULL)
      if (is.character(col) && length(col) == 1L) return(col)
    }
  }
  NULL
}

#' The read-out card, next to the pointer and inside the plot
#'
#' It opens away from the pointer towards the larger side of the plot,
#' right or left and below or above, and never wider than the room left
#' on that side: a bslib card scrolls whatever runs past its edge, and
#' on a 300 px phone plot a fixed-width card would have.
#' @keywords internal
#' @noRd
plot_hover_card <- function(text, event, width = NA_real_, height = NA_real_) {
  x <- as.numeric(event$coords_css$x %||% NA_real_)
  y <- as.numeric(event$coords_css$y %||% NA_real_)
  if (!is.finite(x) || !is.finite(y)) return(NULL)
  gap <- 12
  edge <- 4
  left_side <- !is.finite(width) || x <= width / 2
  horizontal <- if (left_side) {
    room <- if (is.finite(width)) width - x - gap - edge else PLOT_HOVER_MAX_WIDTH_PX
    sprintf("left:%dpx;max-width:%dpx", round(x + gap),
            round(min(PLOT_HOVER_MAX_WIDTH_PX, room)))
  } else {
    sprintf("right:%dpx;max-width:%dpx", round(width - x + gap),
            round(min(PLOT_HOVER_MAX_WIDTH_PX, x - gap - edge)))
  }
  vertical <- if (is.finite(height) && y > height / 2) {
    sprintf("bottom:%dpx", round(height - y + gap))
  } else {
    sprintf("top:%dpx", round(y + gap))
  }
  rows <- text$rows %||% character(0)
  rows <- rows[!is.na(rows)]
  htmltools::tags$div(
    class = "plot-hover-tip",
    style = paste(horizontal, vertical, sep = ";"),
    htmltools::tags$div(class = "plot-hover-title", text$title),
    if (length(rows)) htmltools::tags$table(
      class = "plot-hover-rows",
      htmltools::tags$tbody(lapply(seq_along(rows), function(i) {
        # Numbers line up on the right in the monospace face; words (a
        # group's name) read as text.
        num <- grepl("^[-+]?[0-9.]+(e[-+]?[0-9]+)?$", rows[[i]])
        htmltools::tags$tr(htmltools::tags$th(names(rows)[[i]]),
                           htmltools::tags$td(class = if (num) "num", rows[[i]]))
      }))
    ),
    if (!is.null(text$note) && nzchar(text$note))
      htmltools::tags$div(class = "plot-hover-note", text$note)
  )
}

# ---- number formats ---------------------------------------------------

hover_num <- function(x, digits = 2L) {
  x <- suppressWarnings(as.numeric(x))
  if (!length(x) || is.na(x[[1L]])) return("—")
  # round() first: -0.03 at one digit is "0.0", not "-0.0".
  x <- round(x[[1L]], digits) + 0
  formatC(x, format = "f", digits = digits)
}

# 0.012, 0.0034, 3.2e-05: two significant digits, scientific below 1e-3
# as in the result tables.
hover_p <- function(p) {
  p <- suppressWarnings(as.numeric(p))
  if (!length(p) || is.na(p[[1L]])) return("—")
  p <- p[[1L]]
  if (p < 1e-3) formatC(p, format = "e", digits = 1L)
  else format(signif(p, 2L), scientific = FALSE)
}

# ---- what each figure's card says --------------------------------------

#' PCA: the sample, its group, and where it sits
#' @keywords internal
#' @noRd
pca_hover_text <- function(row, p) {
  group_col <- plot_aes_column(p, "colour")
  rows <- character(0)
  if (!is.null(group_col) && group_col %in% names(row)) {
    g <- as.character(row[[group_col]][[1L]])
    rows[[group_col]] <- if (is.na(g)) "—" else g
  }
  rows[["PC1"]] <- hover_num(row$PC1, 1L)
  rows[["PC2"]] <- hover_num(row$PC2, 1L)
  list(title = as.character(row$sample_id[[1L]]), rows = rows)
}

# The gene, and where a gene has several protein-gene pairs (several
# points) the pair's id as well, which is how the result table tells
# them apart: "TP53 (P04637-2)".
hover_feature_name <- function(row, p) {
  sym <- as.character(row$feature_symbol[[1L]])
  id <- as.character(row$feature_id[[1L]])
  if (is.na(sym) || !nzchar(sym)) return(id)
  if (sum(p$data$feature_symbol %in% sym) <= 1L || is.na(id)) return(sym)
  if (grepl(sym, id, fixed = TRUE)) id else sprintf("%s (%s)", sym, id)
}

# An axis title as drawn: labs() sets it; otherwise the column's name.
hover_axis_title <- function(p, aes) {
  lab <- tryCatch(p$labels[[aes]], error = function(e) NULL)
  if (is.character(lab) && length(lab) == 1L && nzchar(lab)) lab
  else plot_aes_column(p, aes) %||% aes
}

# A colour class in the legend's words, without the legend's count:
# "up in both (14)" is "up in both".
hover_legend_word <- function(p, value) {
  value <- as.character(value)
  sc <- tryCatch(p$scales$get_scales("colour"), error = function(e) NULL)
  brk <- tryCatch(sc$breaks, error = function(e) NULL)
  lab <- tryCatch(sc$labels, error = function(e) NULL)
  if (is.character(brk) && is.character(lab) && length(brk) == length(lab) &&
      value %in% brk) {
    return(sub(" \\([0-9,]+\\)$", "", lab[[match(value, brk)]]))
  }
  value
}

#' Effect in each layer: the gene, both log2FCs, and its class
#' @keywords internal
#' @noRd
effect_pair_hover_text <- function(row, p) {
  rows <- character(0)
  rows[[hover_axis_title(p, "x")]] <- hover_num(row$effect_a)
  rows[[hover_axis_title(p, "y")]] <- hover_num(row$effect_b)
  cls <- as.character(row$.class[[1L]])
  note <- if (is.na(cls) || identical(cls, "background")) {
    "not a hit in both layers"
  } else {
    paste("hit in both layers ·", hover_legend_word(p, cls))
  }
  list(title = hover_feature_name(row, p), rows = rows, note = note)
}

#' Correlation per gene: the gene, its r and adjusted p
#' @keywords internal
#' @noRd
correlation_hover_text <- function(row, p) {
  type <- as.character(row$effect_type[[1L]] %||% NA_character_)
  key <- switch(if (is.na(type)) "" else type,
                spearman_r = "Spearman r", pearson_r = "Pearson r",
                "correlation")
  rows <- character(0)
  rows[[key]] <- hover_num(row$effect)
  rows[["adjusted p"]] <- hover_p(row$adj_p_value)
  list(title = hover_feature_name(row, p), rows = rows)
}

#' Top hits: the gene and its effect in each layer
#'
#' The dumbbell's data is one row per gene and layer, so the card
#' gathers the gene's rows whichever of its dots is hovered.
#' @keywords internal
#' @noRd
top_hits_hover_text <- function(row, p) {
  label <- as.character(row$.label[[1L]])
  d <- p$data[as.character(p$data$.label) == label, , drop = FALSE]
  eff <- hover_axis_title(p, "x")
  rows <- vapply(d$effect, hover_num, character(1))
  names(rows) <- sprintf("%s (%s)", eff, as.character(d$layer))
  list(title = label, rows = rows)
}

#' ActivePathways: the pathway in full, its adjusted p and direction
#' @keywords internal
#' @noRd
active_pathways_hover_text <- function(row, p) {
  rows <- c(`adjusted p` = hover_p(row$adj_p_value))
  if (".dir" %in% names(row)) rows[["direction"]] <- as.character(row$.dir[[1L]])
  ev <- as.character(row$.evidence[[1L]] %||% NA_character_)
  found <- c(shared = "both layers", unique = "one layer",
             combined = "only combined")[ev]
  if (!is.na(found)) rows[["found by"]] <- unname(found)
  list(title = hover_pathway_name(row$feature_symbol[[1L]]), rows = rows)
}

# The pathway as the dot plot's axis names it, but whole: the axis drops
# the collection prefix and the underscores and then shortens a long
# name; the card shows the full name. The raw id,
# "HALLMARK_INFLAMMATORY_RESPONSE", did not match the row the reader
# pointed at, and with no spaces it broke mid-word on a phone.
hover_pathway_name <- function(x) {
  x <- as.character(x)
  if (is.na(x)) return("—")
  pretty <- tryCatch(utils::getFromNamespace("prettify_gene_set_name", "omicsCore"),
                     error = function(e) function(x) trimws(gsub("_", " ", x, fixed = TRUE)))
  pretty(x)
}
