# A download menu for every figure: PNG at 300 dpi, PDF and SVG.
#
# Biologists take figures to slides and papers, and the only way out of
# the app was a right-click on a 96 dpi screen image -- too coarse to
# print, and drawn at whatever width the card had, on a phone the shrunk
# phone layout. Each figure card now carries a small download button in
# its header, and the file is drawn afresh from the figure's own ggplot
# at a print size, with the desktop layout.
#
# Use (a figure built by a reactive, as every card's figure is):
#   UI      bslib::card_header(<title>, <sub>,
#                              plot_download_ui(ns("pca_download")))
#   server  plot_download_server("pca_download", pca_plot,
#                                filename = function() c(project, "pca", layer),
#                                width_in = 7, height_in = 5)
#
# The menu shows only while there is a figure to save: before a result,
# in the demo, or when a figure is only a sentence saying why it is
# empty, the header has no button rather than one that saves an error.
# The wiring test (test-ui-server-wiring.R) checks that every figure an
# app module draws has its menu, "<output id>_download".

# The formats offered, in the menu's order, with what the menu says.
PLOT_DOWNLOAD_FORMATS <- c(png = "PNG (300 dpi)", pdf = "PDF", svg = "SVG")

PLOT_DOWNLOAD_DPI <- 300

# A figure's size when it says nothing else: a full page width less the
# margins, at the proportions of the cards.
PLOT_DOWNLOAD_WIDTH_IN <- 7
PLOT_DOWNLOAD_HEIGHT_IN <- 4.5

#' The download button of a figure card
#'
#' A placeholder that the server fills with the menu while there is a
#' figure to save (see plot_download_server()). Goes in the card header,
#' after the title; styles.scss pushes it to the right.
#'
#' @param id The namespaced id of the download, `ns("<output>_download")`.
#' @keywords internal
#' @noRd
plot_download_ui <- function(id) {
  shiny::uiOutput(shiny::NS(id, "menu"), inline = TRUE, class = "plot-download")
}

#' The menu itself: an icon button that opens the three formats
#'
#' A Bootstrap dropdown: Enter, Space or the arrow keys open it from the
#' keyboard and the arrows move between the formats; Escape closes it.
#' It opens with a fixed position, so a card that clips its contents
#' (every bslib card does) does not cut it off.
#' @keywords internal
#' @noRd
plot_download_menu <- function(ns, label = "Download figure") {
  htmltools::tags$div(
    class = "dropdown plot-download-menu",
    htmltools::tags$button(
      type = "button",
      class = "btn btn-sm plot-download-toggle",
      `data-bs-toggle` = "dropdown",
      `data-bs-popper-config` = '{"strategy":"fixed"}',
      `aria-expanded` = "false",
      `aria-haspopup` = "true",
      `aria-label` = label,
      title = label,
      bsicons::bs_icon("download", a11y = "deco")
    ),
    htmltools::tags$ul(
      class = "dropdown-menu dropdown-menu-end plot-download-items",
      htmltools::tags$li(htmltools::tags$h6(class = "dropdown-header", label)),
      lapply(names(PLOT_DOWNLOAD_FORMATS), function(fmt) {
        htmltools::tags$li(shiny::downloadLink(
          ns(fmt), PLOT_DOWNLOAD_FORMATS[[fmt]], class = "dropdown-item"))
      })
    )
  )
}

#' Serve a figure card's downloads
#'
#' Registers `png`, `pdf` and `svg` download handlers under `id`, and the
#' menu that offers them.
#'
#' @param id The un-namespaced id the UI used, `"<output>_download"`.
#' @param plot_reactive A reactive returning the figure's ggplot (or
#'   patchwork) -- the one its card draws, before fit_to_width(): a file
#'   is drawn with the desktop layout even when the user is on a phone.
#' @param filename_stem A function returning the parts of the file name,
#'   e.g. `c(project, "volcano", comparison)`. Each part is made safe for a
#'   file name (project_slug(), as for the report's downloads), empty parts
#'   are dropped, and the parts are joined by "_".
#' @param width_in,height_in The size of the file in inches: a number, or
#'   a function returning one (a figure as tall as its rows).
#' @param available Optional reactive: FALSE hides the menu. By default
#'   the menu shows while `plot_reactive()` gives a figure that is more
#'   than an "empty" sentence (plot_ready()). Give it when building the
#'   figure only to check costs too much (the volcano), or to hide the
#'   menu in the demo.
#' @param label What the button says to a screen reader and in its
#'   tooltip; a card with two figures names which one.
#' @return The availability reactiveVal, invisibly.
#' @keywords internal
#' @noRd
plot_download_server <- function(id, plot_reactive, filename_stem,
                                 width_in = PLOT_DOWNLOAD_WIDTH_IN,
                                 height_in = PLOT_DOWNLOAD_HEIGHT_IN,
                                 available = NULL, label = "Download figure") {
  shiny::moduleServer(id, function(input, output, session) {
    # Set by an observer, and only when it changes, so the menu is not
    # drawn again (and closed under the user's pointer) every time the
    # figure is.
    shown <- shiny::reactiveVal(FALSE)
    shiny::observe({
      ok <- isTRUE(tryCatch(
        if (is.null(available)) plot_ready(plot_reactive) else available(),
        error = function(e) FALSE))
      if (!identical(ok, shiny::isolate(shown()))) shown(ok)
    })

    output$menu <- shiny::renderUI({
      if (!shown()) return(NULL)
      plot_download_menu(session$ns, label)
    })
    # Drawn even while its placeholder is hidden: an empty placeholder can
    # read as hidden, and would then never be filled.
    shiny::outputOptions(output, "menu", suspendWhenHidden = FALSE)

    lapply(names(PLOT_DOWNLOAD_FORMATS), function(fmt) {
      output[[fmt]] <- shiny::downloadHandler(
        filename = function() {
          paste0(plot_download_name(tryCatch(filename_stem(), error = function(e) NULL)),
                 ".", fmt)
        },
        content = function(file) {
          p <- plot_reactive()
          shiny::req(is_drawable_plot(p))
          save_plot_file(p, file, fmt,
                         width_in = plot_download_size(width_in, PLOT_DOWNLOAD_WIDTH_IN),
                         height_in = plot_download_size(height_in, PLOT_DOWNLOAD_HEIGHT_IN))
        }
      )
      # The links sit in a closed menu, which Shiny reads as hidden; a
      # suspended link gets no address until the menu opens, and is
      # disabled -- skipped by the arrow keys -- for the moment after.
      shiny::outputOptions(output, fmt, suspendWhenHidden = FALSE)
    })
    invisible(shown)
  })
}

#' Draw a figure to a file
#'
#' PNG at 300 dpi (ragg when it is installed, otherwise cairo), PDF
#' through cairo, which embeds its fonts, and SVG through cairo, which
#' turns text into outlines -- svglite would keep it as text but is not
#' installed. Cairo is the one that draws "≥", "·" and Chinese
#' labels, from whichever installed font has them; the plain pdf() device
#' replaced them with dots. The background is white, not transparent: a
#' transparent PNG on a dark slide loses its black text.
#' @keywords internal
#' @noRd
save_plot_file <- function(p, file, format = c("png", "pdf", "svg"),
                           width_in = PLOT_DOWNLOAD_WIDTH_IN,
                           height_in = PLOT_DOWNLOAD_HEIGHT_IN) {
  format <- match.arg(format)
  device <- switch(format,
                   png = plot_png_device(),
                   pdf = grDevices::cairo_pdf,
                   svg = grDevices::svg)
  ggplot2::ggsave(file, plot = p, device = device, width = width_in,
                  height = height_in, units = "in", dpi = PLOT_DOWNLOAD_DPI,
                  bg = "white", limitsize = FALSE)
  invisible(file)
}

# ragg draws text a little more evenly than cairo's PNG and is what
# ggplot2 itself prefers. Looked up rather than called with `::`: it is
# not a dependency, and the cairo PNG is a fine figure too.
plot_png_device <- function() {
  if (has_pkg("ragg")) return(getExportedValue("ragg", "agg_png"))
  function(filename, width, height, units = "in", res = PLOT_DOWNLOAD_DPI,
           bg = "white", ...) {
    grDevices::png(filename, width = width, height = height, units = units,
                   res = res, bg = bg, type = "cairo")
  }
}

#' Whether a reactive gives a figure worth saving
#'
#' FALSE while it errors or waits (`req()`), and for the "No pathways to
#' plot." kind of figure, which is a sentence and not a result.
#' @keywords internal
#' @noRd
plot_ready <- function(plot_reactive) {
  p <- tryCatch(plot_reactive(), error = function(e) NULL)
  is_drawable_plot(p) && !is_empty_plot(p)
}

is_drawable_plot <- function(p) inherits(p, c("ggplot", "patchwork"))

# omicsCore's empty_plot(): no data, no mapping, nothing but text
# annotations on a blank theme.
is_empty_plot <- function(p) {
  if (!inherits(p, "ggplot") || inherits(p, "patchwork")) return(FALSE)
  no_data <- is.null(p$data) || inherits(p$data, "waiver") ||
    (is.data.frame(p$data) && !nrow(p$data))
  layers <- p$layers
  no_data && !length(p$mapping) && length(layers) > 0L &&
    all(vapply(layers, function(l) isFALSE(l$inherit.aes) &&
                 inherits(l$geom, "GeomText"), logical(1)))
}

# A size given as a number or as a function returning one.
plot_download_size <- function(x, default) {
  v <- tryCatch(if (is.function(x)) x() else x, error = function(e) NULL)
  v <- suppressWarnings(as.numeric(v))
  if (length(v) != 1L || !is.finite(v) || v <= 0) default else v
}

#' The file name (without extension) from its parts
#'
#' `c("My project", "volcano", "TreatA_vs_Control")` gives
#' `"My_project_volcano_TreatA_vs_Control"`. Each part goes through
#' project_slug(), as the report's downloads do -- the store's rule, which
#' keeps Chinese names whole and takes out slashes and the characters
#' Windows refuses -- and the whole is kept under 150 bytes, which every
#' file system takes.
#' @keywords internal
#' @noRd
plot_download_name <- function(parts) {
  parts <- vapply(as.list(unlist(parts, use.names = FALSE)), function(x) {
    if (length(x) != 1L || is.na(x)) return("")
    s <- project_slug(as.character(x))
    if (is.na(s)) "" else s
  }, character(1))
  parts <- parts[nzchar(parts)]
  if (!length(parts)) return("figure")
  out <- paste(parts, collapse = "_")
  out <- truncate_bytes(out, 150L)
  out <- gsub("[._]+$", "", out)
  if (nzchar(out)) out else "figure"
}

# The project's name for a file name, or nothing for the demo (which has
# no download) and an unnamed project.
plot_download_project <- function(proj) {
  if (is.null(proj)) return(NULL)
  proj$name %||% NULL
}
