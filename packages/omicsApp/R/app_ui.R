#' Build the omicsApp top-level UI
#'
#' Phase 2 slice 2A scaffold. Returns a [bslib::page_sidebar()] with:
#'   * dark vertical nav listing 7 views (Workflow + Output groups),
#'   * a header with a project picker, omics-layer tabs, and action
#'     buttons,
#'   * a hidden `tabsetPanel` in the main pane — one tabPanel per
#'     view, each rendered once on mount so switching is instant.
#'
#' All view content for 2A is a placeholder card; real content lands
#' in slices 2C-2F.
#'
#' @return A `shiny.tag` page object.
#'
#' @keywords internal
#' @noRd
app_ui <- function() {
  bslib::page_sidebar(
    title = app_header(),
    theme = app_theme(),
    fillable = FALSE,
    window_title = "omicsApp",
    lang = "en",
    sidebar = bslib::sidebar(
      id = "main_sidebar",
      # Always open on a desktop; on a phone a toggle, so the navigation
      # is not stacked below the whole page (it sat at y = 2695).
      open = list(desktop = "always", mobile = "closed"),
      width = 232,
      bg = "#161A26",
      padding = c(14, 12, 14, 12),
      gap = 2,
      app_sidebar_content()
    ),
    # A spinner on any output that is recomputing and a pulse along the
    # top while the server is busy: a plot being redrawn no longer looks
    # like a plot that has stopped.
    shiny::useBusyIndicators(spinners = TRUE, pulse = TRUE),
    app_favicon(),
    shiny::tabsetPanel(
      id = "view",
      type = "hidden",
      selected = "project",
      shiny::tabPanelBody("project",     project_view_ui("project")),
      shiny::tabPanelBody("import",      import_view_ui("import")),
      shiny::tabPanelBody("qc",          qc_view_ui("qc")),
      shiny::tabPanelBody("diff",        diff_view_ui("diff")),
      shiny::tabPanelBody("enrich",      enrich_view_ui("enrich")),
      shiny::tabPanelBody("integration", integration_view_ui("integration")),
      shiny::tabPanelBody("report",      report_view_ui("report"))
    )
  )
}

# The brand hexagon as the tab icon, inline. Without a declared icon
# every page load asked for /favicon.ico, which the app does not serve:
# a 404 in every browser console and every server log.
app_favicon <- function() {
  svg <- paste0(
    "<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 32 32'>",
    "<defs><linearGradient id='g' x1='0' y1='0' x2='1' y2='1'>",
    "<stop offset='0' stop-color='#3D52B8'/><stop offset='1' stop-color='#1FBF9E'/>",
    "</linearGradient></defs>",
    "<path fill='url(#g)' d='M16 1.5 29 9v14L16 30.5 3 23V9z'/></svg>")
  htmltools::tags$head(htmltools::tags$link(
    rel = "icon", type = "image/svg+xml",
    href = paste0("data:image/svg+xml,", utils::URLencode(svg, reserved = TRUE))))
}

# ---- sidebar contents -------------------------------------------------

# Sidebar = brand block + nav list (Workflow / Output) + status footer.
# `nav_item()` emits an `actionLink` styled with .nav-item; the active
# class is toggled from the server via `shinyjs::removeClass()` /
# `addClass()` whenever the current view changes.
app_sidebar_content <- function() {
  htmltools::tagList(
    htmltools::tags$div(
      class = "sidebar-brand",
      htmltools::tags$div(
        class = "app-brand-logo",
        htmltools::HTML("&#x2B22;")  # hexagon
      ),
      htmltools::tags$div(
        htmltools::tags$div(class = "app-brand-name", "omicsApp"),
        htmltools::tags$div(
          class = "app-brand-version",
          paste0("v", utils::packageVersion("omicsApp"))
        )
      )
    ),

    htmltools::tags$div(class = "nav-section", "Workflow"),
    nav_item("nav_project",     "Project",        bsicons::bs_icon("house"),            active = TRUE),
    nav_item("nav_import",      "Import",         bsicons::bs_icon("upload")),
    nav_item("nav_qc",          "Quality Control",bsicons::bs_icon("check-circle")),
    nav_item("nav_diff",        "Differential",   bsicons::bs_icon("graph-up-arrow")),
    nav_item("nav_enrich",      "Enrichment",     bsicons::bs_icon("diagram-3")),
    nav_item("nav_integration", "Integration",    bsicons::bs_icon("intersect")),

    htmltools::tags$div(class = "nav-section", "Output"),
    nav_item("nav_report",      "Report",         bsicons::bs_icon("file-earmark-text")),

    htmltools::tags$div(
      class = "sidebar-footer",
      htmltools::tags$span(class = "status-dot"),
      htmltools::tags$span(
        paste0("omicsCore ", utils::packageVersion("omicsCore"), " \u00B7 ready")
      )
    )
  )
}

# An actionLink styled as a sidebar nav row. We hand-roll it (rather
# than using a Shiny `tabsetPanel`) so the active state lives in a
# reactiveVal — letting us drive both the nav highlight and the hidden
# tabsetPanel from a single source of truth.
nav_item <- function(input_id, label, icon = NULL, active = FALSE) {
  cls <- c("nav-item", if (isTRUE(active)) "active")
  shiny::actionLink(
    inputId = input_id,
    label = htmltools::tagList(
      if (!is.null(icon)) icon,
      htmltools::tags$span(label)
    ),
    class = paste(cls, collapse = " ")
  )
}

# ---- header contents --------------------------------------------------

app_header <- function() {
  # The header shows just the project status. Earlier mockups had
  # Export / Run / Settings buttons here, but they had no server
  # bindings — the Re-run button on each analysis view and the
  # downloads in the Report view cover those affordances already.
  shiny::uiOutput("project_picker")
}
