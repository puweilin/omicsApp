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
    sidebar = bslib::sidebar(
      id = "main_sidebar",
      open = "always",
      width = 232,
      bg = "#161A26",
      padding = c(14, 12, 14, 12),
      gap = 2,
      app_sidebar_content()
    ),
    htmltools::tags$div(
      class = "app-main",
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
  )
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
  htmltools::tagList(
    htmltools::tags$div(
      class = "project-picker",
      bsicons::bs_icon("collection"),
      htmltools::tags$div(
        htmltools::tags$span(class = "label-sm", "Project"),
        htmltools::tags$span("(no project loaded)")
      ),
      bsicons::bs_icon("chevron-down")
    ),
    htmltools::tags$div(
      class = "omics-tabs",
      htmltools::tags$button(
        class = "omics-tab active", `data-omics` = "proteomics",
        htmltools::tags$span(class = "dot"), "Proteomics"
      ),
      htmltools::tags$button(
        class = "omics-tab", `data-omics` = "rnaseq",
        htmltools::tags$span(class = "dot"), "RNA-seq"
      ),
      htmltools::tags$button(
        class = "omics-tab", `data-omics` = "integration",
        htmltools::tags$span(class = "dot"), "Integrated"
      )
    ),
    htmltools::tags$div(
      class = "app-header-actions",
      shiny::actionButton(
        "header_export", "Export",
        icon = bsicons::bs_icon("download"),
        class = "btn btn-outline-secondary"
      ),
      shiny::actionButton(
        "header_run", "Run analysis",
        icon = bsicons::bs_icon("play-fill"),
        class = "btn btn-primary"
      ),
      shiny::actionButton(
        "header_settings", label = NULL,
        icon = bsicons::bs_icon("gear"),
        class = "btn-icon",
        title = "Settings"
      )
    )
  )
}
