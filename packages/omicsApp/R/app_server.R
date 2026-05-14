#' Top-level server for omicsApp
#'
#' Phase 2 slice 2A scaffold. Wires the 7 hand-rolled sidebar
#' `actionLink`s to a `reactiveVal("current_view")` that drives both:
#'   * the hidden [shiny::tabsetPanel()] in `app_ui()`, via
#'     [shiny::updateTabsetPanel()], and
#'   * the `.active` class on the matching nav item, via
#'     [shinyjs::addClass()] / [shinyjs::removeClass()].
#'
#' All view modules are mounted on session start so their UI is
#' rendered once. The 2A stubs have empty servers, but the call sites
#' below stay in place so 2B-2F can fill them in with real reactive
#' wiring without touching this file.
#'
#' @param input,output,session Shiny session triplet.
#'
#' @keywords internal
#' @noRd
app_server <- function(input, output, session) {
  # ---- view router ----------------------------------------------------
  views <- c("project", "import", "qc", "diff",
             "enrich", "integration", "report")
  current_view <- shiny::reactiveVal("project")

  set_view <- function(v) {
    if (!v %in% views) return(invisible(NULL))
    if (identical(current_view(), v)) return(invisible(NULL))
    current_view(v)
  }

  shiny::observeEvent(input$nav_project,     set_view("project"))
  shiny::observeEvent(input$nav_import,      set_view("import"))
  shiny::observeEvent(input$nav_qc,          set_view("qc"))
  shiny::observeEvent(input$nav_diff,        set_view("diff"))
  shiny::observeEvent(input$nav_enrich,      set_view("enrich"))
  shiny::observeEvent(input$nav_integration, set_view("integration"))
  shiny::observeEvent(input$nav_report,      set_view("report"))

  shiny::observeEvent(current_view(), {
    v <- current_view()
    shiny::updateTabsetPanel(session, "view", selected = v)
    for (other in setdiff(views, v)) {
      shinyjs::removeClass(id = paste0("nav_", other), class = "active")
    }
    shinyjs::addClass(id = paste0("nav_", v), class = "active")
  }, ignoreInit = TRUE)

  # ---- view modules ---------------------------------------------------
  project_view_server("project")
  import_view_server("import")
  qc_view_server("qc")
  diff_view_server("diff")
  enrich_view_server("enrich")
  integration_view_server("integration")
  report_view_server("report")
}
