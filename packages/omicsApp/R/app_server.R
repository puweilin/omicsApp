#' Top-level server for omicsApp
#'
#' Owns the shared app state:
#'   * `current_view` (slice 2A) — which view tab is active.
#'   * `current_project` (slice 3B) — the live `omics_project`, or
#'     `NULL` if the user hasn't imported anything yet. Every view
#'     reads this and falls back to fixture data when it is `NULL`.
#'
#' All view modules are mounted on session start. Modules that need
#' the project state receive it as a reactive argument (slice 3B
#' wires this for project + import; later slices add qc, diff,
#' enrich, integration, report).
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

  # ---- shared project state ------------------------------------------
  # NULL = no user data yet, views fall back to example_*().
  # Non-NULL = an omics_project built from one or more confirmed
  # omics_input objects emitted by the Import view.
  current_project <- shiny::reactiveVal(NULL)

  # ---- view modules ---------------------------------------------------
  project_view_server("project", current_project = current_project)
  imported_input <- import_view_server("import")
  qc_view_server("qc", current_project = current_project)
  diff_bundle <- diff_view_server("diff", current_project = current_project)
  enrich_bundle <- enrich_view_server("enrich", diff_bundle = diff_bundle)
  integration_view_server("integration",
                          current_project = current_project,
                          diff_bundle     = diff_bundle)
  report_view_server("report")

  # Every time the Import view confirms a fresh omics_input, fold it
  # into the project under a tag derived from its omics_type. Repeated
  # imports of the same omics_type replace the existing layer; new
  # omics_types extend the project. The Project view re-renders
  # automatically because it observes `current_project`.
  shiny::observe({
    inp <- imported_input()
    if (is.null(inp)) return()
    proj <- current_project()
    tag <- inp$omics_type %||% "experiment"
    if (is.null(proj)) {
      proj <- omicsCore::omics_project(
        name        = "User project",
        experiments = stats::setNames(list(inp), tag)
      )
    } else {
      # `add_experiment()` rejects duplicate tags, so drop any
      # existing layer of this omics_type first.
      if (tag %in% omicsCore::experiment_tags(proj)) {
        proj <- omicsCore::remove_experiment(proj, tag)
      }
      proj <- omicsCore::add_experiment(proj, name = tag, input = inp)
    }
    current_project(proj)
  })
}

# ---- internal helpers ------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a
