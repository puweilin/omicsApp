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

  # Bumped whenever an import replaces an existing omics layer. The
  # analysis views own their own bundles, so `app_server()` cannot clear
  # them directly — this is the signal they watch. Without it a bundle
  # computed on the replaced data is re-attached by the observer below
  # and ends up in both the report and the autosave, presented as if it
  # belonged to the data that just arrived.
  layer_generation <- shiny::reactiveVal(0L)

  # ---- header: project picker ----------------------------------------
  output$project_picker <- shiny::renderUI({
    proj <- current_project()
    if (is.null(proj)) {
      name <- "(no project loaded)"
      icon <- bsicons::bs_icon("collection")
    } else {
      name <- proj$name %||% "Unnamed project"
      icon <- bsicons::bs_icon("collection-fill")
    }
    htmltools::tags$div(
      class = "project-picker",
      icon,
      htmltools::tags$div(
        htmltools::tags$span(class = "label-sm", "Project"),
        htmltools::tags$span(name)
      )
    )
  })

  # ---- view modules ---------------------------------------------------
  # The Project view's per-layer "View" link lands on QC with that layer
  # selected. The tag travels through a reactiveVal rather than through
  # set_view() because the two are separate concerns: which view is on
  # screen, and which layer that view should be showing.
  requested_layer <- shiny::reactiveVal(NULL)
  project_view_server("project", current_project = current_project,
                      on_view_layer = function(tag) {
                        requested_layer(tag)
                        set_view("qc")
                      })
  imported_input <- import_view_server("import",
                        current_project = current_project)
  qc_bundle <- qc_view_server("qc", current_project = current_project,
                              invalidate = layer_generation,
                              requested_layer = requested_layer)
  diff_bundle <- diff_view_server("diff", current_project = current_project,
                                  invalidate = layer_generation)
  enrich_bundle <- enrich_view_server("enrich", diff_bundle = diff_bundle,
                                      invalidate = layer_generation)
  integration_bundle <- integration_view_server("integration",
                          current_project = current_project,
                          diff_bundle     = diff_bundle,
                          invalidate      = layer_generation)
  report_view_server("report", current_project = current_project)

  # Attach analysis bundles to the project whenever any bundle changes,
  # so that export_report() can walk names(project$bundles). We assign
  # back through `current_project(proj)` to propagate the change — a
  # direct `proj$bundles <- b` mutation would NOT invalidate the
  # reactiveVal, leaving the report view stuck on stale bundle status.
  shiny::observe({
    proj <- current_project()
    if (is.null(proj)) return()
    b <- list()
    qb <- qc_bundle()
    if (!is.null(qb)) b$qc <- qb
    db <- diff_bundle()
    if (!is.null(db)) b$diff <- db
    eb <- enrich_bundle()
    if (!is.null(eb)) b$enrich <- eb
    ib <- integration_bundle()
    if (!is.null(ib)) b$integration <- ib
    if (identical(proj$bundles, b)) return()
    proj$bundles <- b
    current_project(proj)
  })

  # Every time the Import view confirms a fresh omics_input, fold it
  # into the project under a tag derived from its omics_type. Repeated
  # imports of the same omics_type replace the existing layer; new
  # omics_types extend the project. The Project view re-renders
  # automatically because it observes `current_project`.
  #
  # `priority = 10` ensures this fires before the bundle-attach observer
  # in the same flush cycle, so bundles always attach to the freshest
  # project rather than the old one.
  shiny::observeEvent(imported_input(), {
    inp <- imported_input()
    if (is.null(inp)) return()
    proj <- current_project()
    tag <- inp$omics_type %||% "experiment"
    if (is.null(proj)) {
      proj <- omicsCore::omics_project(
        name        = "User project",
        experiments = stats::setNames(list(inp), tag)
      )
      proj$bundles <- list()
    } else {
      # `add_experiment()` rejects duplicate tags, so drop any
      # existing layer of this omics_type first. Bundles computed on
      # the replaced layer are also dropped since their input is now
      # gone — leaving them around shows stale results in the report.
      #
      # Clearing `proj$bundles` is not enough on its own: the analysis
      # views still hold what they computed and the bundle-attach
      # observer would put it straight back. Bumping the generation is
      # what actually tells them to let go.
      if (tag %in% omicsCore::experiment_tags(proj)) {
        proj <- omicsCore::remove_experiment(proj, tag)
        if (!is.null(proj$bundles)) proj$bundles <- list()
        layer_generation(shiny::isolate(layer_generation()) + 1L)
      }
      proj <- omicsCore::add_experiment(proj, name = tag, input = inp)
      if (is.null(proj$bundles)) proj$bundles <- list()
    }
    current_project(proj)
  }, priority = 10L)

  # ---- autosave -------------------------------------------------------
  # This is what makes an aggressive container-recycling policy safe: a
  # reaped session costs the user a page reload, not an afternoon of
  # analysis. Failures are swallowed inside `store_autosave()` — a full
  # disk must not interrupt the analysis in progress.
  #
  # Wired to the project *state* rather than to each analysis-completion
  # *event*: `current_project` is updated by the bundle-attach observer
  # above, so an observer firing on the same flush as a finished bundle
  # would read the project from before that bundle was folded in, and
  # persist a snapshot missing the very result that triggered it.
  wire_autosave(current_project)
}

# ---- internal helpers ------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a
