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

  # How many times a project has been made current in this session: opened,
  # restored, or first imported. Separate from project_id, which is the
  # project's own and travels with the file.
  #
  # Comparing project_id alone cannot see the case that matters most --
  # opening an *earlier save of the same project*. Nothing about the identity
  # or the data has changed, so the views kept the results they held and the
  # write-back put them into the file that had just been opened, over the
  # results it was opened for.
  #
  # It stays out of anything saved to disk. Written into a result's record it
  # would make that result invalid the next time the file was opened, which
  # is the opposite of what any of this is for.
  activation <- shiny::reactiveVal(0L)
  bump_activation <- function() {
    activation(shiny::isolate(activation()) + 1L)
    layer_generation(shiny::isolate(layer_generation()) + 1L)
  }

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
                      },
                      on_activate = bump_activation)
  imported_input <- import_view_server("import",
                        current_project = current_project)
  # The QC view returns what the user applied, not what it previewed:
  # list(tag, bundle, params, base_version) or NULL.
  qc_applied <- qc_view_server("qc", current_project = current_project,
                               invalidate = layer_generation,
                               requested_layer = requested_layer,
                               activation = activation)
  # The bundle half, for the report and for `project$bundles`.
  qc_bundle <- shiny::reactive(qc_applied()$bundle)

  # Applying QC changes what every downstream analysis means, so it is
  # written into the project rather than kept beside it: set_preprocessing()
  # records the cleaned layer against the import version it came from, and
  # effective_input() is then the single answer to "which data is this
  # layer". Results computed on the previous version are retired here, in
  # the project itself -- clearing the views alone left a restored project's
  # saved differential in place, where the report and the export still read
  # it.
  shiny::observeEvent(qc_applied(), {
    x <- qc_applied()
    proj <- current_project()
    if (is.null(x) || is.null(proj)) return()
    proj <- omicsCore::set_preprocessing(
      proj, x$tag,
      cleaned = x$bundle$results$cleaned_input,
      params = x$params,
      base_input_version = x$base_version
    )
    proj <- omicsCore::retire_invalid(proj, reason = "QC applied")
    current_project(proj)
    layer_generation(shiny::isolate(layer_generation()) + 1L)
  }, priority = 15L)
  diff_view <- diff_view_server("diff", current_project = current_project,
                                invalidate = layer_generation,
                                activation = activation)
  diff_bundle <- diff_view$bundle
  # Enrichment reads the same thresholds the hit table is read at.
  # Without this it applied its own -- adjusted p at 0.05, no
  # fold-change bound -- and enriched a different set of genes from the
  # ones on screen, silently.
  enrich_bundle <- enrich_view_server("enrich", diff_bundle = diff_bundle,
                                      current_project = current_project,
                                      activation = activation,
                                      diff_thresholds = diff_view$thresholds,
                                      diff_layer = diff_view$layer,
                                      invalidate = layer_generation)
  integration_view <- integration_view_server("integration",
                          current_project = current_project,
                          diff_bundle     = diff_bundle,
                          invalidate      = layer_generation,
                          activation      = activation)
  integration_bundle <- integration_view$bundle
  # The integration's second-layer differential. It has no panel, so it
  # reaches the project through the registry rather than through a view
  # slot -- an integration that names an upstream the project does not
  # hold is judged invalid as soon as it arrives.
  integration_secondary <- integration_view$secondary
  report_view_server("report", current_project = current_project)

  # ---- which data the views' results belong to -------------------------
  # A layer's identity is its upload fingerprint, or its shape when it
  # was built without one. The import observer below bumps the
  # generation for the replacement it performs itself; this covers the
  # arrivals it never sees -- Open, Restore, a dropped layer, the first
  # import after a demo run -- where the project changes under the
  # views and what they hold was computed on data that is no longer in
  # it. Found by the acceptance test: a project restored into a fresh
  # session lost its saved analyses within one flush, and a project
  # opened while the views still held another project's results was
  # handed those results as its own.
  layer_identity <- function(proj) {
    if (is.null(proj) || !length(proj$experiments)) return(character(0))
    vapply(proj$experiments, function(e) {
      e$source_fingerprint %||% paste(dim(e$expr_mat), collapse = "x")
    }, character(1))
  }
  seen_identity <- shiny::reactiveVal(NULL)
  shiny::observeEvent(current_project(), {
    now <- layer_identity(current_project())
    before <- seen_identity()
    seen_identity(now)
    if (is.null(before)) {
      # First project of the session. Whatever the views hold is about
      # the demo, and must not be attached to the user's data.
      if (length(now)) layer_generation(shiny::isolate(layer_generation()) + 1L)
      return()
    }
    gone <- setdiff(names(before), names(now))
    common <- intersect(names(before), names(now))
    changed <- common[before[common] != now[common]]
    if (length(gone) || length(changed)) {
      layer_generation(shiny::isolate(layer_generation()) + 1L)
    }
  }, priority = 20L, ignoreNULL = FALSE, ignoreInit = TRUE)

  # Views publish their results into the project, so that
  # export_report() can walk names(project$bundles). We assign back
  # through `current_project(proj)` to propagate the change -- a direct
  # `proj$bundles <- b` mutation would NOT invalidate the reactiveVal,
  # leaving the report view stuck on stale bundle status.
  #
  # Starts from the project's own bundles rather than from an empty
  # list: a project can arrive already carrying results (Restore and
  # Open bring them from disk), and rebuilding the list from the views
  # alone threw those away the moment they were loaded. A view that
  # lets go of a result on its own -- a change of layer, not a change
  # of data -- still takes it out of the project; a view emptied by a
  # generation bump does not, because whoever bumped has already
  # decided what the project keeps.
  #
  # Priority -10: after every view has reacted to a generation bump, so
  # this never sees a result computed on a layer that has just been
  # replaced. That ordering used to rest on which reactiveVal the
  # import observer happened to assign first.
  shiny::observe({
    proj <- current_project()
    if (is.null(proj)) return()
    layer_generation()
    views <- list(
      qc          = qc_bundle(),
      diff        = diff_bundle(),
      enrich      = enrich_bundle(),
      integration = integration_bundle()
    )
    pid <- proj$metadata$project_id
    b <- proj$bundles %||% list()
    analyses <- proj$analyses %||% list()

    for (nm in names(views)) {
      v <- views[[nm]]
      if (is.null(v)) next
      pb <- v$provenance_binding
      # A result that cannot say which data it came from is not written into
      # the project. Being present used to be enough, and `%||% pid` then
      # filled the gap with whichever project happened to be open -- not a
      # fallback but the mislabelling this exists to prevent.
      if (is.null(pb) || !identical(pb$project_id, pid)) next
      # Nor is one whose data has already moved on. The view will be told to
      # let go by the generation bump; until then it must not be filed, or
      # the retirement below would remove it and the next flush would put it
      # straight back.
      if (!omicsCore::result_inputs_current(pb, proj)) next
      b[[nm]] <- v
      analyses[[pb$analysis_id]] <- v
    }

    # Registered but not shown: nothing displays the integration's second
    # differential, and the integration depends on it.
    sec <- integration_secondary()
    if (!is.null(sec)) {
      spb <- sec$provenance_binding
      if (!is.null(spb) && identical(spb$project_id, pid) &&
          omicsCore::result_inputs_current(spb, proj)) {
        analyses[[spb$analysis_id]] <- sec
      }
    }

    proj$bundles  <- b
    proj$analyses <- analyses
    # Removal is decided from the project's own record of what each result
    # was computed on, not from whether a view happens to be holding
    # something. The previous rule required a view to have held a non-empty
    # bundle before it would drop the project's, so a project restored from
    # disk with a saved differential and an empty view kept that
    # differential through every change of data -- and the report, the
    # export and the autosave all went on reading it.
    proj <- omicsCore::retire_invalid(proj)
    if (identical(proj, shiny::isolate(current_project()))) return()
    current_project(proj)
  }, priority = -10L)

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
  autosave_status <- shiny::reactiveVal(NULL)
  wire_autosave(current_project, status = autosave_status)

  # Checked once, at startup, rather than discovered at the first save.
  # A read-only mount is a deployment fault, and finding out about it
  # after an afternoon's work is the worst moment.
  store_health <- store_writable()

  output$autosave_status <- shiny::renderUI({
    if (!isTRUE(store_health$ok)) {
      return(htmltools::tags$div(
        class = "autosave-status is-error",
        bsicons::bs_icon("exclamation-triangle-fill"),
        htmltools::tags$span("Not saving \u2014 the project store is not writable"),
        htmltools::tags$span(class = "label-sm", store_health$message)))
    }
    st <- autosave_status()
    if (is.null(st)) {
      return(htmltools::tags$div(class = "autosave-status is-idle",
                                 htmltools::tags$span("No analyses to save yet")))
    }
    if (isTRUE(st$ok)) {
      htmltools::tags$div(
        class = "autosave-status is-ok",
        bsicons::bs_icon("cloud-check"),
        htmltools::tags$span(sprintf("Saved %s", format(st$time, "%H:%M:%S"))))
    } else {
      htmltools::tags$div(
        class = "autosave-status is-error",
        bsicons::bs_icon("exclamation-triangle-fill"),
        htmltools::tags$span("Not saved"),
        htmltools::tags$span(class = "label-sm", st$message))
    }
  })
}

# ---- internal helpers ------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a
