#' Project view module
#'
#' Top of the workflow funnel: tells the user which project they're
#' looking at, how many samples and features it covers, and what
#' has been run on it.
#'
#' Slice 3B: the view now reacts to a shared `current_project`
#' reactiveVal. When it is `NULL` (no user data yet) we render the
#' built-in `example_project()` and flag it with a "demo project"
#' label; when it is non-NULL the view re-renders against the
#' real project the Import view has confirmed.
#'
#' Reference markup: `omicsApp/mockup/index.html:462-562`.
#'
#' @param id Module namespace id.
#' @keywords internal
#' @noRd
project_view_ui <- function(id) {
  ns <- shiny::NS(id)
  htmltools::tagList(
    shiny::uiOutput(ns("header")),
    # First thing on the first page: what this app does, the order to do
    # it in, and a project to try it on. Replaced by a progress checklist
    # once a project is loaded.
    shiny::uiOutput(ns("guide")),
    shiny::uiOutput(ns("stats")),
    shiny::uiOutput(ns("body")),
    shiny::uiOutput(ns("storage"))
  )
}

#' @rdname project_view_ui
#' @param current_project Reactive (or reactiveVal) yielding the
#'   live `omics_project` or `NULL`.
#' @param on_view_layer Called with an experiment tag when its "View"
#'   link is clicked. The default does nothing, which keeps the module
#'   usable on its own; the app supplies a callback that switches to the
#'   QC view and asks it for that layer.
#' @keywords internal
#' @noRd
project_view_server <- function(id, current_project = shiny::reactiveVal(NULL),
                                on_view_layer = function(tag) NULL,
                                navigate = NULL) {
  shiny::moduleServer(id, function(input, output, session) {

    # Render against the live project when present, otherwise the
    # built-in demo. We keep a single source of truth (resolved())
    # so all three uiOutputs render in lockstep.
    resolved <- shiny::reactive({
      proj <- current_project()
      if (is.null(proj)) {
        list(project = example_project(), is_demo = TRUE)
      } else {
        list(project = proj, is_demo = FALSE)
      }
    })

    output$header <- shiny::renderUI({
      r <- resolved()
      view_header(
        title    = "Project overview",
        subtitle = htmltools::tagList(
          r$project$name,
          htmltools::HTML(" &middot; "),
          htmltools::tags$span(
            class = "muted",
            if (r$is_demo) "demo data (built-in)"
            else sprintf("%d layer%s loaded",
                         length(r$project$experiments),
                         if (length(r$project$experiments) == 1L) "" else "s")
          )
        )
      )
    })

    output$stats <- shiny::renderUI({
      r <- resolved()
      experiments <- r$project$experiments
      n_experiments <- length(experiments)

      if (n_experiments == 0L) {
        return(htmltools::tags$div(
          class = "stat-grid",
          stat_card(
            label  = "Experiments",
            value  = 0L,
            trend  = "import a file to populate",
            accent = "brand"
          )
        ))
      }

      n_samples_per_exp  <- vapply(experiments,
                                   function(x) ncol(x$expr_mat), integer(1))
      n_features_per_exp <- vapply(experiments,
                                   function(x) nrow(x$expr_mat), integer(1))
      total_features <- sum(n_features_per_exp)

      htmltools::tags$div(
        class = "stat-grid",
        stat_card(
          label  = "Experiments",
          value  = n_experiments,
          trend  = paste(vapply(experiments,
                                project_omics_label,
                                character(1)),
                         collapse = " + "),
          accent = "brand"
        ),
        stat_card(
          label = "Samples",
          value = sum(n_samples_per_exp),
          trend = paste(sprintf("%s: %d",
                                vapply(experiments, project_omics_label,
                                       character(1)),
                                n_samples_per_exp),
                        collapse = " \u00B7 "),
          mono  = TRUE
        ),
        stat_card(
          label = "Features (total)",
          value = format(total_features, big.mark = ","),
          trend = paste(sprintf("%s: %s",
                                vapply(experiments,
                                       project_omics_label,
                                       character(1)),
                                format(n_features_per_exp, big.mark = ",")),
                        collapse = " \u00B7 "),
          mono  = TRUE
        ),
        stat_card(
          label  = "Analyses",
          value  = if (r$is_demo) 0L else length(r$project$bundles),
          trend  = if (r$is_demo) "load a project to run them"
                   else if (length(r$project$bundles))
                     paste(names(r$project$bundles), collapse = " \u00B7 ")
                   else "none yet \u2014 start with Quality control",
          accent = if (!r$is_demo && length(r$project$bundles)) "ok" else "brand"
        )
      )
    })

    output$body <- shiny::renderUI({
      r <- resolved()
      htmltools::tags$div(
        class = "row-grid r-7-5",
        project_experiments_card(r$project$experiments, ns = session$ns),
        project_activity_card(r$project, is_demo = r$is_demo)
      )
    })

    # ---- getting started ---------------------------------------------
    output$guide <- shiny::renderUI({
      proj <- current_project()
      if (is.null(proj)) return(welcome_card(session$ns))
      workflow_card(proj, session$ns, can_navigate = is.function(navigate))
    })

    shiny::observeEvent(input$load_tutorial, {
      proj <- current_project()
      if (!is.null(proj) && length(proj$experiments)) {
        # Replacing a user's project is not something a tutorial button
        # does without asking.
        shiny::showModal(shiny::modalDialog(
          title = "Replace the current project with the example?",
          htmltools::tags$p("The example project replaces what is loaded now. ",
                            "Save your project first if you want to keep it."),
          footer = htmltools::tagList(
            shiny::modalButton("Cancel"),
            shiny::actionButton(session$ns("confirm_tutorial"), "Load example",
                                class = "btn btn-primary")),
          easyClose = TRUE))
        return()
      }
      current_project(tutorial_project())
      shiny::showNotification(
        "Example project loaded: proteomics + RNA-seq, Control vs TreatA / TreatB.",
        type = "message")
    })

    shiny::observeEvent(input$confirm_tutorial, {
      shiny::removeModal()
      current_project(tutorial_project())
    })

    shiny::observeEvent(input$go_next, {
      proj <- current_project()
      if (is.null(proj) || !is.function(navigate)) return()
      navigate(WORKFLOW_STEPS$id[workflow_progress(proj)$next_step])
    })

    # The welcome card's button and the checklist's Import row share this id.
    shiny::observeEvent(input$go_import, {
      if (is.function(navigate)) navigate("import")
    })

    for (step in setdiff(WORKFLOW_STEPS$id, "import")) local({
      target <- step
      shiny::observeEvent(input[[paste0("go_", target)]], {
        if (is.function(navigate)) navigate(target)
      }, ignoreInit = TRUE)
    })

    # One observer per experiment row, registered the first time a row
    # at that position exists. Registering inside an observer rather
    # than up front is what lets the number of rows change;
    # `registered_rows` stops a second registration when the project is
    # replaced, which would otherwise fire the callback once per
    # duplicate.
    pending_drop <- shiny::reactiveVal(NULL)
    drop_layer_confirm <- function(session, tag) {
      pending_drop(tag)
      shiny::showModal(shiny::modalDialog(
        title = sprintf("Remove layer '%s'?", tag),
        htmltools::tags$p(
          "The imported matrix, its metadata and every result computed ",
          "on it go with it. Other layers in this project are untouched."
        ),
        htmltools::tags$p(
          class = "muted",
          "Saved projects on disk are not affected until you save again."
        ),
        footer = htmltools::tagList(
          shiny::modalButton("Cancel"),
          shiny::actionButton(session$ns("confirm_drop_layer"), "Remove",
                              class = "btn btn-danger")
        ),
        easyClose = TRUE
      ))
    }

    registered_rows <- new.env(parent = emptyenv())
    shiny::observe({
      n <- length(resolved()$project$experiments)
      for (i in seq_len(n)) {
        key <- paste0("view_layer_", i)
        if (!is.null(registered_rows[[key]])) next
        registered_rows[[key]] <- TRUE
        local({
          idx <- i
          shiny::observeEvent(input[[paste0("view_layer_", idx)]], {
            # Read the tag at click time: the project may have been
            # replaced since the link was drawn, and the row now
            # belongs to a different layer.
            tags_now <- names(resolved()$project$experiments)
            if (idx <= length(tags_now)) on_view_layer(tags_now[[idx]])
          }, ignoreInit = TRUE)

          shiny::observeEvent(input[[paste0("drop_layer_", idx)]], {
            proj <- current_project()
            if (is.null(proj)) return()
            tags_now <- names(proj$experiments)
            if (idx > length(tags_now)) return()
            drop_layer_confirm(session, tags_now[[idx]])
          }, ignoreInit = TRUE)
        })
      }
    })

    # Confirmed, unlike the View link next to it: removing a layer throws
    # away an import and every result computed on it, and the two links
    # are one word apart.
    shiny::observeEvent(input$confirm_drop_layer, {
      proj <- current_project()
      tag <- shiny::isolate(pending_drop())
      shiny::removeModal()
      if (is.null(proj) || is.null(tag) || !tag %in% names(proj$experiments)) {
        return()
      }
      proj$experiments[[tag]] <- NULL
      # A link naming a layer that is gone would pair samples to nothing.
      if (!is.null(proj$sample_link) && nrow(proj$sample_link) > 0L) {
        proj$sample_link <- proj$sample_link[proj$sample_link$tag != tag, ,
                                             drop = FALSE]
      }
      current_project(proj)
      pending_drop(NULL)
      shiny::showNotification(sprintf("Removed layer '%s'.", tag),
                              type = "message")
    })

    # ---- saved-project store -----------------------------------------
    # `store_tick` is bumped after every mutation so the card re-reads
    # the directory. Reading it inside the reactives below is what makes
    # them invalidate; the value itself is never used.
    store_tick <- shiny::reactiveVal(0L)
    bump_store <- function() {
      store_tick(shiny::isolate(store_tick()) + 1L)
    }

    saved_projects <- shiny::reactive({
      store_tick()
      list_saved_projects()
    })

    autosave_stamp <- shiny::reactive({
      store_tick()
      autosave_mtime()
    })

    output$storage <- shiny::renderUI({
      saved <- saved_projects()
      stamp <- autosave_stamp()
      have_project <- !is.null(current_project())
      choices <- if (nrow(saved) == 0L) character(0) else saved$slug

      htmltools::tags$div(
        class = "row-grid r-7-5",
        bslib::card(
          bslib::card_header(
            htmltools::tags$h3(class = "card-title", "My projects"),
            htmltools::tags$span(class = "card-sub", usage_label())
          ),
          bslib::card_body(
            if (length(choices) == 0L) {
              htmltools::tags$div(
                class = "muted", style = "font-size:13px;padding:4px 0 10px",
                "No saved projects yet. Save the current one to keep it ",
                "across sessions."
              )
            } else {
              htmltools::tagList(
                shiny::selectInput(session$ns("saved_pick"), label = NULL,
                                   choices = choices, selectize = FALSE,
                                   size = min(6L, length(choices))),
                htmltools::tags$div(
                  style = "display:flex;gap:8px",
                  shiny::actionButton(session$ns("open_project"), "Open",
                                      class = "btn btn-sm btn-primary"),
                  shiny::actionButton(session$ns("delete_project"), "Delete",
                                      class = "btn btn-sm btn-outline-danger")
                ),
                saved_projects_table(saved)
              )
            }
          )
        ),
        bslib::card(
          bslib::card_header(
            htmltools::tags$h3(class = "card-title", "Save / restore")
          ),
          bslib::card_body(
            shiny::textInput(session$ns("save_name"), "Project name",
                             placeholder = "e.g. cheek_G2_vs_G1"),
            shiny::checkboxInput(session$ns("save_overwrite"),
                                 "Overwrite if it exists", value = FALSE),
            shiny::actionButton(session$ns("save_project"), "Save as\u2026",
                                class = "btn btn-sm btn-primary"),
            if (!have_project) {
              htmltools::tags$div(
                class = "muted", style = "font-size:12px;padding-top:8px",
                "The built-in demo cannot be saved \u2014 import a file or load the example project first."
              )
            },
            if (!is.null(stamp)) {
              htmltools::tagList(
                htmltools::tags$hr(),
                htmltools::tags$div(
                  class = "muted", style = "font-size:12px;padding-bottom:6px",
                  sprintf("Autosave from %s",
                          format(stamp, "%Y-%m-%d %H:%M"))
                ),
                shiny::actionButton(session$ns("restore_autosave"),
                                    "Restore last session",
                                    class = "btn btn-sm btn-outline-primary")
              )
            }
          )
        )
      )
    })

    shiny::observeEvent(input$save_project, {
      proj <- current_project()
      if (is.null(proj)) {
        shiny::showNotification(
          "Nothing to save yet \u2014 import a file first.", type = "warning")
        return()
      }
      res <- store_save_project(
        proj,
        slug      = project_slug(input$save_name %||% ""),
        overwrite = isTRUE(input$save_overwrite)
      )
      shiny::showNotification(res$message,
                              type = if (isTRUE(res$ok)) "message" else "error",
                              duration = if (isTRUE(res$ok)) 5 else NULL)
      if (isTRUE(res$ok)) bump_store()
    })

    shiny::observeEvent(input$open_project, {
      res <- store_load_project(input$saved_pick %||% NA_character_)
      shiny::showNotification(res$message,
                              type = if (isTRUE(res$ok)) "message" else "error",
                              duration = if (isTRUE(res$ok)) 5 else NULL)
      if (isTRUE(res$ok)) current_project(res$project)
    })

    # Confirmed, as removing a layer is: a deleted project file does not
    # come back.
    shiny::observeEvent(input$delete_project, {
      slug <- input$saved_pick %||% NA_character_
      if (is.na(slug) || !nzchar(slug)) return()
      shiny::showModal(shiny::modalDialog(
        title = sprintf("Delete saved project '%s'?", slug),
        htmltools::tags$p("The saved file is removed from the server. ",
                          "This cannot be undone."),
        footer = htmltools::tagList(
          shiny::modalButton("Cancel"),
          shiny::actionButton(session$ns("confirm_delete_project"), "Delete",
                              class = "btn btn-danger")),
        easyClose = TRUE))
    })

    shiny::observeEvent(input$confirm_delete_project, {
      shiny::removeModal()
      res <- store_delete_project(input$saved_pick %||% NA_character_)
      shiny::showNotification(res$message,
                              type = if (isTRUE(res$ok)) "message" else "error",
                              duration = if (isTRUE(res$ok)) 5 else NULL)
      if (isTRUE(res$ok)) bump_store()
    })

    shiny::observeEvent(input$restore_autosave, {
      proj <- store_read_autosave()
      if (is.null(proj)) {
        shiny::showNotification("No readable autosave found.", type = "error")
        return()
      }
      current_project(proj)
      shiny::showNotification("Restored the last autosaved session.",
                              type = "message")
    })

    # Restore on arrival rather than on a click. The autosave is the
    # user's own last state, and asking them to ask for it every login
    # meant landing on the built-in demo each morning -- a page that is
    # useful once and then noise.
    #
    # Deliberately narrow: only when nothing is loaded, and only once
    # per session, so it cannot overwrite an import that happened first
    # or fire again after the user clears the project.
    #
    # The button stays, because it is now the way back to the autosave
    # after the project has been changed in-session.
    shiny::observe({
      if (!is.null(current_project())) return()
      proj <- tryCatch(store_read_autosave(), error = function(e) NULL)
      if (is.null(proj)) return()
      current_project(proj)
      shiny::showNotification("Restored your last session.", type = "message")
    }) |> shiny::bindEvent(TRUE, once = TRUE, ignoreInit = FALSE)
  })
}

# ---- internal helpers ------------------------------------------------

# The order the views are meant to be used in, and what each is for, in
# one place so the welcome card and the checklist cannot disagree.
WORKFLOW_STEPS <- data.frame(
  id = c("import", "qc", "diff", "enrich", "integration", "report"),
  label = c("Import", "Quality control", "Differential", "Enrichment",
            "Integration", "Report"),
  desc = c(
    "Upload a workbook or CSV per omics layer; confirm the detected matrix, metadata and features.",
    "Look for missing values, outlier samples and batch structure before testing anything.",
    "Compare groups: one control against one or more treatments, with optional covariates.",
    "Ask which pathways the differential hits fall in (ORA / GSEA, MSigDB).",
    "Put two layers side by side: do RNA and protein change together?",
    "Download the report, the result tables and an R script that reproduces them."
  ),
  bundle = c(NA, "qc", "diff", "enrich", "integration", NA),
  tip = c(
    "The example is already imported: two layers, 12 samples each, sample ids differ but a donor column pairs them.",
    "Colour the PCA by 'group' \u2014 TreatA and TreatB separate from Control along PC1/PC2.",
    "Leave Control as the reference and keep both treatments selected, then press Run analysis. Switch between the two comparisons with 'Showing'.",
    "With 'TreatA vs Control' shown, ORA on Hallmark should find INFLAMMATORY_RESPONSE at the top.",
    "Concordance repeats the contrast on RNA-seq. Try 'Sample-level correlation' too \u2014 the donor column pairs the samples.",
    "Everything run so far is in the report and the script."
  ),
  stringsAsFactors = FALSE
)

is_tutorial_project <- function(proj) {
  isTRUE(startsWith(proj$name %||% "", "Tutorial"))
}

welcome_card <- function(ns) {
  steps <- WORKFLOW_STEPS
  bslib::card(
    class = "welcome-card",
    bslib::card_header(
      htmltools::tags$h3(class = "card-title", "Welcome to omicsApp"),
      htmltools::tags$span(class = "card-sub",
                           "proteomics + transcriptomics, from upload to report")
    ),
    bslib::card_body(
      htmltools::tags$p(
        "omicsApp takes one or more omics layers from the same study through ",
        "quality control, differential analysis, pathway enrichment and ",
        "RNA\u2013protein integration. The pages show built-in demo data ",
        "until you load a project."
      ),
      htmltools::tags$ol(
        class = "workflow-list",
        lapply(seq_len(nrow(steps)), function(i) {
          htmltools::tags$li(
            htmltools::tags$strong(steps$label[i]), " \u2014 ", steps$desc[i])
        })
      ),
      htmltools::tags$div(
        class = "welcome-actions",
        shiny::actionButton(ns("load_tutorial"), "Try the example project",
                            icon = shiny::icon("graduation-cap"),
                            class = "btn btn-primary"),
        shiny::actionButton(ns("go_import"), "Import my own data",
                            icon = shiny::icon("upload"),
                            class = "btn btn-outline-primary")
      ),
      htmltools::tags$div(
        class = "muted", style = "font-size:12px;margin-top:8px",
        "The example is small and synthetic (two layers, 252 genes, a control ",
        "and two treatments), and every step of the workflow finds something in it. ",
        "Templates for your own files are on the Import page."
      )
    )
  )
}

# Which steps a project has been through, and which comes next: the
# first not done, skipping Integration when there is only one layer.
workflow_progress <- function(proj) {
  steps <- WORKFLOW_STEPS
  have <- names(proj$bundles %||% list())
  n_layers <- length(proj$experiments)
  done <- vapply(seq_len(nrow(steps)), function(i) {
    if (steps$id[i] == "import") return(n_layers > 0L)
    if (steps$id[i] == "report") return(FALSE)
    !is.na(steps$bundle[i]) && steps$bundle[i] %in% have
  }, logical(1))
  skip <- steps$id == "integration" & n_layers < 2L
  nxt <- which(!done & !skip)[1L]
  if (is.na(nxt)) nxt <- nrow(steps)
  list(done = done, skip = skip, next_step = nxt)
}

workflow_card <- function(proj, ns, can_navigate = TRUE) {
  steps <- WORKFLOW_STEPS
  prog <- workflow_progress(proj)
  done <- prog$done
  skip <- prog$skip
  nxt <- prog$next_step
  tutorial <- is_tutorial_project(proj)

  items <- lapply(seq_len(nrow(steps)), function(i) {
    state <- if (done[i]) "done" else if (i == nxt) "active" else "pending"
    desc <- if (skip[i]) "needs a second layer" else NULL
    htmltools::tags$div(
      class = "workflow-step",
      step_item(i, steps$label[i], desc, state = state),
      if (can_navigate) {
        htmltools::tags$button(
          id = ns(paste0("go_", steps$id[i])), type = "button",
          class = "btn btn-sm btn-link action-button", "Open")
      }
    )
  })
  bslib::card(
    class = "workflow-card",
    bslib::card_header(
      htmltools::tags$h3(class = "card-title",
                         if (tutorial) "Tutorial" else "Workflow"),
      htmltools::tags$span(class = "card-sub",
                           sprintf("next: %s", steps$label[nxt]))
    ),
    bslib::card_body(
      htmltools::tags$div(class = "workflow-steps", items),
      htmltools::tags$div(
        class = "workflow-next",
        htmltools::tags$div(
          htmltools::tags$strong(sprintf("Step %d \u00B7 %s", nxt, steps$label[nxt])),
          htmltools::tags$div(class = "muted",
                              if (tutorial) steps$tip[nxt] else steps$desc[nxt])
        ),
        if (can_navigate) {
          shiny::actionButton(ns("go_next"),
                              sprintf("Go to %s \u2192", steps$label[nxt]),
                              class = "btn btn-primary btn-sm")
        }
      )
    )
  )
}

# Compact listing under the project picker: size and last-modified for
# each saved `.omp`, so a user can tell two similarly-named projects
# apart before opening one.
saved_projects_table <- function(saved) {
  if (!is.data.frame(saved) || nrow(saved) == 0L) return(NULL)
  htmltools::tags$table(
    class = "tbl",
    style = "margin-top:12px",
    htmltools::tags$thead(
      htmltools::tags$tr(
        htmltools::tags$th("Project"),
        htmltools::tags$th(class = "num", "Size"),
        htmltools::tags$th("Modified")
      )
    ),
    htmltools::tags$tbody(
      lapply(seq_len(nrow(saved)), function(i) {
        htmltools::tags$tr(
          htmltools::tags$td(
            htmltools::tags$span(class = "text-mono", saved$slug[[i]])),
          htmltools::tags$td(class = "num",
                             sprintf("%.1f MB", saved$size_mb[[i]])),
          htmltools::tags$td(
            class = "muted",
            format(saved$modified[[i]], "%Y-%m-%d %H:%M"))
        )
      })
    )
  )
}

# Human-readable display label for an omics_input layer.
project_omics_label <- function(x) {
  switch(x$omics_type,
         proteomics = "Proteomics",
         rnaseq     = "RNA-seq",
         x$omics_type)
}

#' @param ns The module's namespace function. The "View" cell is an
#'   `actionLink` when it is supplied, and plain text when it is not --
#'   the card is also rendered in contexts with no session to click in.
#'   Links are keyed by row position, not by tag: a tag is user-supplied
#'   text and would make an unsafe input id.
#' @noRd
project_experiments_card <- function(experiments, ns = NULL) {
  bslib::card(
    bslib::card_header(
      htmltools::tags$h3(class = "card-title", "Experiments"),
      htmltools::tags$span(
        class = "card-sub",
        sprintf("%d layer%s loaded", length(experiments),
                if (length(experiments) == 1L) "" else "s")
      )
    ),
    bslib::card_body(
      if (length(experiments) == 0L) {
        htmltools::tags$div(
          class = "muted",
          style = "font-size:13px;padding:8px 0",
          "No experiments yet. Use the Import view to upload a file."
        )
      } else {
        htmltools::tags$table(
          class = "tbl",
          htmltools::tags$thead(
            htmltools::tags$tr(
              htmltools::tags$th("Tag"),
              htmltools::tags$th("Omics"),
              htmltools::tags$th(class = "num", "Samples"),
              htmltools::tags$th(class = "num", "Features"),
              htmltools::tags$th("Status"),
              htmltools::tags$th("")
            )
          ),
          htmltools::tags$tbody(
            lapply(seq_along(experiments), function(i) {
              tag <- names(experiments)[i]
              exp <- experiments[[i]]
              htmltools::tags$tr(
                htmltools::tags$td(htmltools::tags$span(class = "text-mono", tag)),
                htmltools::tags$td(project_omics_label(exp)),
                htmltools::tags$td(class = "num", ncol(exp$expr_mat)),
                htmltools::tags$td(class = "num",
                                   format(nrow(exp$expr_mat), big.mark = ",")),
                htmltools::tags$td(pill("ready", kind = "ok")),
                htmltools::tags$td(
                  if (is.null(ns)) {
                    htmltools::tags$span(class = "muted", "View \u2192")
                  } else {
                    htmltools::tagList(
                      shiny::actionLink(ns(paste0("view_layer_", i)),
                                        "View \u2192"),
                      # Per layer, because a project is usually only
                      # wrong in one of them: re-importing a mistaken
                      # RNA-seq layer should not cost the proteomics
                      # work sitting next to it.
                      htmltools::tags$span(class = "muted", " \u00b7 "),
                      shiny::actionLink(ns(paste0("drop_layer_", i)),
                                        "Remove",
                                        class = "text-danger")
                    )
                  }
                )
              )
            })
          )
        )
      }
    )
  )
}

# Activity card. Static bullets for the demo project; a real activity
# log would need persistence (out of scope for Phase 3). For a user
# project we show a single "imported N experiments" bullet instead.
project_activity_card <- function(project, is_demo = TRUE) {
  bullet <- function(dot_var, title, meta) {
    htmltools::tags$div(
      style = "display:flex;gap:12px;padding:8px 0;border-bottom:1px dashed var(--border)",
      htmltools::tags$div(
        style = sprintf(
          "width:8px;height:8px;border-radius:50%%;background:var(%s);margin-top:6px;flex:none",
          dot_var
        )
      ),
      htmltools::tags$div(
        htmltools::tags$div(style = "font-size:13px;font-weight:500", title),
        htmltools::tags$div(class = "muted", style = "font-size:12px", meta)
      )
    )
  }
  bslib::card(
    bslib::card_header(htmltools::tags$h3(class = "card-title", "Recent activity")),
    bslib::card_body(
      style = "padding-top:8px",
      if (isTRUE(is_demo)) {
        # Nothing has happened yet, and the card says so. It used to list
        # a limma run "just now" that no view had performed -- the
        # Differential page then opened empty, contradicting it.
        htmltools::tags$div(
          class = "muted", style = "font-size:12px",
          "Nothing yet. Load the example project or import a file, and the ",
          "layers and analyses appear here."
        )
      } else {
        experiments <- project$experiments
        if (length(experiments) == 0L) {
          htmltools::tags$div(
            class = "muted",
            style = "font-size:12px",
            "Nothing yet."
          )
        } else {
          htmltools::tagList(
            lapply(names(experiments), function(tag) {
              exp <- experiments[[tag]]
              bullet(
                if (exp$omics_type == "rnaseq") "--brand-500" else "--ok",
                sprintf("Imported %s experiment",
                        project_omics_label(exp)),
                sprintf("tag = %s \u00B7 %d samples \u00B7 %d features",
                        tag, ncol(exp$expr_mat), nrow(exp$expr_mat))
              )
            }),
            lapply(names(project$bundles %||% list()), function(nm) {
              b <- project$bundles[[nm]]
              what <- switch(nm, qc = "Quality control", diff = "Differential",
                             enrich = "Enrichment", integration = "Integration", nm)
              prm <- if (is.list(b)) b$params else NULL
              detail <- paste(c(prm$method, prm$comparison, prm$type,
                                prm$database),
                              collapse = " \u00B7 ")
              bullet("--accent-500", what,
                     gsub("_vs_", " vs ", if (nzchar(detail)) detail else "done"))
            })
          )
        }
      }
    )
  )
}
