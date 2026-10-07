# The saved-project store in the project view: My projects, Save /
# restore, and restoring the last session.
#
# Split out of mod_project_view.R. The server half is a plain function
# called from inside project_view_server()'s moduleServer(), not a
# module, so the ids stay "project-save_project", "project-open_project".
project_saved_server <- function(input, output, session, current_project, set_project) {
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
    # The project takes the name it is saved under; it stayed "User
    # project" in the header, the report and the script.
    nm <- trimws(input$save_name %||% "")
    if (nzchar(nm) && !identical(proj$name, nm)) {
      proj$name <- nm
      current_project(proj)
    }
    res <- store_save_project(
      proj,
      slug      = project_slug(input$save_name %||% ""),
      overwrite = isTRUE(input$save_overwrite)
    )
    # The name as typed, not its file-safe slug.
    shiny::showNotification(if (isTRUE(res$ok) && nzchar(nm)) sprintf("Saved '%s'.", nm) else res$message,
                            type = if (isTRUE(res$ok)) "message" else "error",
                            duration = if (isTRUE(res$ok)) 5 else NULL)
    if (isTRUE(res$ok)) bump_store()
  })

  shiny::observeEvent(input$open_project, {
    res <- store_load_project(input$saved_pick %||% NA_character_)
    shiny::showNotification(res$message,
                            type = if (isTRUE(res$ok)) "message" else "error",
                            duration = if (isTRUE(res$ok)) 5 else NULL)
    if (isTRUE(res$ok)) set_project(res$project)
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

  restore_from <- function(path = NULL) {
    proj <- store_read_autosave(path = path)
    if (is.null(proj)) {
      shiny::showNotification("No readable autosave found.", type = "error")
      return()
    }
    set_project(proj)
    shiny::showNotification("Restored the last autosaved session.",
                            type = "message")
  }

  # One snapshot per browser session, so with two tabs open there are
  # two to choose from; the newest is not necessarily the one wanted.
  shiny::observeEvent(input$restore_autosave, {
    snaps <- list_autosaves()
    if (nrow(snaps) <= 1L) return(restore_from())
    labels <- sprintf("%s \u00B7 %s layer(s) \u00B7 saved %s",
                      ifelse(is.na(snaps$name), "(unnamed)", snaps$name),
                      ifelse(is.na(snaps$n_layers), "?", snaps$n_layers),
                      format(snaps$modified, "%Y-%m-%d %H:%M"))
    shiny::showModal(shiny::modalDialog(
      title = "Restore which session?",
      shiny::radioButtons(session$ns("autosave_pick"), label = NULL,
                          choices = stats::setNames(basename(snaps$path), labels),
                          selected = basename(snaps$path)[[1L]]),
      htmltools::tags$p(class = "muted",
                        "Each browser tab keeps its own snapshot. Restoring replaces what is loaded now."),
      footer = htmltools::tagList(
        shiny::modalButton("Cancel"),
        shiny::actionButton(session$ns("confirm_restore_autosave"), "Restore",
                            class = "btn btn-primary")),
      easyClose = TRUE))
  })

  shiny::observeEvent(input$confirm_restore_autosave, {
    shiny::removeModal()
    restore_from(input$autosave_pick)
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
  invisible()
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
