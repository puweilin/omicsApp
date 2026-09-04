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
                                on_view_layer = function(tag) NULL) {
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
            if (r$is_demo) "demo project (built-in)"
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
          value = unname(n_samples_per_exp[1L]),
          trend = sprintf("%d per omics layer",
                          unname(n_samples_per_exp[1L])),
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
          label  = "Analyses cached",
          value  = if (r$is_demo) 1L else 0L,
          trend  = if (r$is_demo) "1 diff (limma, G2 vs G1)"
                   else "run from the Differential view",
          accent = if (r$is_demo) "ok" else "brand"
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
                "The built-in demo cannot be saved \u2014 import a file first."
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
                              type = if (isTRUE(res$ok)) "message" else "error")
      if (isTRUE(res$ok)) bump_store()
    })

    shiny::observeEvent(input$open_project, {
      res <- store_load_project(input$saved_pick %||% NA_character_)
      shiny::showNotification(res$message,
                              type = if (isTRUE(res$ok)) "message" else "error")
      if (isTRUE(res$ok)) current_project(res$project)
    })

    shiny::observeEvent(input$delete_project, {
      res <- store_delete_project(input$saved_pick %||% NA_character_)
      shiny::showNotification(res$message,
                              type = if (isTRUE(res$ok)) "message" else "error")
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
        htmltools::tagList(
          bullet("--brand-500", "Differential \u00B7 Proteomics \u00B7 limma",
                 "G2 vs G1, age-adjusted \u00B7 just now"),
          bullet("--ok",        "Imported RNA-seq experiment",
                 sprintf("%d samples \u00B7 just now",
                         ncol(example_input("rnaseq")$expr_mat))),
          bullet("--ok",        "Imported Proteomics experiment",
                 sprintf("%d samples \u00B7 just now",
                         ncol(example_proteomics_input()$expr_mat))),
          bullet("--accent-500","Created demo project",
                 "Cheek \u00B7 G2 vs G1 \u00B7 just now")
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
            })
          )
        }
      }
    )
  )
}
