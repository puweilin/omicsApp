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
                                navigate = NULL,
                                replace_project = NULL) {
  shiny::moduleServer(id, function(input, output, session) {

    # A whole new project (Open, Restore, the example): through the
    # caller, which tells the analysis views to let go of what they
    # hold. Setting current_project() alone left them holding the
    # previous project's results whenever the new one was built from the
    # same files, and those results then replaced the ones it was saved
    # with.
    set_project <- function(p) {
      if (is.function(replace_project)) replace_project(p) else current_project(p)
    }

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
            label  = "Layers",
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
          label  = "Layers",
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
    # The welcome card, the workflow checklist and the example project
    # (mod_project_guide.R).
    project_guide_server(input, output, session, current_project, navigate, set_project)

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
      # And the results computed on it: the report and the script
      # carried them as if the layer were still there.
      proj$bundles <- drop_layer_bundles(proj$bundles, tag, proj)
      proj$experiments[[tag]] <- NULL
      # A link naming a layer that is gone would pair samples to nothing.
      if (!is.null(proj$sample_link) && nrow(proj$sample_link) > 0L) {
        proj$sample_link <- proj$sample_link[proj$sample_link$tag != tag, ,
                                             drop = FALSE]
      }
      # Likewise a feature mapping table with a column for that layer.
      if (!is.null(proj$feature_link) && tag %in% names(proj$feature_link)) {
        proj$feature_link <- NULL
      }
      current_project(proj)
      pending_drop(NULL)
      shiny::showNotification(sprintf("Removed layer '%s'.", tag),
                              type = "message")
    })

    # ---- saved-project store -----------------------------------------
    # My projects, Save / restore, and the autosave (mod_project_saved.R).
    project_saved_server(input, output, session, current_project, set_project)
  })
}

# ---- internal helpers ------------------------------------------------

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
      htmltools::tags$h3(class = "card-title", "Omics layers"),
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
          "No layers yet. Use the Import view to upload a file."
        )
      } else {
        htmltools::tags$table(
          class = "tbl",
          htmltools::tags$thead(
            htmltools::tags$tr(
              htmltools::tags$th("Layer"),
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
                sprintf("Imported %s layer",
                        project_omics_label(exp)),
                sprintf("layer %s \u00B7 %d samples \u00B7 %d features",
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

# The bundles computed on a layer: by the layer they record
# (omicsCore::bundle_layer(), read against `project` as it was with the
# layer still in it), and for integration by the pair of layers it names.
drop_layer_bundles <- function(bundles, tag, project) {
  if (!length(bundles)) return(bundles)
  keep <- vapply(bundles, function(b) {
    if (!omicsCore::is_analysis_bundle(b)) return(TRUE)
    if (tag %in% (b$params$experiments %||% character(0))) return(FALSE)
    !identical(omicsCore::bundle_layer(project, b), tag)
  }, logical(1))
  bundles[keep]
}
