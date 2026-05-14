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
    shiny::uiOutput(ns("body"))
  )
}

#' @rdname project_view_ui
#' @param current_project Reactive (or reactiveVal) yielding the
#'   live `omics_project` or `NULL`.
#' @keywords internal
#' @noRd
project_view_server <- function(id, current_project = shiny::reactiveVal(NULL)) {
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
        project_experiments_card(r$project$experiments),
        project_activity_card(r$project, is_demo = r$is_demo)
      )
    })
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

project_experiments_card <- function(experiments) {
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
            lapply(names(experiments), function(tag) {
              exp <- experiments[[tag]]
              htmltools::tags$tr(
                htmltools::tags$td(htmltools::tags$span(class = "text-mono", tag)),
                htmltools::tags$td(project_omics_label(exp)),
                htmltools::tags$td(class = "num", ncol(exp$expr_mat)),
                htmltools::tags$td(class = "num",
                                   format(nrow(exp$expr_mat), big.mark = ",")),
                htmltools::tags$td(pill("ready", kind = "ok")),
                htmltools::tags$td(htmltools::tags$a("View \u2192"))
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
                         ncol(example_input("proteomics")$expr_mat))),
          bullet("--accent-500","Created demo project",
                 "CHISSS \u00B7 Cheek \u00B7 G2 vs G1 \u00B7 just now")
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
