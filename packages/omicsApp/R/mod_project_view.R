#' Project view module
#'
#' Top of the workflow funnel: tells the user which project they're
#' looking at, how many samples and features it covers, and what
#' has been run on it. Slice 2C renders a static demo project built
#' from `example_project()`; later slices will swap that for a
#' reactive project handle.
#'
#' Reference markup: `omicsApp/mockup/index.html:462-562`.
#'
#' @param id Module namespace id.
#' @keywords internal
#' @noRd
project_view_ui <- function(id) {
  ns <- shiny::NS(id)
  project <- example_project()

  experiments <- project$experiments
  n_experiments <- length(experiments)
  n_samples_per_exp <- vapply(experiments, function(x) ncol(x$expr_mat), integer(1))
  n_features_per_exp <- vapply(experiments, function(x) nrow(x$expr_mat), integer(1))
  total_features <- sum(n_features_per_exp)

  htmltools::tagList(
    view_header(
      title    = "Project overview",
      subtitle = htmltools::tagList(
        project$name,
        htmltools::HTML(" &middot; "),
        htmltools::tags$span(class = "muted", "demo project (built-in)")
      )
    ),
    htmltools::tags$div(
      class = "stat-grid",
      stat_card(
        label  = "Experiments",
        value  = n_experiments,
        trend  = paste(vapply(experiments, project_omics_label, character(1)),
                       collapse = " + "),
        accent = "brand"
      ),
      stat_card(
        label = "Samples",
        value = unname(n_samples_per_exp[1L]),
        trend = sprintf("%d per omics layer", unname(n_samples_per_exp[1L])),
        mono  = TRUE
      ),
      stat_card(
        label = "Features (total)",
        value = format(total_features, big.mark = ","),
        trend = paste(sprintf("%s: %s",
                              vapply(experiments, project_omics_label, character(1)),
                              format(n_features_per_exp, big.mark = ",")),
                      collapse = " \u00B7 "),
        mono  = TRUE
      ),
      stat_card(
        label  = "Analyses cached",
        value  = 1L,
        trend  = "1 diff (limma, G2 vs G1)",
        accent = "ok"
      )
    ),
    htmltools::tags$div(
      class = "row-grid r-7-5",
      project_experiments_card(experiments),
      project_activity_card()
    )
  )
}

#' @rdname project_view_ui
#' @keywords internal
#' @noRd
project_view_server <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {
    # Slice 2C: all content is static UI. Reactives land in 2D+.
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
    )
  )
}

# Static activity bullets — a real activity log would need persistence,
# which is well outside this slice's scope. Keeping these literal so the
# Project page renders something close to the mockup.
project_activity_card <- function() {
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
  )
}
