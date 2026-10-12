# The getting-started half of the project view: the welcome card, the
# workflow checklist and the example project.
#
# Split out of mod_project_view.R. The server half is a plain function
# called from inside project_view_server()'s moduleServer(), not a
# module, so the ids stay "project-load_tutorial", "project-go_qc".
project_guide_server <- function(input, output, session, current_project, navigate,
                                 set_project) {
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
    set_project(tutorial_project())
    shiny::showNotification(
      "Example project loaded: proteomics + RNA-seq, Control vs TreatA / TreatB.",
      type = "message")
  })

  shiny::observeEvent(input$confirm_tutorial, {
    shiny::removeModal()
    set_project(tutorial_project())
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
  invisible()
}

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
  # A step is done when its result exists and the user has been to it.
  # Results that appear on their own -- QC on import, Enrichment and
  # Integration after a differential run -- used to tick their steps,
  # and the tutorial skipped the tips for three of its six steps.
  # Projects saved before visits were recorded count results alone.
  visited <- proj$visited_steps
  done <- vapply(seq_len(nrow(steps)), function(i) {
    if (steps$id[i] == "import") return(n_layers > 0L)
    if (steps$id[i] == "report") return(FALSE)
    !is.na(steps$bundle[i]) && steps$bundle[i] %in% have &&
      (is.null(visited) || steps$id[i] %in% visited)
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
