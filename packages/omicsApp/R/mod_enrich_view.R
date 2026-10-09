#' Pathway enrichment view module
#'
#' Slice 3E: takes a live `diff_bundle` from the Differential
#' view and runs `omicsCore::run_enrichment()` against it
#' (ORA or GSEA, against an MSigDB database). The Re-run button
#' gates the call per slice-3 convention. When clusterProfiler is
#' missing we fall back to `example_enrich_table()` and surface a
#' notice with the install hint.
#'
#' Reference markup: `omicsApp/mockup/index.html:918-964`.
#'
#' @param id Module namespace id.
#' @keywords internal
#' @noRd
enrich_view_ui <- function(id) {
  ns <- shiny::NS(id)

  htmltools::tagList(
    shiny::uiOutput(ns("header")),
    shiny::uiOutput(ns("notices")),
    shiny::uiOutput(ns("input_summary")),
    # Parameters across the top rather than down a 3fr rail: three
    # radio groups and a dropdown do not need a quarter of the width,
    # and the dotplot -- whose y axis carries pathway names -- needed
    # every bit of it. The two result cards then split the full width.
    enrich_params_card(ns),
    htmltools::tags$div(
      class = "row-grid r-7-5",
      enrich_dot_card(ns),
      enrich_hits_card(ns)
    ),
    # The pathway picked in the table (mod_enrich_results.R).
    shiny::uiOutput(ns("pathway_detail")),
    # Only when the differential run holds several comparisons.
    shiny::uiOutput(ns("compare_card"))
  )
}

#' @rdname enrich_view_ui
#' @param diff_bundle Reactive (or reactiveVal) yielding the
#'   live differential `analysis_bundle` from the Diff view, or
#'   `NULL`.
#' @keywords internal
#' @noRd
enrich_view_server <- function(id, diff_bundle = shiny::reactiveVal(NULL),
                               diff_thresholds = shiny::reactive(list(
                                 p_cutoff = 0.05, p_preference = "adjusted",
                                 effect_cutoff = NULL)),
                               diff_layer = shiny::reactive(NULL),
                               invalidate = shiny::reactiveVal(0L),
                               navigate = NULL,
                               diff_all = shiny::reactive(NULL),
                               current_project = shiny::reactive(NULL)) {
  shiny::moduleServer(id, function(input, output, session) {

    have_cp <- has_pkg("clusterProfiler")

    # The species the gene sets come from. Suggested from the symbols --
    # mouse genes are written Trp53, human TP53, zebrafish tp53 -- because
    # the default human sets matched no mouse gene, and the result was an
    # empty table with no word of why. Every species msigdbr carries the
    # MSigDB sets to (through orthologs) is offered; mouse and rat write
    # their genes alike, so a rat table is suggested mouse and switched by
    # hand.
    suggested_organism <- shiny::reactive({
      b <- diff_bundle()
      if (is.null(b)) return("Hs")
      guess_organism(b$results$diff_result_df$feature_symbol)
    })
    organism <- shiny::reactive(input$organism %||% suggested_organism())

    # The gene list(s) to test, as the Direction control says -- except
    # for a result with no direction (a spline fit), which has one list
    # only: it is enriched pooled, whatever the control is set to.
    run_direction <- function(bundle) {
      if (diff_undirected(bundle)) return("both")
      input$direction %||% "separate"
    }

    # ---- every comparison side by side --------------------------------
    # Every comparison of the run enriched with the same settings
    # (mod_enrich_compare.R).
    compare <- enrich_compare_server(input, output, session, diff_all, diff_thresholds,
                                     have_cp, organism, run_direction)
    compare_bundle  <- compare$compare_bundle
    compare_error   <- compare$compare_error
    all_comparisons <- compare$all_comparisons
    compare_epoch   <- compare$compare_epoch
    compare_running <- compare$compare_running

    output$ui_organism <- shiny::renderUI({
      sel <- shiny::isolate(input$organism) %||% suggested_organism()
      shiny::selectInput(session$ns("organism"), label = "Species",
                         choices = species_choices(), selected = sel)
    })
    shiny::observeEvent(suggested_organism(), {
      shiny::updateSelectInput(session, "organism", selected = suggested_organism())
    }, ignoreInit = TRUE)

    enrich_bundle <- shiny::reactiveVal(NULL)
    enrich_error  <- shiny::reactiveVal(NULL)
    is_demo       <- shiny::reactiveVal(TRUE)
    # No demo pathways once the user has a project of their own.
    shiny::observeEvent(current_project(), {
      if (is.null(enrich_bundle())) is_demo(is.null(current_project()) || !have_cp)
    }, ignoreNULL = FALSE)

    # The layer the upstream diff was computed on has been replaced, so
    # this enrichment is no longer about anything in the project. Back
    # to the module's own start-up state.
    enrich_epoch <- run_epoch()
    shiny::observeEvent(invalidate(), {
      enrich_epoch$bump()
      compare_epoch$bump()
      enrich_bundle(NULL)
      enrich_error(NULL)
      is_demo(is.null(current_project()))
      compare_bundle(NULL)
    }, ignoreInit = TRUE)

    do_run <- function() {
      my_run <- enrich_epoch$start()
      db_arg <- input$database %||% "hallmark"
      type   <- input$type %||% "ora"
      thr    <- diff_thresholds()
      bundle <- diff_bundle()
      dir_   <- run_direction(bundle)
      if (is.null(bundle) || !have_cp) {
        # Demo fallback: synthetic table re-shaped to match the
        # standardized enrich schema.
        enrich_error(if (!have_cp) {
          paste0("clusterProfiler is not installed; showing the demo ",
                 "fixture. Install with `omicsCore::install_optional",
                 "('enrichment')`.")
        } else NULL)
        # The demo only when there is no project of the user's: its
        # pathways under a restored project read as the user's result.
        is_demo(is.null(current_project()) || !have_cp)
        enrich_bundle(NULL)
        return(invisible())
      }
      set_button_busy("rerun", TRUE)
      run_async(
        # Detached for the same reason as the Differential view: a
        # closure defined here carries this module's whole scope to the
        # worker, and an enrichment runs on a bundle that is already
        # large.
        detached_call(
          function() {
            # ORA takes the thresholds; GSEA does not.
            #
            # GSEA is a ranked-list method: it reads every feature,
            # ordered by effect, and pre-filtering the list is not a
            # stricter version of the analysis but a different and
            # invalid one. run_enrichment() also spends p_cutoff on
            # bounding the *pathway* table, so handing it a threshold
            # chosen for features would quietly narrow which pathways
            # are reported -- two different tests on two different
            # objects, sharing one number.
            if (identical(type, "ora")) {
              omicsCore::run_enrichment(
                diff_bundle   = bundle,
                type          = type,
                database      = db_arg,
                organism      = org,
                direction     = dir_,
                p_cutoff      = thr$p_cutoff,
                p_preference  = thr$p_preference,
                effect_cutoff = thr$effect_cutoff
              )
            } else {
              omicsCore::run_enrichment(
                diff_bundle = bundle,
                type        = type,
                database    = db_arg,
                organism    = org,
                direction   = dir_
              )
            }
          },
          bundle = bundle, type = type, db_arg = db_arg, dir_ = dir_,
          thr = thr, org = organism()
        ),
        on_success = function(result) {
          if (enrich_epoch$is_last_started(my_run)) set_button_busy("rerun", FALSE)
          if (!enrich_epoch$is_current(my_run)) return(invisible())
          enrich_error(NULL)
          is_demo(FALSE)
          enrich_bundle(result)
        },
        on_error = function(msg) {
          if (enrich_epoch$is_last_started(my_run)) set_button_busy("rerun", FALSE)
          if (!enrich_epoch$is_current(my_run)) return(invisible())
          # The error, and nothing under it: the demo's pathways beneath
          # a failure notice read as the user's result.
          enrich_error(msg)
          is_demo(FALSE)
          enrich_bundle(NULL)
        },
        message = "Running pathway enrichment..."
      )
    }

    # Auto-run once on first diff_bundle change so the view lands
    # populated; subsequent updates require Re-run.
    # Also when the upstream result goes away (a layer change clears it):
    # with ignoreNULL the old enrichment stayed on screen, now labelled
    # with the new layer, and went into the project without its diff.
    shiny::observeEvent(diff_bundle(), {
      # A project opened or restored brings its enrichment with it: shown
      # as it was saved, not recomputed at whatever the controls say.
      saved <- current_project()$bundles$enrich
      b <- diff_bundle()
      if (!is.null(saved) && !is.null(b) && is.null(enrich_bundle()) &&
          identical(saved$params$comparison, b$params$comparison) &&
          identical(omicsCore::bundle_layer(current_project(), saved),
                    omicsCore::bundle_layer(current_project(), b))) {
        enrich_epoch$bump()
        restore_controls(saved$params)
        enrich_error(NULL)
        is_demo(FALSE)
        enrich_bundle(saved)
        sc <- current_project()$bundles$enrich_compare
        if (!is.null(sc) && is.null(compare_bundle())) compare_bundle(sc)
        return(invisible())
      }
      do_run()
    }, ignoreInit = TRUE, ignoreNULL = FALSE)

    restore_controls <- function(params) {
      if (!is.null(params$type)) shiny::updateRadioButtons(session, "type", selected = params$type)
      if (!is.null(params$database)) shiny::updateSelectInput(session, "database",
                                                              selected = params$database[[1L]])
      if (!is.null(params$direction)) shiny::updateRadioButtons(session, "direction",
                                                                selected = params$direction)
      org <- species_code(params$organism)
      if (!is.null(org)) shiny::updateSelectInput(session, "organism", selected = org)
    }

    # The result on screen no longer answers the controls: the test,
    # database, species or lists changed, or (for ORA) the thresholds
    # that chose its genes did. Read from the result's own parameters,
    # so a restored result is judged the same way as a fresh one.
    settings_changed <- shiny::reactive({
      b <- enrich_bundle()
      if (is.null(b) || isTRUE(is_demo())) return(FALSE)
      p <- b$params
      type <- input$type %||% "ora"
      same_num <- function(a, b) {
        a <- if (length(a)) as.numeric(a[[1L]]) else NA_real_
        b <- if (length(b)) as.numeric(b[[1L]]) else NA_real_
        (is.na(a) && is.na(b)) || isTRUE(all.equal(a, b))
      }
      changed <- !identical(type, p$type %||% type) ||
        !identical(input$database %||% "hallmark", p$database[[1L]] %||% "hallmark") ||
        !identical(organism(), species_code(p$organism) %||% organism())
      if (identical(p$type, "ora")) {
        thr <- diff_thresholds()
        changed <- changed ||
          !identical(run_direction(diff_bundle()), p$direction %||% "both") ||
          !same_num(thr$p_cutoff, p$p_cutoff) ||
          !identical(thr$p_preference %||% "adjusted", p$p_preference %||% "adjusted") ||
          !same_num(thr$effect_cutoff, p$effect_cutoff)
      }
      changed
    })

    shiny::observeEvent(input$rerun, do_run())
    if (is.function(navigate)) {
      shiny::observeEvent(input$go_next, navigate("integration"))
    }

    # The header, the gene-list summary, the notices, the dot plot and
    # the table (mod_enrich_results.R).
    results <- enrich_results_server(input, output, session, navigate, diff_bundle,
                                     diff_layer, diff_thresholds, current_project, have_cp,
                                     enrich_bundle, enrich_error, is_demo, settings_changed)
    table_data        <- results$table_data
    diff_layer_tag    <- results$diff_layer_tag
    selected_features <- results$selected_features
    plot_bundle       <- results$plot_bundle
    show_p            <- results$show_p
    show_cutoff       <- results$show_cutoff

    # Expose the bundle for slice 3F (report).
    list(
      bundle = shiny::reactive(enrich_bundle()),
      compare = shiny::reactive(compare_bundle())
    )
  })
}

# ---- internal helpers ------------------------------------------------

utils::globalVariables(".data")

enrich_params_card <- function(ns) {
  bslib::card(
    # A dropdown opened from this card renders *inside* it, and a card
    # clips its overflow -- so the Database menu was cut off at the card
    # edge. Only an issue since the controls moved to a top row: down a
    # tall rail there was always card below the menu.
    class = "param-row-card",
    bslib::card_header(
      htmltools::tags$h3(class = "card-title", "Parameters"),
      htmltools::tags$span(class = "card-sub",
                           "test type, database, direction")
    ),
    bslib::card_body(
      htmltools::tags$div(
        class = "param-row",
        param_group(
          "Test",
          help = paste("ORA tests whether the differential hits are over-represented",
                       "in a pathway. GSEA ranks every feature by effect and needs",
                       "no hit threshold."),
          shiny::radioButtons(
            ns("type"), label = NULL,
            choices  = c("ORA" = "ora", "GSEA" = "gsea"),
            selected = "ora", inline = TRUE
          )
        ),
        param_group(
          "Species",
          help = "Which species' gene sets to use. Suggested from how the gene names are written.",
          shiny::uiOutput(ns("ui_organism"))
        ),
        param_group(
          "Database",
          shiny::selectInput(
            ns("database"),
            label = htmltools::tags$span(class = "visually-hidden-label", "Database"),
            choices  = c("MSigDB Hallmark" = "hallmark",
                         "KEGG" = "kegg",
                         "Reactome" = "reactome",
                         "GO Biological Process" = "go_bp",
                         "GO Molecular Function" = "go_mf",
                         "GO Cellular Component" = "go_cc",
                         "WikiPathways" = "wikipathways"),
            selected = "hallmark"
          )
        ),
        param_group(
          "Direction",
          help = paste("Which hits ORA tests. Separately (the default) tests the",
                       "up-regulated and the down-regulated hits as two lists, so",
                       "opposite changes are not mixed in one pathway; pooled tests",
                       "them as one list. GSEA always keeps up and down apart."),
          shiny::radioButtons(
            ns("direction"), label = NULL,
            choices  = ENRICH_DIRECTION_CHOICES,
            selected = "separate", inline = TRUE
          )
        ),
        param_group(
          "Significance",
          # A display threshold, like the Differential view's sliders:
          # applied to the result that is already computed, so switching
          # it costs nothing and re-runs nothing. GSEA over 2695 GO sets
          # keeps 9 at adjusted p and 131 at raw; which of those is the
          # right answer is the reader's call, not ours.
          shiny::radioButtons(
            ns("show_p"), label = NULL,
            choices  = c("adjusted p" = "adjusted", "raw p" = "raw"),
            selected = "adjusted", inline = TRUE
          ),
          shiny::numericInput(
            ns("show_cutoff"),
            label = htmltools::tags$span(class = "visually-hidden-label", "Display p cutoff"),
            value = 0.05, min = 0, max = 1, step = 0.01
          )
        ),
        htmltools::tags$div(
          class = "param-action",
          shiny::actionButton(
            ns("rerun"), "Re-run",
            class = "btn btn-primary"
          )
        )
      )
    )
  )
}

# Supplied to the worker function's environment by detached_call(), which
# codetools cannot see.
utils::globalVariables("bundle")

# The Direction control, in plain words. "both" keeps its old meaning --
# up and down pooled -- so a project saved before "separate" existed
# restores to what it ran.
ENRICH_DIRECTION_CHOICES <- c("Up and down separately" = "separate",
                              "Up only" = "up",
                              "Down only" = "down",
                              "Up and down pooled" = "both")

# A differential result that says whether features change, not which way
# (a global test, a spline fit): it can only be enriched as one list.
diff_undirected <- function(bundle) {
  at <- bundle$results$diff_result_df$analysis_type
  any(at %in% c("anova", "continuous_spline"))
}

# The species offered, as label = code, from omicsCore's own list so the
# menu and the engine cannot drift apart.
species_choices <- function() {
  sp <- omicsCore::enrichment_species()
  stats::setNames(sp$code, sp$label)
}

species_code <- function(organism) {
  if (is.null(organism) || !length(organism) || is.na(organism[[1L]])) return(NULL)
  sp <- omicsCore::enrichment_species()
  hit <- match(organism[[1L]], c(sp$species, sp$code))
  if (is.na(hit)) return(NULL)
  c(sp$code, sp$code)[[hit]]
}

species_label <- function(code) {
  sp <- omicsCore::enrichment_species()
  lab <- sp$label[match(code, sp$code)]
  if (is.na(lab)) code else sub(" \\(.*\\)$", "", lab)
}

species_example <- function(code) {
  switch(code, Hs = "TP53", Mm = "Trp53", Rn = "Tp53", Dr = "tp53",
         Dm = "p53", Sc = "RAD9", Ce = "cep-1", "TP53")
}

# Which species a set of gene symbols most likely comes from, by how the
# names are written: human upper case (TP53), mouse and rat capitalised
# (Trp53 -- the two cannot be told apart this way, so mouse is
# suggested), zebrafish lower case (tp53), worm lower case with a dash and
# number (unc-54), fly with many CG-numbered genes (CG1824), yeast with
# systematic ORF names (YAL001C).
guess_organism <- function(symbols) {
  x <- symbols[is_gene_symbol(symbols)]
  if (length(x) < 10L) return("Hs")
  orf_like <- mean(grepl("^Y[A-P][LR][0-9]{3}[WC](-[A-Z])?$", x))
  if (orf_like > 0.3) return("Sc")
  cg_like <- mean(grepl("^C[GR][0-9]+$", x))
  if (cg_like > 0.2) return("Dm")
  worm_like <- mean(grepl("^[a-z]{2,5}-[0-9]+(\\.[0-9]+)?$", x))
  if (worm_like > 0.5) return("Ce")
  lower_like <- mean(x == tolower(x) & grepl("[a-z]", x))
  mouse_like <- mean(grepl("^[A-Z][a-z0-9]+[a-z0-9-]*$", x) & grepl("[a-z]", x))
  human_like <- mean(x == toupper(x))
  if (lower_like > 0.5 && lower_like > human_like) return("Dr")
  if (mouse_like > 0.5 && mouse_like > human_like) "Mm" else "Hs"
}

# Values that look like gene symbols rather than accessions or blanks.
is_gene_symbol <- function(x) {
  x <- as.character(x)
  !is.na(x) & nzchar(x) & grepl("^[A-Za-z][A-Za-z0-9.-]{0,19}$", x) &
    !grepl("^(ENS[A-Z]*[0-9]{6,}|[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]|N[MRP]_|X[MR]_)", x)
}
