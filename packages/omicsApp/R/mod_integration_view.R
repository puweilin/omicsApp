#' Multi-omics integration view module
#'
#' Slice 3E: when the project contains >=2 experiments AND the
#' Diff view has emitted a bundle, we auto-run a second diff on
#' the secondary layer (with the same contrast as the primary)
#' and call `omicsCore::run_integration(method = "concordance")`.
#' When prerequisites aren't met (single experiment, no diff
#' bundle yet) we fall back to `example_integration_tables()`
#' and flag the header.
#'
#' Re-run is gated by an actionButton per slice-3 convention.
#'
#' Reference markup: `omicsApp/mockup/index.html:967-1055`.
#'
#' @param id Module namespace id.
#' @keywords internal
#' @noRd
integration_view_ui <- function(id) {
  ns <- shiny::NS(id)

  htmltools::tagList(
    shiny::uiOutput(ns("header")),
    shiny::uiOutput(ns("notices")),
    shiny::uiOutput(ns("stats")),
    htmltools::tags$div(
      class = "row-grid r-6-6",
      integration_dual_card(ns),
      integration_scatter_card(ns)
    ),
    htmltools::tags$div(
      class = "row-grid r-6-6",
      integration_shared_card(ns),
      integration_ap_card(ns)
    ),
    htmltools::tags$div(
      style = "margin-top:8px",
      shiny::actionButton(
        ns("rerun"), "Re-run integration",
        class = "btn btn-primary"
      )
    )
  )
}

#' @rdname integration_view_ui
#' @param current_project Reactive yielding the live `omics_project`.
#' @param diff_bundle Reactive yielding the primary differential bundle.
#' @keywords internal
#' @noRd
integration_view_server <- function(id,
                                    current_project = shiny::reactiveVal(NULL),
                                    diff_bundle = shiny::reactiveVal(NULL),
                                    invalidate = shiny::reactiveVal(0L)) {
  shiny::moduleServer(id, function(input, output, session) {

    # Decide whether we can run real concordance integration. We
    # need: (1) a project with >=2 experiments, (2) a diff_bundle
    # whose tag matches one of those experiments, (3) the diff
    # bundle's contrast can be applied to the second experiment
    # (same group_col + control + case levels in its meta_df).
    can_run <- shiny::reactive({
      proj <- current_project()
      bundle <- diff_bundle()
      if (is.null(proj) || is.null(bundle)) return(list(ok = FALSE))
      exps <- proj$experiments
      if (length(exps) < 2L) return(list(ok = FALSE))
      types <- vapply(exps, function(e) e$omics_type %||% "", character(1))
      primary_idx <- which(types == bundle$input_info$omics_type)
      if (length(primary_idx) == 0L) primary_idx <- 1L
      secondary_idx <- setdiff(seq_along(exps), primary_idx[1L])[1L]
      params <- bundle$params
      gc <- params$group_col
      ctrl <- params$control_group
      case <- params$case_group
      if (is.null(gc) || is.null(ctrl) || is.null(case)) {
        return(list(ok = FALSE))
      }
      sec <- exps[[secondary_idx]]
      if (!gc %in% names(sec$meta_df)) return(list(ok = FALSE))
      lv <- unique(as.character(stats::na.omit(sec$meta_df[[gc]])))
      if (!all(c(ctrl, case) %in% lv)) return(list(ok = FALSE))
      list(
        ok            = TRUE,
        primary_tag   = names(exps)[primary_idx[1L]],
        secondary_tag = names(exps)[secondary_idx],
        secondary     = sec,
        group_col     = gc,
        control       = ctrl,
        case          = case
      )
    })

    integration_bundle <- shiny::reactiveVal(NULL)
    integration_error  <- shiny::reactiveVal(NULL)
    is_demo            <- shiny::reactiveVal(TRUE)

    # One of the layers this integration spanned has been replaced, so
    # the pairing it reports no longer exists. Back to the module's own
    # start-up state.
    shiny::observeEvent(invalidate(), {
      integration_bundle(NULL)
      integration_error(NULL)
      is_demo(TRUE)
    }, ignoreInit = TRUE)

    do_run <- function() {
      info <- can_run()
      if (!isTRUE(info$ok)) {
        is_demo(TRUE)
        integration_bundle(NULL)
        integration_error(NULL)
        return(invisible())
      }
      proj <- current_project()
      primary <- diff_bundle()
      # Second diff bundle: same contrast, same method when the
      # backend supports the secondary's omics_type, else "auto".
      sec_method <- if (info$secondary$omics_type ==
                        primary$input_info$omics_type) {
        primary$params$method %||% "auto"
      } else {
        "auto"
      }
      run_async(
        function() {
          sec_diff <- omicsCore::run_diff(
            input         = info$secondary,
            method        = sec_method,
            analysis_type = "group",
            group_col     = info$group_col,
            control_group = info$control,
            case_group    = info$case
          )
          diff_bundles <- stats::setNames(
            list(primary, sec_diff),
            c(info$primary_tag, info$secondary_tag)
          )
          omicsCore::run_integration(
            project      = proj,
            method       = "concordance",
            experiments  = c(info$primary_tag, info$secondary_tag),
            diff_bundles = diff_bundles
          )
        },
        on_success = function(result) {
          integration_error(NULL)
          is_demo(FALSE)
          integration_bundle(result)
        },
        on_error = function(msg) {
          integration_error(msg)
          is_demo(TRUE)
          integration_bundle(NULL)
        },
        message = "Running multi-omics integration..."
      )
    }

    shiny::observeEvent(input$rerun, do_run())
    shiny::observeEvent(can_run(), {
      # Auto-run as soon as prerequisites become true; observers
      # higher up (project change, new diff bundle) drive this.
      # `ignoreInit = TRUE` keeps this quiet on session start —
      # otherwise the first evaluation of can_run() (=FALSE) would
      # fire `is_demo(TRUE)` unnecessarily and, combined with a
      # Re-run click in the same flush, could double-run do_run().
      if (isTRUE(can_run()$ok)) do_run() else is_demo(TRUE)
    }, ignoreInit = TRUE)

    # Data sources shared by stats + plots. `req(b)` guards against
    # the transient NULL window between an async run's start and
    # completion in live mode.
    conc_df <- shiny::reactive({
      if (isTRUE(is_demo())) return(example_integration_tables()$concordance_df)
      b <- integration_bundle()
      shiny::req(b)
      b$results$integration_df
    })

    output$header <- shiny::renderUI({
      info <- can_run()
      demo <- is_demo()
      subtitle <- if (demo) {
        "Proteomics \u00D7 RNA-seq \u00B7 demo fixture \u00B7 60 paired features"
      } else {
        sprintf("%s \u00D7 %s \u00B7 %s",
                info$primary_tag, info$secondary_tag,
                paste(info$case, "vs", info$control))
      }
      view_header(
        title    = "Multi-omics integration",
        subtitle = subtitle,
        actions  = htmltools::tagList(
          htmltools::tags$button(type = "button",
                                 class = "btn btn-ghost",
                                 "Mapping"),
          htmltools::tags$button(type = "button",
                                 class = "btn btn-ghost",
                                 "Export integration bundle")
        )
      )
    })

    output$notices <- shiny::renderUI({
      tagged <- htmltools::tagList()
      err <- integration_error()
      if (!is.null(err)) {
        tagged <- htmltools::tagAppendChild(
          tagged,
          notice(title  = "run_integration failed",
                 detail = err,
                 kind   = "warn")
        )
      }
      info <- can_run()
      if (!isTRUE(info$ok)) {
        proj <- current_project()
        reason <- if (is.null(proj) || length(proj$experiments) < 2L) {
          "Integration needs two experiments in the project."
        } else if (is.null(diff_bundle())) {
          "Run a differential analysis first."
        } else {
          paste0("The second layer doesn't have matching ",
                 "group / control / case levels.")
        }
        tagged <- htmltools::tagAppendChild(
          tagged,
          notice(
            title  = "Showing the demo fixture",
            detail = reason,
            kind   = "info"
          )
        )
      }
      tagged
    })

    output$stats <- shiny::renderUI({
      df <- conc_df()
      shiny::req(nrow(df) > 0L)
      paired_n <- nrow(df)
      # `quadrant` is required by the integration result schema, so both
      # a live bundle and the demo fixture carry it.
      quad_col <- df$quadrant
      up_n   <- sum(quad_col == "up_up",     na.rm = TRUE)
      down_n <- sum(quad_col == "down_down", na.rm = TRUE)
      a_col <- if ("effect_a" %in% names(df)) df$effect_a else df$effect
      b_col <- if ("effect_b" %in% names(df)) df$effect_b else NA_real_
      rho <- if (length(b_col) > 1L && !all(is.na(b_col))) {
        sig <- quad_col %in% c("up_up", "down_down", "up_down", "down_up")
        suppressWarnings(stats::cor(a_col[sig], b_col[sig],
                                    method = "spearman",
                                    use = "pairwise.complete.obs"))
      } else NA_real_
      htmltools::tags$div(
        class = "stat-grid",
        stat_card(
          label  = "Paired features",
          value  = format(paired_n, big.mark = ","),
          trend  = "RNA \u2194 protein mapping",
          accent = "brand", mono = TRUE
        ),
        stat_card(
          label  = "Concordant \u2191",
          value  = up_n,
          trend  = "both layers up",
          accent = "up", mono = TRUE
        ),
        stat_card(
          label  = "Concordant \u2193",
          value  = down_n,
          trend  = "both layers down",
          accent = "down", mono = TRUE
        ),
        stat_card(
          label  = "Spearman \u03C1 (sig)",
          value  = if (is.na(rho)) "\u2014" else sprintf("%.2f", rho),
          trend  = "log2FC RNA vs Protein",
          mono   = TRUE
        )
      )
    })

    # As with the other views, demo and live share one dispatcher.
    plot_bundle <- shiny::reactive({
      if (isTRUE(is_demo())) example_integration_bundle()
      else integration_bundle()
    })

    output$dual <- shiny::renderPlot({
      b <- plot_bundle()
      shiny::req(b)
      omicsCore::plot_integration(b, view = "dual_volcano")
    })

    output$scatter <- shiny::renderPlot({
      b <- plot_bundle()
      shiny::req(b)
      omicsCore::plot_integration(b, view = "effect_pair")
    })

    # ActivePathways panels stay on the synthetic fixture in
    # slice 3E (concordance result doesn't carry per-pathway p
    # tables, and active_pathways needs its own Bioconductor
    # Suggests). Display them from `example_integration_tables`
    # when no live AP run is available.
    ap_df <- shiny::reactive({
      example_integration_tables()$active_pathways_df
    })

    output$shared <- shiny::renderPlot({
      ap <- ap_df()
      long <- data.frame(
        pathway_name = rep(ap$pathway_name, times = 2L),
        omics        = rep(c("A", "B"), each = nrow(ap)),
        neglog10p    = c(-log10(ap$p_a), -log10(ap$p_b)),
        direction    = rep(ap$direction, times = 2L),
        stringsAsFactors = FALSE
      )
      long$pathway_name <- factor(long$pathway_name,
                                  levels = rev(ap$pathway_name))
      long$omics <- factor(long$omics, levels = c("A", "B"))
      ggplot2::ggplot(
        long,
        ggplot2::aes(x = .data$omics,
                     y = .data$pathway_name,
                     size = .data$neglog10p,
                     color = .data$direction)
      ) +
        ggplot2::geom_point() +
        ggplot2::scale_color_manual(values = c(shared = omics_colors$shared,
                                               unique = omics_colors$unique_)) +
        ggplot2::scale_size_continuous(range = c(3, 10)) +
        ggplot2::labs(x = NULL, y = NULL,
                      size = "-log10(p)", color = NULL) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
          panel.grid.minor = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_text(color = omics_colors$fg_dark)
        )
    })

    output$ap_table <- DT::renderDT({
      ap <- ap_df()
      out <- data.frame(
        Pathway          = ap$pathway_name,
        `p (A)`          = signif(ap$p_a, 3),
        `p (B)`          = signif(ap$p_b, 3),
        `p (combined)`   = signif(ap$p_combined, 3),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      DT::datatable(
        out,
        rownames  = FALSE,
        selection = "single",
        options   = list(
          pageLength = 10,
          dom        = "tip",
          scrollX    = TRUE,
          columnDefs = list(list(className = "dt-right",
                                 targets = c(1, 2, 3)))
        )
      )
    }, server = TRUE)

    # Expose the integration bundle for slice 3F.
    shiny::reactive(integration_bundle())
  })
}

# ---- internal helpers ------------------------------------------------

utils::globalVariables(".data")

`%||%` <- function(a, b) if (is.null(a)) b else a

# The quadrant palette and the effect-recovery helpers moved to
# omicsCore alongside `plot_integration()`, which is now the only thing
# that needed them.

integration_dual_card <- function(ns) {
  bslib::card(
    bslib::card_header(
      htmltools::tags$h3(class = "card-title", "Mirrored volcano"),
      htmltools::tags$span(class = "card-sub",
                           "x = log2FC(A) - log2FC(B)")
    ),
    bslib::card_body(
      shiny::plotOutput(ns("dual"), height = "320px")
    )
  )
}

integration_scatter_card <- function(ns) {
  bslib::card(
    bslib::card_header(
      htmltools::tags$h3(class = "card-title", "log2FC concordance"),
      htmltools::tags$span(class = "card-sub",
                           "A vs B \u00B7 dashed = y=x")
    ),
    bslib::card_body(
      shiny::plotOutput(ns("scatter"), height = "320px")
    )
  )
}

integration_shared_card <- function(ns) {
  bslib::card(
    bslib::card_header(
      htmltools::tags$h3(class = "card-title", "Shared pathway dotplot"),
      htmltools::tags$span(class = "card-sub",
                           "rows = pathway \u00B7 cols = omics \u00B7 demo fixture")
    ),
    bslib::card_body(
      shiny::plotOutput(ns("shared"), height = "340px")
    )
  )
}

integration_ap_card <- function(ns) {
  bslib::card(
    bslib::card_header(
      htmltools::tags$h3(class = "card-title",
                         "ActivePathways \u00B7 combined p"),
      htmltools::tags$span(class = "card-sub",
                           "Brown's method \u00B7 demo fixture")
    ),
    bslib::card_body(
      DT::DTOutput(ns("ap_table"))
    )
  )
}
