#' Multi-omics integration view module
#'
#' Pairs the layer the Differential view ran on (the *primary*) with a
#' second layer of the project (the *partner*) and integrates them by one
#' of three methods:
#'
#' * **Fold-change concordance** -- re-runs the primary's contrast on the
#'   partner (same group column, control, case, covariates and, for a
#'   multi-contrast run, the same set of groups in the model) and compares
#'   the two differential results feature by feature.
#' * **Sample-level correlation** -- correlates each shared gene across
#'   the samples the pairing card says are the same person. Needs no
#'   differential result.
#' * **ActivePathways** -- pathway-level evidence merged across the two
#'   differential results (only offered when the package is installed).
#'
#' With no project the view shows the built-in demo fixture, labelled as
#' such. With a project it never does: a failed or impossible run shows
#' why, not somebody else's numbers.
#'
#' @param id Module namespace id.
#' @keywords internal
#' @noRd
integration_view_ui <- function(id) {
  ns <- shiny::NS(id)

  htmltools::tagList(
    shiny::uiOutput(ns("header")),
    shiny::uiOutput(ns("notices")),
    integration_setup_card(ns),
    # Before the results, not after: which sample was treated as which
    # person is an assumption the reader should see before reading
    # anything computed on it.
    integration_pairing_card(ns),
    shiny::uiOutput(ns("stats")),
    shiny::uiOutput(ns("results"))
  )
}

#' @rdname integration_view_ui
#' @param current_project Reactive yielding the live `omics_project`.
#' @param diff_bundle Reactive yielding the primary differential bundle
#'   (one contrast).
#' @param diff_layer Reactive yielding the project tag of the layer the
#'   differential bundle was computed on.
#' @param diff_thresholds Reactive yielding the Differential view's
#'   `list(p_cutoff, p_preference, effect_cutoff)`, so a "hit" means the
#'   same thing in both views.
#' @keywords internal
#' @noRd
integration_view_server <- function(id,
                                    current_project = shiny::reactiveVal(NULL),
                                    diff_bundle = shiny::reactiveVal(NULL),
                                    invalidate = shiny::reactiveVal(0L),
                                    diff_layer = shiny::reactiveVal(NULL),
                                    diff_thresholds = shiny::reactiveVal(NULL),
                                    navigate = NULL) {
  shiny::moduleServer(id, function(input, output, session) {

    method <- shiny::reactive(input$method %||% "concordance")

    # ---- which two layers --------------------------------------------
    # The primary is the layer the diff ran on -- the Differential view
    # says which. It used to be guessed from the bundle's omics type,
    # which picked the wrong layer whenever two layers shared a type.
    layers <- shiny::reactive({
      proj <- current_project()
      if (is.null(proj) || length(proj$experiments) < 2L) return(NULL)
      tags <- names(proj$experiments)
      primary <- diff_layer()
      if (is.null(primary) || !primary %in% tags) {
        b <- diff_bundle()
        types <- vapply(proj$experiments, function(e) e$omics_type %||% "",
                        character(1))
        hit <- if (!is.null(b)) which(types == (b$input_info$omics_type %||% "")) else integer(0)
        primary <- tags[if (length(hit)) hit[1L] else 1L]
      }
      others <- setdiff(tags, primary)
      want <- input$partner
      partner <- if (!is.null(want) && want %in% others) want else others[[1L]]
      list(primary = primary, partner = partner, others = others)
    })

    output$ui_partner <- shiny::renderUI({
      l <- layers()
      if (is.null(l)) return(NULL)
      if (length(l$others) < 2L) {
        return(htmltools::tags$div(
          class = "muted", style = "font-size:12.5px;padding-top:6px",
          sprintf("%s × %s", l$primary, l$partner)))
      }
      shiny::selectInput(session$ns("partner"),
                         label = sprintf("Integrate %s with", l$primary),
                         choices = l$others,
                         selected = shiny::isolate(l$partner))
    })

    # ---- prerequisites ------------------------------------------------
    # Everything a run needs, or the one reason it cannot happen, in
    # words that name the missing piece.
    can_run <- shiny::reactive({
      proj <- current_project()
      if (is.null(proj) || length(proj$experiments) < 2L) {
        return(list(ok = FALSE, reason = "layers"))
      }
      l <- layers()
      sec <- proj$experiments[[l$partner]]
      base <- list(primary_tag = l$primary, secondary_tag = l$partner,
                   secondary = sec)
      if (identical(method(), "correlation")) {
        return(c(list(ok = TRUE), base))
      }
      bundle <- diff_bundle()
      if (is.null(bundle)) return(list(ok = FALSE, reason = "diff"))
      params <- bundle$params
      gc <- params$group_col
      ctrl <- params$control_group
      case <- params$case_group
      spec <- params$contrasts
      if (!is.null(gc) && length(spec) == 1L) {
        # A contrast written as an expression (a pair from "all pairs", or
        # a weighted one): repeated on the partner as the same expressions.
        if (!gc %in% names(sec$meta_df)) {
          return(list(ok = FALSE, reason = "design", detail = sprintf(
            "Layer '%s' has no '%s' column in its sample metadata, so the comparison %s cannot be repeated on it.",
            l$partner, gc, params$comparison)))
        }
        lv <- unique(as.character(stats::na.omit(sec$meta_df[[gc]])))
        ok_spec <- tryCatch({
          omicsCore::contrast_labels(spec, lv)
          TRUE
        }, error = function(e) conditionMessage(e))
        if (!isTRUE(ok_spec)) {
          return(list(ok = FALSE, reason = "design", detail = sprintf(
            "The comparison %s cannot be repeated on '%s': %s",
            params$comparison, l$partner, ok_spec)))
        }
        # The rest of the run's contrasts, where the partner has their
        # groups, so the partner's model holds the same groups.
        all_specs <- Filter(function(s) isTRUE(tryCatch({
          omicsCore::contrast_labels(s, lv); TRUE
        }, error = function(e) FALSE)), params$all_contrasts %||% spec)
        return(c(list(
          ok             = TRUE,
          group_col      = gc,
          control        = ctrl,
          case           = case,
          contrasts      = {
            keep <- unlist(all_specs)
            if (spec %in% keep) keep else c(keep, spec)
          },
          comparison     = params$comparison,
          covariates     = setdiff(params$covariates, setdiff(params$covariates, names(sec$meta_df))),
          dropped_covs   = setdiff(params$covariates, names(sec$meta_df)),
          primary_method = params$method %||% "auto",
          primary_omics  = bundle$input_info$omics_type
        ), base))
      }
      if (is.null(gc) || is.null(ctrl) || is.null(case)) {
        return(list(ok = FALSE, reason = "design",
                    detail = "The differential result is not a group comparison."))
      }
      if (!gc %in% names(sec$meta_df)) {
        return(list(ok = FALSE, reason = "design", detail = sprintf(
          "Layer '%s' has no '%s' column in its sample metadata, so the contrast %s vs %s cannot be repeated on it.",
          l$partner, gc, case, ctrl)))
      }
      lv <- unique(as.character(stats::na.omit(sec$meta_df[[gc]])))
      missing_lv <- setdiff(as.character(c(ctrl, case)), lv)
      if (length(missing_lv)) {
        return(list(ok = FALSE, reason = "design", detail = sprintf(
          "Layer '%s' has no samples in %s of '%s' (it has: %s).",
          l$partner, paste(sprintf("'%s'", missing_lv), collapse = ", "), gc,
          paste(lv, collapse = ", "))))
      }
      # A multi-contrast diff fitted every case group in one model; the
      # partner is fitted the same way, so the two results rest on the
      # same design. Groups the partner lacks are left out of its model.
      all_cases <- intersect(as.character(params$all_case_groups %||% case), lv)
      covs <- params$covariates
      dropped_covs <- setdiff(covs, names(sec$meta_df))
      c(list(
        ok             = TRUE,
        group_col      = gc,
        control        = ctrl,
        case           = case,
        all_cases      = union(as.character(case), all_cases),
        comparison     = params$comparison,
        covariates     = setdiff(covs, dropped_covs),
        dropped_covs   = dropped_covs,
        primary_method = params$method %||% "auto",
        primary_omics  = bundle$input_info$omics_type
      ), base)
    })

    integration_bundle <- shiny::reactiveVal(NULL)
    integration_error  <- shiny::reactiveVal(NULL)
    is_demo            <- shiny::reactive(is.null(current_project()))
    running            <- shiny::reactiveVal(FALSE)
    # The partner's own diff is the expensive half of a concordance run
    # and depends only on the partner and the design, so it is kept and
    # reused when only the thresholds (or the method) change.
    sec_cache <- shiny::reactiveVal(NULL)

    # ---- sample pairing ------------------------------------------------
    # Shown whether or not the current method needs it. Concordance
    # matches features, not samples, so it runs without a pairing -- but
    # a user who gets that far and then cannot correlate should find out
    # here rather than after the run, and should be able to fix it here
    # rather than by re-importing.
    pairing <- shiny::reactive({
      proj <- current_project()
      l <- layers()
      if (is.null(proj) || is.null(l)) return(NULL)
      omicsCore::sample_pairing_preview(proj, l$primary, l$partner)
    })

    output$pairing_table <- DT::renderDT({
      p <- pairing()
      shiny::req(p, nrow(p$pairs) > 0L)
      l <- layers()
      out <- p$pairs
      names(out) <- c("Donor", l$primary, l$partner)
      DT::datatable(out, rownames = FALSE, selection = "none",
                    options = list(pageLength = 8, dom = "tip"))
    }, server = TRUE)

    output$pairing_note <- shiny::renderUI({
      p <- pairing()
      if (is.null(p)) {
        return(notice("Import a second layer to pair samples across omics.",
                      kind = "info"))
      }
      n <- nrow(p$pairs)
      if (n == 0L) {
        return(notice(
          "No pairing found",
          paste0(
            "Sample-level integration needs to know which sample in each ",
            "layer came from the same person. Add a `donor` column to ",
            "each layer's sample metadata and re-import those layers, or ",
            "rename the samples so they share a leading id — RD001-C ",
            "and RD001_Folli both give RD001. Nothing already computed ",
            "has to be re-run: donor is read by Integration and by ",
            "nothing else. Fold-change concordance does not need a pairing."),
          kind = "warn"))
      }
      n_donor <- length(unique(p$pairs$donor_id))
      dup <- n > n_donor
      extra <- if (dup) {
        sprintf(" %d donor(s) have more than one sample in a layer; correlation uses one sample per donor.",
                n - n_donor)
      } else ""
      switch(
        p$source,
        linked = notice(sprintf("%d pairs, saved with this project.%s", n, extra),
                        kind = if (dup) "warn" else "info"),
        donor = notice(sprintf("%d pairs, from the donor column in both layers.%s", n, extra),
                       kind = if (dup) "warn" else "info"),
        sample_id = notice(sprintf("%d pairs, from sample ids that match outright.", n),
                           kind = "info"),
        suggested = notice(
          sprintf("%d pairs, guessed from the sample ids — check them", n),
          paste0("Nothing states that these are the same person; the ids ",
                 "merely share a leading part. Accept the pairing to keep ",
                 "it with the project (sample-level correlation will not ",
                 "use a guess), or state it properly with a donor column."),
          kind = "warn"),
        notice(sprintf("%d pairs.", n), kind = "info")
      )
    })

    # A guess becomes a decision only when someone says so, and then it is
    # saved with the project rather than re-derived every time.
    output$pairing_action <- shiny::renderUI({
      p <- pairing()
      if (is.null(p) || nrow(p$pairs) == 0L) return(NULL)
      if (!identical(p$source, "suggested")) return(NULL)
      shiny::actionButton(session$ns("accept_pairing"), "Accept pairing",
                          class = "btn btn-sm btn-primary")
    })

    shiny::observeEvent(input$accept_pairing, {
      proj <- current_project()
      p <- pairing()
      l <- layers()
      shiny::req(proj, p, l, nrow(p$pairs) > 0L)
      new_rows <- rbind(
        data.frame(tag = l$primary, sample_id = p$pairs$a,
                   donor_id = p$pairs$donor_id, stringsAsFactors = FALSE),
        data.frame(tag = l$partner, sample_id = p$pairs$b,
                   donor_id = p$pairs$donor_id, stringsAsFactors = FALSE)
      )
      # Added to what the project already states, not in place of it: a
      # third layer's saved pairing is not this pair's to throw away.
      old <- proj$sample_link
      if (!is.null(old) && nrow(old)) {
        clash <- paste(old$tag, old$sample_id) %in%
          paste(new_rows$tag, new_rows$sample_id)
        new_rows <- rbind(old[!clash, c("tag", "sample_id", "donor_id"), drop = FALSE],
                          new_rows)
      }
      rownames(new_rows) <- NULL
      proj$sample_link <- new_rows
      current_project(proj)
      shiny::showNotification(
        sprintf("Pairing saved: %d donors.", nrow(p$pairs)),
        type = "message")
    })

    # One of the layers this integration spanned has been replaced, so
    # the pairing it reports no longer exists. Back to the module's own
    # start-up state.
    shiny::observeEvent(invalidate(), {
      integration_bundle(NULL)
      integration_error(NULL)
      sec_cache(NULL)
      runs$last_key <- NULL
    }, ignoreInit = TRUE)

    # ---- running ------------------------------------------------------
    # Everything a result depends on, as one string. The auto-run fires
    # only when this changes. It used to fire on every change of the
    # project -- including the one that files this view's own result into
    # the project -- so every run was followed by a second, identical one.
    run_key <- shiny::reactive({
      info <- can_run()
      if (!isTRUE(info$ok)) return(NULL)
      proj <- current_project()
      fp <- function(tag) {
        e <- proj$experiments[[tag]]
        e$source_fingerprint %||% paste(dim(e$expr_mat), collapse = "x")
      }
      th <- diff_thresholds() %||% list()
      b <- diff_bundle()
      link <- proj$sample_link
      paste(method(), info$primary_tag, fp(info$primary_tag),
            info$secondary_tag, fp(info$secondary_tag),
            if (!identical(method(), "correlation") && !is.null(b))
              paste(b$params$method, paste(b$params$comparison, collapse = ","),
                    nrow(b$results$diff_result_df),
                    sum(b$results$diff_result_df$p_value, na.rm = TRUE)),
            th$p_cutoff, th$p_preference, th$effect_cutoff,
            if (identical(method(), "correlation") && !is.null(link)) nrow(link),
            sep = "|")
    })
    # Plain state, not reactive: which inputs the last run saw, and a
    # counter that lets a newer run win over an older one still in flight.
    runs <- new.env(parent = emptyenv())
    runs$last_key <- NULL
    runs$token <- 0L

    set_busy <- function(busy) {
      running(busy)
      tryCatch(
        if (busy) shinyjs::disable("rerun") else shinyjs::enable("rerun"),
        error = function(e) NULL)
    }

    do_run <- function() {
      info <- can_run()
      if (!isTRUE(info$ok)) {
        integration_bundle(NULL)
        integration_error(NULL)
        return(invisible())
      }
      runs$last_key <- shiny::isolate(run_key())
      proj <- current_project()
      th <- diff_thresholds() %||% list()
      m <- method()
      runs$token <- runs$token + 1L
      token <- runs$token

      primary <- if (!identical(m, "correlation")) diff_bundle() else NULL
      sec_key <- if (!identical(m, "correlation")) {
        paste(info$secondary_tag,
              info$secondary$source_fingerprint %||% paste(dim(info$secondary$expr_mat), collapse = "x"),
              info$group_col, info$control, paste(info$all_cases, collapse = ","),
              paste(info$covariates, collapse = ","), info$primary_method,
              paste(info$contrasts, collapse = ";"),
              sep = "|")
      }
      cached <- sec_cache()
      sec_diff <- if (!is.null(cached) && identical(cached$key, sec_key)) cached$bundle
      # Same engine on the partner when its data allow it (limma on
      # proteomics and on logCPM RNA-seq alike), else the engine
      # run_diff() would pick for that layer.
      sec_method <- if (!identical(m, "correlation") &&
                        info$primary_method %in% omicsCore::applicable_diff_methods(info$secondary)) {
        info$primary_method
      } else "auto"

      set_busy(TRUE)
      run_async(
        # Detached for the same reason as the other views, and it matters
        # most here: this one carries a project with every layer in it.
        detached_call(
          function() {
            tags <- c(info$primary_tag, info$secondary_tag)
            if (identical(m, "correlation")) {
              res <- omicsCore::run_integration(
                project = proj, method = "correlation", experiments = tags)
              return(list(result = res, sec_diff = NULL))
            }
            if (is.null(sec_diff)) {
              sec_full <- if (length(info$contrasts)) {
                omicsCore::run_diff(
                  input         = info$secondary,
                  method        = sec_method,
                  analysis_type = "group",
                  group_col     = info$group_col,
                  contrasts     = info$contrasts,
                  covariates    = if (length(info$covariates)) info$covariates else NULL
                )
              } else {
                omicsCore::run_diff(
                  input         = info$secondary,
                  method        = sec_method,
                  analysis_type = "group",
                  group_col     = info$group_col,
                  control_group = info$control,
                  case_group    = info$all_cases,
                  covariates    = if (length(info$covariates)) info$covariates else NULL
                )
              }
              sec_diff <- omicsCore::select_comparison(
                sec_full,
                if (length(info$contrasts)) info$comparison
                else paste0(info$case, "_vs_", info$control))
            }
            diff_bundles <- stats::setNames(list(primary, sec_diff), tags)
            args <- list(project = proj, method = m, experiments = tags,
                         diff_bundles = diff_bundles,
                         p_cutoff = th$p_cutoff %||% 0.05)
            if (identical(m, "concordance")) {
              args$p_preference <- th$p_preference %||% "adjusted"
              args$effect_cutoff <- th$effect_cutoff %||% 0
            }
            res <- do.call(omicsCore::run_integration, args)
            list(result = res, sec_diff = sec_diff)
          },
          info = info, m = m, th = th, primary = primary, proj = proj,
          sec_diff = sec_diff, sec_method = sec_method
        ),
        on_success = function(out) {
          # A run overtaken by a newer one does not get to overwrite it.
          if (!identical(token, runs$token)) return()
          set_busy(FALSE)
          integration_error(NULL)
          if (!is.null(out$sec_diff)) {
            sec_cache(list(key = sec_key, bundle = out$sec_diff))
          }
          integration_bundle(out$result)
        },
        on_error = function(msg) {
          if (!identical(token, runs$token)) return()
          set_busy(FALSE)
          integration_error(msg)
          integration_bundle(NULL)
        },
        message = "Running multi-omics integration..."
      )
    }

    shiny::observeEvent(input$rerun, {
      sec_cache(NULL)
      do_run()
    })
    shiny::observeEvent(run_key(), {
      key <- run_key()
      if (identical(key, runs$last_key)) return()
      do_run()
    })
    # Prerequisites lost (layer removed, diff cleared): the old result is
    # about something that is no longer on screen.
    shiny::observeEvent(can_run(), {
      if (!isTRUE(can_run()$ok)) {
        integration_bundle(NULL)
        runs$last_key <- NULL
      }
    }, ignoreInit = TRUE)

    # ---- what is shown ------------------------------------------------
    # Demo and live share one dispatcher.
    plot_bundle <- shiny::reactive({
      if (isTRUE(is_demo())) example_integration_bundle()
      else integration_bundle()
    })

    shown_method <- shiny::reactive({
      b <- plot_bundle()
      if (is.null(b)) method() else b$params$method %||% method()
    })

    conc_df <- shiny::reactive({
      b <- plot_bundle()
      shiny::req(b)
      b$results$integration_df
    })

    output$header <- shiny::renderUI({
      info <- can_run()
      subtitle <- if (isTRUE(is_demo())) {
        "Proteomics × RNA-seq · demo data · 60 paired features"
      } else if (isTRUE(info$ok)) {
        what <- switch(method(),
                       correlation = "sample-level correlation",
                       active_pathways = "ActivePathways",
                       gsub("_vs_", " vs ", info$comparison %||%
                              paste0(info$case, "_vs_", info$control)))
        sprintf("%s × %s · %s", info$primary_tag,
                info$secondary_tag, what)
      } else {
        "not run yet"
      }
      view_header(
        title    = "Multi-omics integration",
        subtitle = subtitle,
        actions  = shiny::actionButton(
          session$ns("rerun"),
          if (is.null(integration_bundle())) "Run integration" else "Re-run integration",
          icon = shiny::icon("play"),
          class = "btn btn-primary"
        )
      )
    })

    output$notices <- shiny::renderUI({
      tagged <- htmltools::tagList()
      err <- integration_error()
      if (identical(err, CANCELLED_MESSAGE)) {
        tagged <- htmltools::tagAppendChild(tagged, notice(
          "Cancelled", "Press Run integration to start again.", kind = "info"))
      } else if (!is.null(err)) {
        tagged <- htmltools::tagAppendChild(
          tagged,
          notice(title  = "The integration could not be computed",
                 detail = integration_error_hint(err),
                 kind   = "error",
                 technical = err)
        )
      }
      if (isTRUE(running())) {
        tagged <- htmltools::tagAppendChild(
          tagged, notice("Running…", kind = "info"))
      }
      info <- can_run()
      if (isTRUE(is_demo())) {
        tagged <- htmltools::tagAppendChild(tagged, notice(
          title  = "Showing demo data",
          detail = paste("Import two omics layers (or load the example",
                         "project from the Project view) to integrate your own."),
          kind   = "info"))
      } else if (!isTRUE(info$ok)) {
        reason <- switch(info$reason %||% "",
          layers = "Integration needs two layers in the project. Import a second omics layer.",
          diff   = paste("Fold-change concordance compares two differential results.",
                         "Run a differential analysis first, or switch to sample-level correlation."),
          info$detail %||% "The prerequisites are not met.")
        tagged <- htmltools::tagAppendChild(tagged, notice(
          title = "Nothing to integrate yet", detail = reason, kind = "info"))
        if (identical(info$reason, "diff") && is.function(navigate)) {
          tagged <- htmltools::tagAppendChild(tagged, htmltools::tags$div(
            style = "margin:6px 0 10px",
            shiny::actionButton(session$ns("go_diff"), "Go to Differential →",
                                class = "btn btn-sm btn-ghost")))
        }
      } else if (length(info$dropped_covs)) {
        tagged <- htmltools::tagAppendChild(tagged, notice(
          title = "Covariates not repeated on the partner layer",
          detail = sprintf("'%s' has no %s column, so its differential run is not adjusted for it.",
                           info$secondary_tag,
                           paste(sprintf("'%s'", info$dropped_covs), collapse = ", ")),
          kind = "warn"))
      }
      b <- integration_bundle()
      if (!is.null(b) && length(b$warnings)) {
        for (w in b$warnings) {
          tagged <- htmltools::tagAppendChild(tagged, notice(w, kind = "warn"))
        }
      }
      tagged
    })

    if (is.function(navigate)) {
      shiny::observeEvent(input$go_diff, navigate("diff"))
    }

    output$stats <- shiny::renderUI({
      b <- plot_bundle()
      shiny::req(b)
      df <- b$results$integration_df
      shiny::req(nrow(df) > 0L)
      integration_stat_cards(df, b)
    })

    output$results <- shiny::renderUI({
      m <- shown_method()
      if (!isTRUE(is_demo()) && is.null(integration_bundle())) return(NULL)
      ns <- session$ns
      if (identical(m, "correlation")) {
        htmltools::tagList(
          htmltools::tags$div(class = "row-grid r-6-6",
            integration_plot_card(ns("cor_scatter"),
                                  "Correlation per gene",
                                  "r across paired samples · y = -log10 adj. p"),
            integration_table_card(ns("top_table"), "Top features",
                                   "ranked by adjusted p")))
      } else if (identical(m, "active_pathways")) {
        htmltools::tagList(
          htmltools::tags$div(class = "row-grid r-6-6",
            integration_plot_card(ns("ap_dot"), "ActivePathways",
                                  "shape = which layers found it"),
            integration_table_card(ns("top_table"), "Pathways",
                                   "ranked by adjusted p")))
      } else {
        htmltools::tagList(
          htmltools::tags$div(class = "row-grid r-6-6",
            integration_plot_card(ns("dual"), "Mirrored volcano",
                                  "x = effect(A) - effect(B) · y = combined p"),
            integration_plot_card(ns("scatter"), "Fold-change concordance",
                                  "A vs B · dashed = y=x · coloured = hit in both")),
          htmltools::tags$div(class = "row-grid r-6-6",
            integration_table_card(ns("top_table"), "Features",
                                   "hits in both layers first"),
            if (isTRUE(is_demo())) integration_ap_card(ns)))
      }
    })

    output$dual <- shiny::renderPlot({
      b <- plot_bundle()
      shiny::req(b, identical(b$params$method, "concordance"))
      omicsCore::plot_integration(b, view = "dual_volcano")
    })

    output$scatter <- shiny::renderPlot({
      b <- plot_bundle()
      shiny::req(b, identical(b$params$method, "concordance"))
      omicsCore::plot_integration(b, view = "effect_pair")
    })

    output$cor_scatter <- shiny::renderPlot({
      b <- plot_bundle()
      shiny::req(b, identical(b$params$method, "correlation"))
      omicsCore::plot_integration(b, view = "scatter")
    })

    output$ap_dot <- shiny::renderPlot({
      b <- plot_bundle()
      shiny::req(b, identical(b$params$method, "active_pathways"))
      omicsCore::plot_integration(b, view = "dotplot")
    })

    output$top_table <- DT::renderDT({
      b <- plot_bundle()
      shiny::req(b)
      out <- integration_result_table(b$results$integration_df,
                                      b$params$method %||% "concordance",
                                      b$params$experiments)
      DT::datatable(out, rownames = FALSE, selection = "single",
                    options = list(pageLength = 10, dom = "ftip",
                                   scrollX = TRUE))
    }, server = TRUE)

    # The pathway fixture is only ever part of the demo; a live
    # ActivePathways run draws from its own result.
    ap_df <- shiny::reactive(example_integration_tables()$active_pathways_df)

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

    # Exposed for the report and the project.
    shiny::reactive(integration_bundle())
  })
}

# ---- internal helpers ------------------------------------------------

utils::globalVariables(".data")

`%||%` <- function(a, b) if (is.null(a)) b else a

# A plain-language reading of the errors a run most often ends in.
integration_error_hint <- function(msg) {
  msg <- msg %||% ""
  if (grepl("No shared .* features", msg)) {
    return(paste("The two layers have no gene symbols in common. Check that both",
                 "carry a symbol column (feature_symbol) from the same organism."))
  }
  if (grepl("guess", msg, fixed = TRUE)) {
    return("Accept the suggested sample pairing below, or add a donor column to both layers.")
  }
  if (grepl("No sample pairing|paired samples", msg)) {
    return("Sample-level correlation needs samples paired across the layers; see the pairing card below.")
  }
  if (grepl("ActivePathways", msg, fixed = TRUE)) {
    return("The ActivePathways package is not installed on this server.")
  }
  "See the technical details below."
}

integration_stat_cards <- function(df, bundle) {
  method <- bundle$params$method %||% "concordance"
  exps <- bundle$params$experiments %||% c("A", "B")
  if (identical(method, "correlation")) {
    sig <- df$is_significant %in% TRUE
    return(htmltools::tags$div(
      class = "stat-grid",
      stat_card("Genes correlated", format(nrow(df), big.mark = ","),
                trend = sprintf("%s ↔ %s", exps[1], exps[2]),
                accent = "brand", mono = TRUE),
      stat_card("Paired samples", bundle$params$method_info$n_samples %||% "—",
                trend = sprintf("pairing: %s",
                                bundle$params$method_info$pairing_source %||% "—"),
                mono = TRUE),
      stat_card("Positive (sig)", sum(sig & df$effect > 0, na.rm = TRUE),
                trend = "RNA and protein rise together", accent = "up", mono = TRUE),
      stat_card("Median r", sprintf("%.2f", stats::median(df$effect, na.rm = TRUE)),
                trend = "all genes", mono = TRUE)
    ))
  }
  if (identical(method, "active_pathways")) {
    sig <- df$is_significant %in% TRUE
    return(htmltools::tags$div(
      class = "stat-grid",
      stat_card("Pathways tested", nrow(df), accent = "brand", mono = TRUE),
      stat_card("Significant", sum(sig), mono = TRUE),
      stat_card("Shared by both", sum(sig & df$direction %in% "shared"),
                accent = "up", mono = TRUE),
      stat_card("Combined only", sum(sig & df$direction %in% "combined"),
                trend = "found only when the layers are merged", mono = TRUE)
    ))
  }
  # Concordance. The counts are of features that are hits in *both*
  # layers: every feature has some sign pair, so counting quadrants over
  # all of them reported half of an unrelated background as "concordant".
  both <- if (all(c("significant_a", "significant_b") %in% names(df))) {
    df$significant_a %in% TRUE & df$significant_b %in% TRUE
  } else {
    df$is_significant %in% TRUE
  }
  quad <- df$quadrant
  up_n <- sum(both & quad %in% "up_up")
  down_n <- sum(both & quad %in% "down_down")
  disc_n <- sum(both & quad %in% c("up_down", "down_up"))
  rho <- if (all(c("effect_a", "effect_b") %in% names(df)) && sum(both) >= 3L) {
    suppressWarnings(stats::cor(df$effect_a[both], df$effect_b[both],
                                method = "spearman", use = "pairwise.complete.obs"))
  } else NA_real_
  htmltools::tags$div(
    class = "stat-grid",
    stat_card(
      label  = "Paired features",
      value  = format(nrow(df), big.mark = ","),
      trend  = sprintf("%s ↔ %s · %d hits in both", exps[1], exps[2], sum(both)),
      accent = "brand", mono = TRUE
    ),
    stat_card(label = "Concordant ↑", value = up_n,
              trend = "hit and up in both layers", accent = "up", mono = TRUE),
    stat_card(label = "Concordant ↓", value = down_n,
              trend = "hit and down in both layers", accent = "down", mono = TRUE),
    stat_card(
      label = "Discordant",
      value = disc_n,
      trend = if (is.na(rho)) "hit in both, opposite signs"
              else sprintf("opposite signs · Spearman ρ of hits %.2f", rho),
      accent = if (disc_n > 0L) "warn" else "ok",
      mono = TRUE
    )
  )
}

integration_result_table <- function(df, method, experiments) {
  exps <- experiments %||% c("A", "B")
  if (identical(method, "concordance") &&
      all(c("effect_a", "effect_b") %in% names(df))) {
    both <- if (all(c("significant_a", "significant_b") %in% names(df))) {
      df$significant_a %in% TRUE & df$significant_b %in% TRUE
    } else df$is_significant %in% TRUE
    ord <- order(!both, df$p_value, na.last = TRUE)
    d <- df[ord, , drop = FALSE]
    out <- data.frame(
      Feature = d$feature_symbol,
      a = round(d$effect_a, 3),
      b = round(d$effect_b, 3),
      Quadrant = d$quadrant,
      `Hit in both` = ifelse(both[ord], "yes", ""),
      `Combined p` = signif(d$p_value, 3),
      check.names = FALSE, stringsAsFactors = FALSE)
    names(out)[2:3] <- paste0("Effect (", exps[1:2], ")")
    return(out)
  }
  d <- df[order(df$adj_p_value, na.last = TRUE), , drop = FALSE]
  data.frame(
    Feature = d$feature_symbol,
    Effect = round(d$effect, 3),
    `Adj. p` = signif(d$adj_p_value, 3),
    Direction = d$direction,
    check.names = FALSE, stringsAsFactors = FALSE)
}

integration_setup_card <- function(ns) {
  choices <- c("Fold-change concordance" = "concordance",
               "Sample-level correlation" = "correlation")
  if (has_pkg("ActivePathways")) {
    choices <- c(choices, "ActivePathways (pathways)" = "active_pathways")
  }
  bslib::card(
    bslib::card_header(
      htmltools::tags$h3(class = "card-title", "Setup"),
      htmltools::tags$span(class = "card-sub", "which layers, and how")
    ),
    bslib::card_body(
      htmltools::tags$div(
        class = "row-grid r-6-6",
        htmltools::tags$div(
          htmltools::tags$label(class = "control-label", "Layers"),
          shiny::uiOutput(ns("ui_partner"))
        ),
        shiny::radioButtons(
          ns("method"), label = "Method", choices = choices,
          selected = "concordance", inline = TRUE
        )
      ),
      htmltools::tags$div(
        class = "muted", style = "font-size:12px",
        paste("Concordance repeats the Differential view's contrast on the",
              "second layer and compares the two results gene by gene, at the",
              "same thresholds. Correlation compares the two layers sample by",
              "sample and needs the pairing below.")
      )
    )
  )
}

integration_pairing_card <- function(ns) {
  bslib::card(
    bslib::card_header(
      htmltools::tags$h3(class = "card-title", "Sample pairing"),
      htmltools::tags$span(class = "card-sub",
                           "which sample is which person"),
      shiny::uiOutput(ns("pairing_action"), inline = TRUE)
    ),
    bslib::card_body(
      shiny::uiOutput(ns("pairing_note")),
      # Folded: the note says how many pairs and where they came from,
      # which is what most readers need; the rows are one click away for
      # the reader who wants to check them.
      htmltools::tags$details(
        class = "pairing-details",
        htmltools::tags$summary("Show the pairs"),
        DT::DTOutput(ns("pairing_table"))
      )
    )
  )
}

integration_plot_card <- function(output_id, title, sub) {
  bslib::card(
    bslib::card_header(
      htmltools::tags$h3(class = "card-title", title),
      htmltools::tags$span(class = "card-sub", sub)
    ),
    bslib::card_body(shiny::plotOutput(output_id, height = "320px"))
  )
}

integration_table_card <- function(output_id, title, sub) {
  bslib::card(
    bslib::card_header(
      htmltools::tags$h3(class = "card-title", title),
      htmltools::tags$span(class = "card-sub", sub)
    ),
    bslib::card_body(DT::DTOutput(output_id))
  )
}

integration_ap_card <- function(ns) {
  bslib::card(
    bslib::card_header(
      htmltools::tags$h3(class = "card-title",
                         "ActivePathways · combined p"),
      htmltools::tags$span(class = "card-sub",
                           "Brown's method · demo data")
    ),
    bslib::card_body(
      DT::DTOutput(ns("ap_table"))
    )
  )
}
