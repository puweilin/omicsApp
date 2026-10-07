# The Parameters card of the differential view: its UI, and the server
# side that fills it.
#
# Split out of mod_diff_view.R. The server half is a plain function
# called from inside diff_view_server()'s moduleServer(), not a module of
# its own, so every input and output keeps the id it has always had
# ("diff-layer", "diff-ui_contrast"); the reactives it builds are handed
# back and kept under the same names in the module.

# Builds the layer, method, contrast and covariate controls, and returns
# the reactives the rest of the view reads: the active layer, the default
# contrast, the design mode, and the levels of the chosen group column.
diff_params_server <- function(input, output, session, current_project, preferred_layer) {
  # Active experiment: whichever layer the user picked, else
  # proteomics, else the first of any kind.
  #
  # The layer is a control rather than a fixed choice because the
  # engines available depend on it: applicable_diff_methods() offers
  # deseq2 and edger only for rnaseq raw counts. Pinned to proteomics,
  # as this was, those two could never be reached from the interface
  # at all -- the gate was right and there was no way to the side of
  # it where it opens.
  #
  # The demo resolves through example_project() for the same reason,
  # rather than being handed a fixed proteomics input: its rnaseq
  # layer is raw counts, so it is the one place a user can see the
  # method list change without importing anything.
  active <- shiny::reactive({
    proj <- current_project()
    is_demo <- is.null(proj)
    if (is_demo) proj <- example_project()
    exps <- proj$experiments
    if (length(exps) == 0L) {
      return(list(input = NULL, tag = NULL, is_demo = TRUE))
    }
    want <- input$layer %||% preferred_layer()
    tag <- if (!is.null(want) && want %in% names(exps)) {
      want
    } else {
      default_layer_tag(exps)
    }
    list(input = exps[[tag]], tag = tag, is_demo = is_demo)
  })

  output$ui_layer <- shiny::renderUI({
    proj <- current_project() %||% example_project()
    tags_avail <- names(proj$experiments)
    if (length(tags_avail) == 0L) return(NULL)
    # isolate(), because active() reads input$layer and this output
    # writes it. Reading it here closed the loop: re-rendering the
    # control re-sent its value, which invalidated active(), which
    # re-rendered the control.
    sel <- shiny::isolate(input$layer) %||% shiny::isolate(preferred_layer())
    if (is.null(sel) || !sel %in% tags_avail) {
      sel <- default_layer_tag(proj$experiments)
    }
    shiny::selectInput(
      session$ns("layer"), label = "Omics layer",
      choices = tags_avail, selected = sel
    )
  })

  # Method dropdown, restricted to the engines whose assumptions the
  # active layer meets. DESeq2 handed continuous intensities does not
  # error: it rounds them to integers and reports p-values for a
  # negative-binomial model the data never fitted. Nothing downstream
  # can tell that apart from a real result, so the guard has to be
  # here, at the point of choosing.
  output$ui_method <- shiny::renderUI({
    a <- active()
    inp <- a$input
    shiny::req(inp)   # a project with no layers has no methods to offer
    choices <- omicsCore::applicable_diff_methods(inp)
    # Keep the user's pick. active() changes whenever the project does
    # -- including when a finished run is attached to it -- so the
    # control re-rendered after every run and snapped back to "auto".
    sel <- shiny::isolate(input$method)
    if (is.null(sel) || !sel %in% choices) sel <- "auto"
    shiny::selectInput(session$ns("method"),
                       label = htmltools::tags$span(class = "visually-hidden-label", "Method"),
                       choices = choices, selected = sel)
  })

  output$method_note <- shiny::renderUI({
    a <- active()
    inp <- a$input
    shiny::req(inp)
    dropped <- setdiff(omicsCore::SUPPORTED_DIFF_METHODS,
                       omicsCore::applicable_diff_methods(inp))
    if (length(dropped) == 0L) return(NULL)
    htmltools::tags$div(
      class = "muted", style = "font-size:11.5px;padding-top:4px",
      sprintf("%s hidden: not valid for %s data.",
              paste(dropped, collapse = ", "),
              inp$assay_type %||% inp$omics_type)
    )
  })

  # The contrast the controls will default to, computed from the data
  # rather than read back off them.
  #
  # do_run() fires once as soon as active() settles, and at that
  # moment the controls do not exist yet: they are renderUI output,
  # and Shiny has not been round the loop. Reading input$group_col
  # there gets NULL, which used to surface as "Pick a group column
  # with distinct Control and Case levels" on a view the user had
  # only just opened. Deriving the default in one place means the
  # first run and the rendered controls cannot disagree about it.
  default_contrast <- shiny::reactive({
    meta <- active()$input$meta_df
    if (is.null(meta) || !ncol(meta)) {
      return(list(group_col = NULL, control = NULL, case = NULL,
                  levels = character(0), candidates = character(0)))
    }
    cands <- grouping_candidates(meta)
    # No column with two replicated groups. Offer the columns that at
    # least repeat a value -- never one that names each sample once,
    # which is the sample ID and cannot be compared.
    if (length(cands) == 0L) cands <- grouping_candidates(meta, min_per_level = 1L,
                                                          replicated = TRUE)
    # What was stated at import comes first; the guess is the fallback.
    design <- tryCatch(omicsCore::study_design(active()$input), error = function(e) NULL)
    if (!is.null(design)) cands <- c(design$group_col, setdiff(cands, design$group_col))
    # NULL, not cands[1L]: that is NA on an empty vector, and the
    # "no column splits the samples" message below never fired.
    gc <- if (length(cands)) cands[1L] else NULL
    if (is.null(gc)) {
      return(list(group_col = NULL, control = NULL, case = NULL,
                  levels = character(0), candidates = character(0)))
    }
    lv <- sort(unique(as.character(stats::na.omit(meta[[gc]]))))
    ctrl <- if (!is.null(design$reference) && identical(gc, design$group_col)) {
      design$reference
    } else default_control_level(lv)
    n_per <- table(as.character(meta[[gc]]))
    list(
      group_col  = gc,
      control    = ctrl,
      # Every other group against the control: with one control and
      # several treatments that is the design people mean, and with two
      # groups it is the one contrast there is. A group of one sample
      # cannot be tested and is left out of the default.
      case       = if (length(lv) >= 2L) setdiff(lv[n_per[lv] >= 2L], ctrl) else NULL,
      levels     = lv,
      candidates = cands
    )
  })

  # Groups or a continuous variable (dose, time, age). Numeric columns
  # with many values had no way in: they are not groups, and the view
  # offered nothing else.
  continuous_cols <- shiny::reactive({
    meta <- active()$input$meta_df
    if (is.null(meta)) return(character(0))
    names(meta)[vapply(meta, function(col) {
      v <- suppressWarnings(as.numeric(as.character(col)))
      sum(!is.na(v)) >= 4L && length(unique(stats::na.omit(v))) >= 3L
    }, logical(1))]
  })
  design_mode <- shiny::reactive(input$design_mode %||% "groups")
  # A trend's effect is a slope, not a fold change; the cutoff says so.
  shiny::observeEvent(design_mode(), {
    shiny::updateNumericInput(
      session, "fc_cut",
      label = if (identical(design_mode(), "continuous")) "|slope| cutoff" else "|log2FC| cutoff")
  }, ignoreInit = TRUE)

  output$ui_design_mode <- shiny::renderUI({
    if (!length(continuous_cols())) return(NULL)
    shiny::radioButtons(session$ns("design_mode"), label = "Compare",
                        choices = c("groups" = "groups",
                                    "a continuous variable" = "continuous"),
                        selected = shiny::isolate(input$design_mode) %||% "groups",
                        inline = TRUE)
  })

  output$ui_continuous <- shiny::renderUI({
    if (!identical(design_mode(), "continuous")) return(NULL)
    cols <- continuous_cols()
    sel <- shiny::isolate(input$continuous_col)
    if (is.null(sel) || !sel %in% cols) sel <- cols[[1L]]
    htmltools::tagList(
      shiny::selectInput(session$ns("continuous_col"), "Variable (dose, time, age...)",
                         choices = cols, selected = sel),
      shiny::radioButtons(session$ns("continuous_model"), label = "Model",
                          choices = c("linear trend" = "linear",
                                      "any smooth change (spline)" = "spline"),
                          selected = shiny::isolate(input$continuous_model) %||% "linear",
                          inline = TRUE),
      htmltools::tags$div(
        class = "muted", style = "font-size:11.5px;margin-top:-6px",
        paste("Each feature is tested for a trend with the variable. The",
              "effect shown is Spearman's rho; a spline asks whether it",
              "changes at all, without a direction."))
    )
  })

  # Columns that pair samples across the groups: one sample per group
  # for each patient, donor, mouse. Offered whenever such a column
  # exists, and chosen for the user when its name says it is one --
  # a before/after design analysed unpaired found nothing at all.
  pairing_candidates <- shiny::reactive({
    meta <- active()$input$meta_df
    gc <- input$group_col %||% default_contrast()$group_col
    if (is.null(meta) || is.null(gc) || !gc %in% names(meta)) return(character(0))
    g <- as.character(meta[[gc]])
    names(meta)[vapply(names(meta), function(nm) {
      if (identical(nm, gc)) return(FALSE)
      b <- as.character(meta[[nm]])
      ok <- !is.na(b) & !is.na(g)
      if (sum(ok) < 4L) return(FALSE)
      per_block <- table(b[ok])
      nb <- length(per_block)
      # Repeated, but not one level per sample and not the group again.
      nb >= 2L && nb < sum(ok) && all(per_block >= 2L) &&
        # The same block appears in more than one group.
        mean(tapply(g[ok], b[ok], function(x) length(unique(x))) >= 2L) >= 0.5
    }, logical(1))]
  })

  output$ui_paired <- shiny::renderUI({
    if (identical(design_mode(), "continuous")) return(NULL)
    cands <- pairing_candidates()
    if (!length(cands)) return(NULL)
    sel <- shiny::isolate(input$paired_col)
    if (is.null(sel) || !sel %in% c("", cands)) {
      hinted <- cands[grepl(PAIRING_COL_HINTS, tolower(cands))]
      sel <- if (length(hinted)) hinted[[1L]] else ""
    }
    htmltools::tagList(
      shiny::selectInput(session$ns("paired_col"),
                         label = "Paired by (same patient / donor / animal)",
                         choices = c("not paired" = "", cands), selected = sel),
      htmltools::tags$div(
        class = "muted", style = "font-size:11.5px;margin-top:-6px",
        paste("Each pair is compared with itself, which removes the",
              "differences between individuals."))
    )
  })

  # Group column dropdown: any meta_df column with >= 2 unique
  # non-NA values (continuous columns like `age` are excluded
  # for the simple "control vs case" UI in this slice).
  output$ui_group_col <- shiny::renderUI({
    d <- default_contrast()
    if (!length(d$candidates)) return(NULL)
    # Keep what the user picked. This output re-renders whenever the
    # layer or the project changes, and re-rendering with the default
    # threw away their choice -- they selected condition / G1 / G2,
    # ran it, and the control came back reading `label`.
    sel <- shiny::isolate(input$group_col)
    if (is.null(sel) || !sel %in% d$candidates) sel <- d$group_col
    shiny::selectInput(session$ns("group_col"),
                       label    = "Group column",
                       choices  = d$candidates,
                       selected = sel)
  })

  # Reactive level set for the chosen group column.
  levels_ <- shiny::reactive({
    a <- active()
    meta <- a$input$meta_df
    gc <- input$group_col %||% default_contrast()$group_col
    if (is.null(gc) || !(gc %in% names(meta))) return(character(0))
    sort(unique(as.character(stats::na.omit(meta[[gc]]))))
  })

  output$ui_contrast <- shiny::renderUI({
    if (identical(design_mode(), "continuous")) return(NULL)
    lv <- levels_()
    if (length(lv) < 2L) {
      return(htmltools::tags$div(
        class = "muted",
        style = "font-size:12px",
        "Pick a group column with at least two levels."
      ))
    }
    keep <- function(current, fallback) {
      current <- shiny::isolate(current)
      current <- current[current %in% lv]
      if (!length(current)) fallback else current
    }
    d <- default_contrast()
    meta <- active()$input$meta_df
    gc_now <- input$group_col %||% d$group_col
    n_per <- table(as.character(meta[[gc_now]]))
    small <- lv[n_per[lv] < 3L]
    ctrl <- keep(input$control,
                 if (!is.null(d$control) && d$control %in% lv &&
                     identical(input$group_col %||% d$group_col, d$group_col))
                   d$control else default_control_level(lv))[1L]
    mode <- shiny::isolate(input$contrast_mode) %||% "control"
    htmltools::tagList(
      # Three designs, one fit each: every treatment against the control
      # (the common case), every pair of groups, or comparisons written
      # out -- "(TreatA + TreatB)/2 - Control", "TreatB - TreatA".
      if (length(lv) > 2L) {
        shiny::radioButtons(
          session$ns("contrast_mode"), label = "Compare",
          choices = c("each vs control" = "control",
                      "all pairs" = "pairwise",
                      "custom" = "custom"),
          selected = mode, inline = TRUE)
      },
      shiny::selectInput(session$ns("control"),
                         label = "Control (reference)", choices = lv,
                         selected = ctrl),
      shiny::conditionalPanel(
        sprintf("!input['%s'] || input['%s'] == 'control'",
                session$ns("contrast_mode"), session$ns("contrast_mode")),
        shiny::selectizeInput(
          session$ns("case"),
          label = if (length(lv) > 2L) "Compare against control (one or more)" else "Case",
          # The control is not something to compare with itself.
          choices = setdiff(lv, ctrl), multiple = TRUE,
          selected = setdiff(keep(input$case, setdiff(lv[n_per[lv] >= 2L], ctrl)), ctrl),
          options = list(plugins = list("remove_button")))
      ),
      if (length(small)) {
        htmltools::tags$div(
          class = "muted", style = "font-size:11.5px;color:var(--warn)",
          sprintf("Few samples: %s. A group needs at least 2 samples to be tested and 3 for a reliable estimate.",
                  paste(sprintf("%s (n = %d)", small, as.integer(n_per[small])), collapse = ", ")))
      },
      if (length(lv) > 2L) shiny::conditionalPanel(
        sprintf("input['%s'] == 'custom'", session$ns("contrast_mode")),
        shiny::textAreaInput(
          session$ns("custom_contrasts"),
          label = "Comparisons, one per line",
          value = shiny::isolate(input$custom_contrasts) %||% "",
          placeholder = paste(
            c(paste(contrast_token(lv[3L]), "-", contrast_token(lv[2L])),
              sprintf("(%s + %s)/2 - %s", contrast_token(lv[2L]),
                      contrast_token(lv[3L]), contrast_token(lv[1L]))),
            collapse = "\n"),
          rows = 3),
        htmltools::tags$div(
          class = "muted", style = "font-size:11.5px;margin-top:-6px",
          "Group names that are not plain words go in backticks: `Drug A` - Control.")
      ),
      if (length(lv) > 2L) {
        htmltools::tags$div(
          class = "muted", style = "font-size:11.5px;margin-top:-6px",
          paste("Every group compared is fitted in one model, so each",
                "comparison borrows strength from all of their samples."))
      }
    )
  })

  # The control is not a case: keep the case list free of it when the
  # reference changes.
  shiny::observeEvent(input$control, {
    lv <- levels_()
    shiny::updateSelectizeInput(session, "case", choices = setdiff(lv, input$control),
                                selected = setdiff(input$case, input$control))
  }, ignoreInit = TRUE)

  output$ui_covariates <- shiny::renderUI({
    a <- active()
    meta <- a$input$meta_df
    if (identical(input$method, "ttest")) {
      return(htmltools::tags$div(class = "muted", style = "font-size:12px",
                                 "The t-test does not adjust for covariates; choose limma or lm to add them."))
    }
    gc <- if (identical(design_mode(), "continuous")) input$continuous_col %||% ""
          else input$group_col %||% ""
    cands <- setdiff(names(meta), c(gc, "sample_id", input$paired_col))
    # Kept across re-renders, for the reason given at ui_method.
    sel <- intersect(shiny::isolate(input$covariates), cands)
    shiny::selectizeInput(
      session$ns("covariates"),
      label    = NULL,
      choices  = cands,
      multiple = TRUE,
      selected = if (length(sel)) sel,
      options  = list(placeholder = "optional, e.g. age")
    )
  })

  list(active = active, default_contrast = default_contrast,
       continuous_cols = continuous_cols, design_mode = design_mode,
       pairing_candidates = pairing_candidates, levels_ = levels_)
}

diff_params_card <- function(ns) {
  bslib::card(
    class = "params-card",
    bslib::card_header(
      htmltools::tags$h3(class = "card-title", "Parameters"),
      htmltools::tags$span(class = "card-sub",
                           "design + thresholds")
    ),
    bslib::card_body(
      htmltools::tags$div(
        class = "param-stack",
        param_group(
          "Layer",
          # Which engines are on offer follows from this: deseq2 and
          # edger need rnaseq raw counts, so without a way to change
          # layer they were unreachable.
          shiny::uiOutput(ns("ui_layer"))
        ),
        param_group(
          "Method",
          help = paste("auto picks limma for proteomics and log-scale RNA-seq,",
                       "DESeq2 for raw counts. Engines that do not fit the",
                       "layer's data are hidden."),
          # Rendered server-side: which engines are valid depends on the
          # active layer's assay, and offering an invalid one produces a
          # complete, plausible, meaningless result table.
          shiny::uiOutput(ns("ui_method")),
          shiny::uiOutput(ns("method_note"))
        ),
        param_group(
          "Contrast",
          help = paste("The control is the reference every other group is compared",
                       "with. Choose several groups to fit them in one model;",
                       "'Showing' then picks the comparison on screen."),
          shiny::uiOutput(ns("ui_design_mode")),
          shiny::uiOutput(ns("ui_continuous")),
          shiny::uiOutput(ns("ui_group_col")),
          shiny::uiOutput(ns("ui_paired")),
          shiny::uiOutput(ns("ui_contrast")),
          shiny::uiOutput(ns("ui_comparison"))
        ),
        param_group(
          "Covariates",
          help = paste("Variables to adjust for (age, sex, batch). The group effect",
                       "is then estimated holding them constant. Leave empty for",
                       "an unadjusted comparison. For before/after or matched",
                       "samples use 'Paired by' above instead."),
          shiny::uiOutput(ns("ui_covariates"))
        ),
        param_group(
          "Thresholds",
          help = paste("Which features count as hits. They filter the table and the",
                       "counts, and are what Enrichment and Integration use; they",
                       "do not re-run the model."),
          # Which p to threshold on is the user's call, not ours. An
          # exploratory screen on 50 proteins and a confirmatory one on
          # 20,000 genes want different answers, and forcing adj.P made
          # the first look empty.
          shiny::radioButtons(
            ns("p_kind"), label = "Significance on",
            choices = c("adjusted p" = "adj", "raw p" = "raw"),
            selected = "adj", inline = TRUE
          ),
          # Ticks off: on a rail this narrow their labels ran together
          # ("0.020.04"); the handle shows the value.
          shiny::sliderInput(
            ns("fdr_cut"), label = "p cutoff",
            min = 0, max = 0.2, value = 0.05, step = 0.005, ticks = FALSE
          ),
          # A box, not a slider. The slider stepped 0.05, which cannot
          # express log2(1.2) = 0.263 -- so the fold change most often
          # wanted here was one of the few the control could not reach.
          shiny::numericInput(
            ns("fc_cut"), label = "|log2FC| cutoff",
            value = round(log2(1.2), 3), min = 0, max = 10, step = 0.05
          ),
          htmltools::tags$div(
            class = "muted", style = "font-size:11.5px;margin-top:-6px",
            sprintf("%.3f = %.2gx fold change", log2(1.2), 1.2)
          ),
          shinyWidgets::materialSwitch(
            ns("label_top"), label = "Label top 20",
            value = FALSE, status = "primary", right = TRUE
          )
        ),
        # Pinned to the bottom of the window while the rail scrolls: at
        # the end of the parameters it sat at y = 1,470 on a 900 px
        # screen, below the fold, under the controls it acts on.
        htmltools::tags$div(
          class = "run-sticky",
          shiny::uiOutput(ns("run_button"))
        )
      )
    )
  )
}

# A group name as a contrast expression needs it.
contrast_token <- function(x) {
  x <- as.character(x %||% "B")
  if (identical(make.names(x), x)) x else paste0("`", x, "`")
}

# The level a control group is usually called, so the default reference
# is the control rather than whichever label sorts first ("DMSO" before
# "Drug", but "Treated" after "Control" only by luck).
CONTROL_LEVEL_PATTERN <- paste0(
  "^(ctrl|control|controls|con|dmso|vehicle|veh|mock|sham|untreated|",
  "placebo|baseline|wt|wild[ _-]?type|normal|healthy|nc|neg|negative|",
  "blank|t0|0h|day0|d0|before|pre|pretreatment|pre[ _-]?treatment|",
  "pbs|saline|sinc|sictrl|shctrl|scramble|scrambled|",
  "\u5bf9\u7167|\u5bf9\u7167\u7ec4|\u7a7a\u767d|\u7a7a\u767d\u7ec4|",
  "\u6b63\u5e38|\u6b63\u5e38\u7ec4|\u9634\u6027\u5bf9\u7167)$")

default_control_level <- function(levels) {
  if (!length(levels)) return(NULL)
  hit <- grepl(CONTROL_LEVEL_PATTERN, trimws(tolower(levels)))
  if (!any(hit)) {
    hit <- grepl("(^|[ _-])(ctrl|control|vehicle|dmso|wt|before|pre|baseline)([ _-]|$)",
                 tolower(levels)) | grepl("\u5bf9\u7167", levels)
  }
  levels[if (any(hit)) which(hit)[1L] else 1L]
}

# Which layer the view lands on when the user has not chosen one:
# proteomics if present, else the first.
default_layer_tag <- function(experiments) {
  if (!length(experiments)) return(NULL)
  types <- vapply(experiments, function(e) e$omics_type %||% "", character(1))
  i <- which(types == "proteomics")
  names(experiments)[if (length(i)) i[1L] else 1L]
}

# Columns of `meta_df` that could name a contrast, best first.
#
# The old rule was "not numeric, and at least two distinct values",
# which a sample identifier satisfies perfectly: one level per sample.
# On a real workbook the first such column was `label`, whose values are
# the sample names, so the view defaulted to a contrast of one sample
# against one other. limma cannot fit that -- no residual degrees of
# freedom -- and reported it as "Partial NA coefficients for 2294
# probe(s)" and an empty result, which is not a sentence anyone can act
# on.
#
# The real requirement is replication: every level needs at least two
# samples, or the level cannot be tested. That single condition
# excludes identifiers exactly, without having to guess from names.
GROUP_COL_HINTS <- c("group", "condition", "treatment", "arm", "status",
                     "genotype", "diet", "disease", "phenotype", "diagnosis",
                     "cohort", "timepoint", "time_point", "visit",
                     "\u5206\u7ec4", "\u7ec4\u522b", "\u7ec4", "\u5904\u7406")
# Columns that describe samples but are rarely the comparison: ranked last.
NUISANCE_COL_RE <- paste0("^(batch|sex|gender|cage|plate|lane|run|replicate|rep|",
                          "patient|donor|subject|individual|animal|mouse|pair|",
                          "\u6279\u6b21|\u6027\u522b)([ _.-]?(id|no))?$")
# Columns that pair samples, preselected in 'Paired by'.
PAIRING_COL_HINTS <- paste0("^(patient|donor|subject|individual|animal|mouse|pair|",
                            "participant|case)([ _.-]?(id|no))?$|^(\u60a3\u8005|\u4e2a\u4f53)")

grouping_candidates <- function(meta, min_per_level = 2L, replicated = FALSE) {
  if (is.null(meta) || !ncol(meta)) return(character(0))
  usable <- vapply(names(meta), function(nm) {
    col <- meta[[nm]]
    # A numeric column is a grouping when it has a few repeated values:
    # 0/1 coding, dose levels. Age or BMI has as many values as samples.
    if (is.numeric(col) && length(unique(stats::na.omit(col))) > 6L) return(FALSE)
    counts <- table(as.character(col), useNA = "no")
    length(counts) >= 2L && min(counts) >= min_per_level &&
      (!replicated || max(counts) >= 2L)
  }, logical(1))
  cands <- names(meta)[usable]
  if (!length(cands)) return(character(0))

  # Fewest levels first: a two-level column is the contrast someone
  # almost always means. Conventional names win over the count, since a
  # column called `condition` is a stated intent and a level count is an
  # inference.
  n_levels <- vapply(cands, function(nm) {
    length(unique(stats::na.omit(as.character(meta[[nm]]))))
  }, integer(1))
  low <- tolower(cands)
  hinted <- low %in% GROUP_COL_HINTS | grepl("group|condition|treat", low)
  nuisance <- grepl(NUISANCE_COL_RE, low)
  cands[order(!hinted, nuisance, n_levels)]
}
