# What Re-run starts in the differential view.
#
# Split out of mod_diff_view.R and called from inside its module server,
# so the ids and the module's state stay where they were. It defines no
# outputs or observers of its own: it returns do_run() and
# run_continuous(), which the Re-run observer calls.
diff_run_server <- function(input, output, session, active, default_contrast,
                            design_mode, continuous_cols, settings_now, set_busy,
                            diff_epoch, diff_error, ran_with, diff_bundle) {
  do_run <- function() {
    a <- active()
    snap <- shiny::isolate(settings_now())
    if (is.null(a$input)) {
      diff_error("This project has no layers to analyse.")
      return(invisible())
    }
    d <- default_contrast()
    method     <- input$method %||% "auto"
    # Fall back to the derived default, not to nothing: on the first
    # run the controls have not been rendered yet, and refusing then
    # put an error on a view the user had only just opened.
    group_col  <- input$group_col %||% d$group_col
    control    <- input$control   %||% d$control
    case       <- input$case      %||% d$case
    case       <- setdiff(case, control)
    covariates <- input$covariates
    mode       <- input$contrast_mode %||% "control"
    contrasts  <- switch(
      mode,
      pairwise = "pairwise",
      custom   = {
        lines <- trimws(strsplit(input$custom_contrasts %||% "", "\n")[[1L]])
        lines[nzchar(lines)]
      },
      NULL)
    if (identical(design_mode(), "continuous")) {
      return(run_continuous(a, method, covariates))
    }
    paired_col <- input$paired_col
    if (!length(paired_col) || !nzchar(paired_col)) paired_col <- NULL
    if (identical(mode, "custom") && !length(contrasts)) {
      diff_error("Write at least one comparison, e.g. \"TreatB - TreatA\".")
      return(invisible())
    }
    if (is.null(group_col) && !length(d$candidates)) {
      diff_error(paste("No column of the sample metadata splits the samples",
                       "into groups. Add one (e.g. 'group' = Control / Treated)",
                       "to the sample sheet and import again."))
      return(invisible())
    }
    if (is.null(contrasts) &&
        (is.null(group_col) || is.null(control) || !length(case))) {
      diff_error("Pick a group column, a control group, and at least one group distinct from the control to compare with it.")
      return(invisible())
    }
    set_busy(TRUE)
    my_run <- diff_epoch$start()
    run_async(
      # Detached, so the worker receives the input and the six
      # parameters rather than this module's whole scope. Defined
      # inline it carried the previous bundle, the project and the
      # demo fixtures with it -- 527 MB of "globals", which future
      # refused to export, from inside an observer, which greyed the
      # page.
      detached_call(
        function() {
          omicsCore::run_diff(
            input         = inp,
            method        = method,
            analysis_type = "group",
            group_col     = group_col,
            control_group = control,
            case_group    = if (is.null(contrasts)) case,
            covariates    = covariates,
            paired_col    = paired_col,
            contrasts     = contrasts
          )
        },
        inp        = a$input,
        method     = method,
        group_col  = group_col,
        control    = control,
        case       = case,
        covariates = if (length(covariates)) covariates else NULL,
        paired_col = paired_col,
        contrasts  = contrasts
      ),
      on_success = function(bundle) {
        if (diff_epoch$is_last_started(my_run)) set_busy(FALSE)
        if (!diff_epoch$is_current(my_run)) return(invisible())
        diff_error(NULL)
        ran_with(snap)
        diff_bundle(with_layer(bundle, a))
      },
      on_error = function(msg) {
        if (diff_epoch$is_last_started(my_run)) set_busy(FALSE)
        if (!diff_epoch$is_current(my_run)) return(invisible())
        # The previous result stays, labelled as the previous one: a
        # mistyped contrast should not lose it, and unlabelled beside
        # the error it read as the answer to the new settings.
        diff_error(msg)
      },
      message = "Running differential analysis..."
    )
  }

  run_continuous <- function(a, method, covariates) {
    snap <- shiny::isolate(settings_now())
    col <- input$continuous_col %||% continuous_cols()[1L]
    if (is.null(col) || is.na(col)) {
      diff_error("Pick a numeric column to test for a trend.")
      return(invisible())
    }
    model <- input$continuous_model %||% "linear"
    set_busy(TRUE)
    my_run <- diff_epoch$start()
    run_async(
      detached_call(
        function() {
          args <- list(input = inp, method = method, analysis_type = "continuous",
                       continuous_col = col, covariates = covariates)
          # The spline is limma's; a linear trend runs on any engine.
          if (identical(model, "spline")) {
            args$method <- "limma"
            args$model <- "spline"
          }
          do.call(omicsCore::run_diff, args)
        },
        inp = a$input, method = method, col = col, model = model,
        covariates = if (length(covariates)) covariates else NULL
      ),
      on_success = function(bundle) {
        if (diff_epoch$is_last_started(my_run)) set_busy(FALSE)
        if (!diff_epoch$is_current(my_run)) return(invisible())
        diff_error(NULL)
        ran_with(snap)
        diff_bundle(with_layer(bundle, a))
      },
      on_error = function(msg) {
        if (diff_epoch$is_last_started(my_run)) set_busy(FALSE)
        if (!diff_epoch$is_current(my_run)) return(invisible())
        diff_error(msg)
      },
      message = "Testing for trends..."
    )
  }

  list(do_run = do_run, run_continuous = run_continuous)
}
