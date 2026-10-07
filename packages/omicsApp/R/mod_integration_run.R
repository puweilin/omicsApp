# The prerequisites and the run of the integration view.
#
# Split out of mod_integration_view.R. These are plain functions called
# from inside integration_view_server()'s moduleServer(), not modules, so
# the ids and the module's state (the result, the partner's cached diff,
# the last run's key) stay where they were.

# Everything a run needs, or the one reason it cannot happen, in words
# that name the missing piece. Returns the reactive.
integration_can_run <- function(current_project, layers, method, diff_bundle) {
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
    gc_primary <- params$group_col
    # The partner's own name for the column: "Group" on one layer and
    # "group" on the other blocked the run, although the partner's
    # recorded design said which column it was.
    gc <- match_partner_col(gc_primary, sec)
    if (is.null(gc)) gc <- gc_primary
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
        primary_omics  = bundle$input_info$omics_type,
        paired_col     = match_partner_col(params$paired_col, sec)
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
      primary_omics  = bundle$input_info$omics_type,
      paired_col     = match_partner_col(params$paired_col, sec)
    ), base)
  })
  can_run
}

# The key a result depends on, and do_run(), which the module's observers
# call. Defines no outputs or observers of its own.
integration_run_server <- function(input, output, session, current_project, diff_bundle,
                                   diff_thresholds, method, can_run, runs, sec_cache,
                                   set_busy, integration_bundle, integration_error) {
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
    fl <- proj$feature_link
    paste(method(), info$primary_tag, fp(info$primary_tag),
          info$secondary_tag, fp(info$secondary_tag),
          if (!identical(method(), "correlation") && !is.null(b))
            paste(b$params$method, paste(b$params$comparison, collapse = ","),
                  nrow(b$results$diff_result_df),
                  sum(b$results$diff_result_df$p_value, na.rm = TRUE)),
          th$p_cutoff, th$p_preference, th$effect_cutoff,
          if (identical(method(), "correlation") && !is.null(link)) nrow(link),
          # A new or removed mapping table changes which features meet.
          if (!is.null(fl)) paste(nrow(fl), paste(names(fl), collapse = ","),
                                  attr(fl, "source")$path %||% "",
                                  sum(nchar(unlist(fl, use.names = FALSE)), na.rm = TRUE)),
          sep = "|")
  })

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
            paste(info$covariates, collapse = ","), info$paired_col %||% "", info$primary_method,
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
                covariates    = if (length(info$covariates)) info$covariates else NULL,
                paired_col    = info$paired_col
              )
            } else {
              omicsCore::run_diff(
                input         = info$secondary,
                method        = sec_method,
                analysis_type = "group",
                group_col     = info$group_col,
                control_group = info$control,
                case_group    = info$all_cases,
                covariates    = if (length(info$covariates)) info$covariates else NULL,
                paired_col    = info$paired_col
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
          if (identical(m, "active_pathways")) args$organism <- ap_org
          res <- do.call(omicsCore::run_integration, args)
          list(result = res, sec_diff = sec_diff)
        },
        info = info, m = m, th = th, primary = primary, proj = proj,
        sec_diff = sec_diff, sec_method = sec_method,
        # The gene sets' species, read from how the gene names are
        # written, as the Enrichment view suggests it: mouse genes
        # against human sets lost every gene whose name differs
        # between the two (Trp53 / TP53).
        ap_org = if (!is.null(primary))
          guess_organism(primary$results$diff_result_df$feature_symbol) else "Hs"
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

  list(run_key = run_key, do_run = do_run)
}

# A metadata column of the partner layer matching `col`: the same name,
# the same name in another case, or -- for the group column -- the
# column the partner's recorded study design names.
match_partner_col <- function(col, partner) {
  if (is.null(col) || !nzchar(col)) return(NULL)
  nms <- names(partner$meta_df)
  if (col %in% nms) return(col)
  ci <- nms[tolower(nms) == tolower(col)]
  if (length(ci)) return(ci[[1L]])
  d <- tryCatch(omicsCore::study_design(partner), error = function(e) NULL)
  if (!is.null(d$group_col) && d$group_col %in% nms) return(d$group_col)
  NULL
}
