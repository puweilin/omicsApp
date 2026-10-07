# Render a project's analysis history as a runnable R script.
#
# Every analysis in omicsCore is a call to one of a handful of entry
# points with scalar or enumerated arguments, and `new_analysis_bundle()`
# already records those arguments in `params`. So the provenance needed
# to reconstruct a run is captured the moment it happens -- this file
# only renders it.
#
# The contract is that the script *reproduces*, not that it illustrates.
# A script that looks right but computes something else is worse than no
# script at all, because a reader will trust it. Two consequences:
#
#   * arguments are emitted as they were resolved, not as they were
#     requested. `run_diff(method = "auto")` becomes `method = "limma"`,
#     because that is the engine that ran.
#   * anything that cannot be reconstructed faithfully is marked with a
#     NOTE comment in the script rather than guessed at.

#' Arguments of `fn` that may appear in a generated call
#'
#' Bundle `params` also carries derived values (`comparison`,
#' `method_info`) that are outputs rather than arguments. Intersecting
#' with the formals drops them without needing a per-analysis exclusion
#' list.
#'
#' @param fn A function.
#'
#' @return Character vector of argument names, first argument and `...`
#'   removed, in signature order.
#' @keywords internal
#' @noRd
script_arg_names <- function(fn) {
  nms <- names(formals(fn))
  setdiff(nms[-1L], "...")
}

#' Render an R value as source code
#'
#' @param x Value to render.
#'
#' @return A single string, or `NA_character_` when `x` is not something
#'   that round-trips through `deparse()` as a literal.
#' @keywords internal
#' @noRd
render_value <- function(x) {
  if (is.null(x)) return("NULL")
  if (!is.atomic(x) || length(x) == 0L) return(NA_character_)
  out <- paste(deparse(x, width.cutoff = 500L), collapse = "")
  # Whether this is usable source is decided by trying it, not by
  # pattern-matching the text. An earlier version rejected anything
  # containing "<", which quietly dropped legitimate arguments -- a
  # group label like "<30" is an ordinary thing for a cohort study to
  # carry, and losing it produced an incomplete call rather than a
  # wrong one only by luck.
  #
  # `out` is the deparse of a value already in hand, so evaluating it
  # introduces nothing the caller did not already have; baseenv() keeps
  # it away from anything else.
  ok <- tryCatch(identical(eval(parse(text = out), envir = baseenv()), x),
                 error = function(e) FALSE)
  if (!isTRUE(ok)) return(NA_character_)
  out
}

#' Render one call with aligned, signature-ordered arguments
#'
#' @param fn_name Function name to emit.
#' @param first First (unnamed) argument, already rendered.
#' @param params Named list of recorded parameters.
#' @param arg_names Argument names to consider, in signature order.
#' @param assign_to Variable name to assign the result to, or `NULL`.
#'
#' @return A list with `lines` and `notes`.
#' @keywords internal
#' @noRd
render_call <- function(fn_name, first, params, arg_names, assign_to = NULL) {
  notes <- character(0)
  keep <- list()
  for (nm in arg_names) {
    if (!nm %in% names(params)) next
    value <- params[[nm]]
    if (is.null(value)) next
    rendered <- if (inherits(value, "script_code")) unclass(value) else render_value(value)
    if (is.na(rendered)) {
      notes <- c(notes, sprintf(
        "%s: `%s` was not a literal and is omitted; the call below is incomplete.",
        fn_name, nm))
      next
    }
    keep[[nm]] <- rendered
  }

  prefix <- if (is.null(assign_to)) "" else paste0(assign_to, " <- ")
  if (length(keep) == 0L) {
    return(list(lines = sprintf("%s%s(%s)", prefix, fn_name, first),
                notes = notes))
  }

  pad <- max(nchar(names(keep)))
  arg_lines <- vapply(names(keep), function(nm) {
    sprintf("  %-*s = %s", pad, nm, keep[[nm]])
  }, character(1))
  arg_lines <- paste0(c(paste0("  ", first), arg_lines), ",")
  arg_lines[length(arg_lines)] <- sub(",$", "", arg_lines[length(arg_lines)])

  list(
    lines = c(sprintf("%s%s(", prefix, fn_name), arg_lines, ")"),
    notes = notes
  )
}

# Code written into a call as it stands (a variable, a list of them).
script_code <- function(x) structure(x, class = "script_code")

# The run_diff() call that made a bundle, followed by the selection when
# the bundle is one contrast taken out of a shared fit: re-running only
# its two groups would give different p-values, since the shared fit
# pools the variance of every group in the model.
diff_call_lines <- function(params, input_var, var) {
  all_cases <- params$all_case_groups
  all_specs <- params$all_contrasts
  shown <- params$comparison
  selected <- FALSE
  if (length(shown) == 1L && length(all_specs) > 1L) {
    params$contrasts <- all_specs
    params$case_group <- NULL
    selected <- TRUE
  } else if (length(shown) == 1L && length(all_cases) > 1L) {
    params$case_group <- all_cases
    selected <- TRUE
  }
  call <- render_call("run_diff", input_var, params,
                      script_arg_names(run_diff), assign_to = var)
  lines <- call$lines
  if (selected) {
    lines <- c(lines, sprintf("%s <- select_comparison(%s, %s)", var, var,
                              render_value(shown)))
  }
  list(lines = lines, notes = call$notes, selected = selected)
}

# A data frame as the code that builds it, or NA when a column is not a
# literal.
render_data_frame <- function(df) {
  cols <- vapply(names(df), function(nm) {
    v <- render_value(as.vector(df[[nm]]))
    if (is.na(v)) NA_character_ else sprintf("%s = %s", render_name(nm), v)
  }, character(1))
  if (anyNA(cols)) return(NA_character_)
  sprintf("data.frame(%s, stringsAsFactors = FALSE)", paste(cols, collapse = ", "))
}

# The lines that rebuild the feature link an integration ran with, and
# the argument that passes it on. Nothing for a run that matched symbols
# (or one saved before links existed). A link read from a file is read
# from the archived copy; a small one that came from no file is written
# out; anything else gets a placeholder that fails where it stands,
# rather than a script that quietly matches symbols instead.
feature_link_script_lines <- function(src, experiments) {
  none <- list(lines = character(0), notes = character(0), arg = NULL)
  if (is.null(src) || identical(src$source, "symbol")) return(none)
  arg <- script_code("feature_link")
  if (!is.null(src$path)) {
    call <- render_call("read_feature_link",
                        render_value(file.path("raw", basename(src$path))),
                        params = list(columns = src$columns),
                        arg_names = "columns", assign_to = "feature_link")
    return(list(lines = c("# Which feature of one layer is which feature of the other:",
                          call$lines),
                notes = call$notes, arg = arg))
  }
  code <- if (is.data.frame(src$table)) render_data_frame(src$table) else NA_character_
  if (!is.na(code)) {
    return(list(lines = c("# Which feature of one layer is which feature of the other:",
                          paste0("feature_link <- ", code)),
                notes = character(0), arg = arg))
  }
  tags <- experiments %||% c("layer_a", "layer_b")
  cols <- stats::setNames(rep("<column>", length(tags)), tags)
  list(
    lines = c("# Which feature of one layer is which feature of the other:",
              sprintf("feature_link <- read_feature_link(%s, columns = %s)",
                      render_value("<path-to-feature-link-file>"), render_value(cols))),
    notes = sprintf(paste(
      "run_integration() matched features with a feature link (%s rows) that",
      "was not archived; fill in its file and columns below."),
      format(src$n_rows %||% NA_integer_)),
    arg = arg)
}

# The figures plot_integration() draws for each method.
INTEGRATION_FIGURES <- list(
  correlation = "scatter",
  concordance = c("dual_volcano", "effect_pair", "quadrant"),
  active_pathways = "dotplot"
)

# Text bound for a comment. A newline in a project name ended the
# comment and handed the rest of the name to the parser as code -- and
# an apostrophe in what followed opened a string that ran to the next
# one, so the script did not parse.
comment_text <- function(x) {
  gsub("[[:cntrl:]]+", " ", as.character(x))
}

# A tag as a list name in generated code: written bare when R would
# read it back as the same symbol, and in backticks otherwise.
render_name <- function(x) {
  if (identical(make.names(x), x)) return(x)
  paste0("`", gsub("`", "\\`", x, fixed = TRUE), "`")
}

section <- function(title) {
  c("", sprintf("# ---- %s %s", title,
                strrep("-", max(0L, 62L - nchar(title)))))
}

as_notes <- function(notes) {
  if (length(notes) == 0L) return(character(0))
  paste0("# NOTE: ", comment_text(notes))
}

# Find the experiment a bundle was computed on (bundle_layer()), and
# degrade to the first experiment, with a note elsewhere, when it is
# about none of them.
resolve_tag <- function(project, bundle) {
  tags <- names(project$experiments)
  if (length(tags) == 0L) return(NULL)
  hit <- bundle_layer(project, bundle)
  if (is.na(hit)) tags[[1L]] else hit
}

# One variable per layer. Through make.names(): a tag is whatever the
# user called the layer, and `input_rna seq` is not a symbol R can
# read back. Unique, so two tags that differ only in punctuation do
# not collapse onto one variable.
script_input_vars <- function(tags) {
  if (length(tags) == 1L) return(stats::setNames("input", tags))
  stats::setNames(make.names(paste0("input_", tags), unique = TRUE), tags)
}

#' Export a project's analysis history as a runnable R script
#'
#' Renders the calls that produced a project's analyses, in the order
#' they depend on each other, as an R script that reproduces them.
#' Arguments appear as they were *resolved*: an analysis run with
#' `method = "auto"` is emitted with the engine that auto selected, so
#' the script and the report agree.
#'
#' Whatever cannot be reconstructed faithfully is flagged with a `NOTE`
#' comment in the script rather than approximated — a script that reads
#' correctly but computes something else is worse than none.
#'
#' The input line points at the archived upload when the project records
#' one (see `source_path` on [omics_input()]); otherwise it emits a
#' placeholder path for the reader to fill in.
#'
#' @param project An [`omics_project`][is_omics_project()].
#' @param path Optional file to write to. The lines are returned either
#'   way.
#' @param include_plots Whether to append the plotting calls. On by
#'   default: the Shiny views draw through these same functions, so the
#'   figures the script produces are the figures the user was looking
#'   at. Set `FALSE` for a script that only recomputes the numbers.
#'
#' @return Character vector of script lines, invisibly when `path` is
#'   given.
#' @export
#' @family persistence
#' @examples
#' \dontrun{
#'   export_script(project, "reproduce.R")
#' }
export_script <- function(project, path = NULL, include_plots = TRUE) {
  if (!is_omics_project(project)) {
    stop("`project` must be an `omics_project`.")
  }
  assert_string(path, "path", allow_null = TRUE)
  assert_flag(include_plots, "include_plots")

  notes <- character(0)
  lines <- c(
    "# Reproducibility script generated by omicsCore::export_script().",
    "#",
    sprintf("# Project : %s", comment_text(project$name %||% "(unnamed)")),
    sprintf("# Written : %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    sprintf("# Versions: omicsCore %s | %s",
            as.character(utils::packageVersion("omicsCore")),
            R.version.string),
    "#",
    "# Arguments are shown as they were resolved at run time, so this",
    "# script performs the same computation the report describes.",
    "",
    "library(omicsCore)"
  )

  # ---- inputs ---------------------------------------------------------
  experiments <- project$experiments %||% list()
  n_exp <- length(experiments)
  input_vars <- character(0)
  if (n_exp > 0L) {
    lines <- c(lines, section("Input"))
    input_vars <- script_input_vars(names(experiments))
    for (tag in names(experiments)) {
      exp <- experiments[[tag]]
      var <- input_vars[[tag]]
      # A layer merged from Salmon / RSEM / kallisto files is read back
      # the way it was built, from the archived files under the names
      # they were uploaded with (those names are the sample names).
      if (!is.null(exp$quant_source$paths)) {
        qs <- exp$quant_source
        call <- render_call(
          "read_quant_files", render_value(file.path("raw", basename(qs$paths))),
          params = list(file_names = qs$names,
                        sample_sheet = if (!is.null(exp$sample_sheet_path))
                          file.path("raw", basename(exp$sample_sheet_path))),
          arg_names = script_arg_names(read_quant_files),
          assign_to = var)
        last <- length(call$lines)
        call$lines[last] <- paste0(call$lines[last], "$input")
        lines <- c(lines, "# read_quant_files() returns the merged input and its import report.",
                   call$lines)
        notes <- c(notes, call$notes)
        if (length(exp$excluded_samples)) {
          lines <- c(lines,
                     "# Samples excluded in QC:",
                     sprintf("%s <- subset_omics(%s, samples = setdiff(colnames(%s$expr_mat), %s))",
                             var, var, var, render_value(as.character(exp$excluded_samples))))
        }
        if (!is.null(exp$design$group_col)) {
          dc <- render_call("set_study_design", var,
                            params = exp$design[c("group_col", "reference")],
                            arg_names = script_arg_names(set_study_design),
                            assign_to = var)
          lines <- c(lines, dc$lines)
          notes <- c(notes, dc$notes)
        }
        next
      }
      src <- exp$source_path %||% NA_character_
      if (is.na(src)) {
        notes <- c(notes, sprintf(
          "the file behind '%s' was not archived; fill in the path below.",
          tag))
        src <- sprintf("<path-to-%s-file>", tag)
      } else {
        src <- file.path("raw", basename(src))
      }
      # A layer normalized in the app is read as the file holds it and
      # normalized again here. Reading it with the normalized label
      # skipped the transform and ran every step on linear values.
      norm <- exp$normalization
      file_assay <- if (!is.null(norm)) norm$from_assay_type %||% "raw_intensity"
                    else exp$assay_type
      if (is.null(norm) && !is.null(exp$normalized_mat)) {
        notes <- c(notes, sprintf(paste(
          "'%s' was normalized in the app before the normalization was",
          "recorded; add the normalize_omics() call it used after reading."), tag))
      }
      call <- render_call(
        "read_omics", render_value(src),
        params = list(omics_type = exp$omics_type,
                      assay_type = file_assay,
                      sheet_roles = exp$sheet_roles,
                      orientation = exp$orientation,
                      sample_sheet = if (!is.null(exp$sample_sheet_path))
                        file.path("raw", basename(exp$sample_sheet_path))),
        arg_names = script_arg_names(read_omics),
        assign_to = var
      )
      # `read_omics()` hands back the parsed input alongside its import
      # report, so the analysis functions want the `$input` element.
      # Emitting the bare call produced a script that read correctly and
      # then failed on the first `run_*()`.
      last <- length(call$lines)
      call$lines[last] <- paste0(call$lines[last], "$input")
      lines <- c(lines,
                 "# read_omics() returns the parsed input and its import report.",
                 call$lines)
      notes <- c(notes, call$notes)
      if (!is.null(norm)) {
        # `center` only when it was used, so a script for an uncentred
        # layer reads as it always has (and runs on older omicsCore).
        norm_args <- norm[c("method", "offset")]
        if (!identical(norm$center %||% "none", "none")) norm_args$center <- norm$center
        nc <- render_call("normalize_omics", var,
                          params = norm_args,
                          arg_names = script_arg_names(normalize_omics),
                          assign_to = var)
        lines <- c(lines, nc$lines)
        notes <- c(notes, nc$notes)
      }
      if (length(exp$excluded_samples)) {
        lines <- c(lines,
                   "# Samples excluded in QC:",
                   sprintf("%s <- subset_omics(%s, samples = setdiff(colnames(%s$expr_mat), %s))",
                           var, var, var, render_value(as.character(exp$excluded_samples))))
      }
      if (!is.null(exp$design$group_col)) {
        dc <- render_call("set_study_design", var,
                          params = exp$design[c("group_col", "reference")],
                          arg_names = script_arg_names(set_study_design),
                          assign_to = var)
        lines <- c(lines, dc$lines)
        notes <- c(notes, dc$notes)
      }
    }
  }

  bundles <- project$bundles %||% list()

  # A project can carry analyses without carrying the data they ran on:
  # `omics_project(experiments = list())` is legal, and a hand-built or
  # partially-restored project can reach here. Emitting the analysis
  # calls anyway would produce a script that references an `input` no
  # line defines and fails somewhere inside run_diff(). Name the gap and
  # fail at an obvious line instead.
  if (n_exp == 0L && length(bundles) > 0L) {
    notes <- c(notes, paste(
      "this project carries analyses but no imported data, so the script",
      "cannot run as written -- supply `input` yourself."))
    lines <- c(lines, section("Input"),
               "# The project held no experiment to read from.",
               "input <- NULL  # <- supply the omics_input these analyses used")
    input_vars[["__missing__"]] <- "input"
  }

  emit <- function(key, title, fn_name, fn, first, var) {
    bundle <- bundles[[key]]
    if (is.null(bundle)) return(NULL)
    call <- render_call(fn_name, first, bundle$params %||% list(),
                        script_arg_names(fn), assign_to = var)
    lines <<- c(lines, section(title), call$lines)
    notes <<- c(notes, call$notes)
  }

  input_for <- function(key) {
    tag <- resolve_tag(project, bundles[[key]])
    if (is.null(tag)) "input" else input_vars[[tag]] %||% "input"
  }

  if (!is.null(bundles$qc)) {
    emit("qc", "Quality control", "run_qc", run_qc, input_for("qc"), "qc")
  }
  diff_var_for_enrich <- "diff"
  if (!is.null(bundles$diff)) {
    # One contrast taken out of a shared fit is reproduced as that fit
    # followed by the selection. Re-running only its two groups would
    # give different p-values: the shared fit pools the variance of every
    # group in the model.
    shown <- bundles$diff$params$comparison
    dc <- diff_call_lines(bundles$diff$params, input_for("diff"), "diff")
    lines <- c(lines, section("Differential analysis"), dc$lines)
    notes <- c(notes, dc$notes)
    # A bundle holding every contrast of a run: the downstream steps each
    # ran on one of them, so that one is taken out for them by name.
    if (length(shown) > 1L &&
        (!is.null(bundles$enrich) || !is.null(bundles$integration))) {
      pick <- bundles$enrich$params$comparison %||%
        bundles$diff$params$shown_comparison %||% shown[[1L]]
      lines <- c(lines,
                 "# The steps below ran on one of these comparisons:",
                 sprintf("diff_shown <- select_comparison(diff, %s)",
                         render_value(pick[[1L]])))
      diff_var_for_enrich <- "diff_shown"
    }
  }
  if (!is.null(bundles$gsva)) {
    if (isTRUE(bundles$gsva$params$gene_sets_supplied)) {
      # The sets were the caller's own; they are in the bundle, not in
      # any database the script could name.
      bundles$gsva$params$database <- NULL
      notes <- c(notes, paste(
        "run_gsva() ran on gene sets supplied by hand; pass them as",
        "`gene_sets =` (they are in the bundle's results$gsva_gene_sets)."))
    }
    emit("gsva", "Gene-set variation", "run_gsva", run_gsva,
         input_for("gsva"), "gsva")
  }
  if (!is.null(bundles$enrich)) {
    emit("enrich", "Pathway enrichment", "run_enrichment", run_enrichment,
         diff_var_for_enrich, "enrich")
    if (is.null(bundles$diff)) {
      notes <- c(notes,
                 "enrichment ran on a differential result the project no longer holds.")
    }
  }
  if (!is.null(bundles$enrich_compare)) {
    # compare_enrichment() hands its settings to run_enrichment() through
    # `...`, so the arguments rendered are run_enrichment()'s.
    call <- render_call("compare_enrichment", "diff",
                        bundles$enrich_compare$params[
                          setdiff(names(bundles$enrich_compare$params), "comparison")],
                        script_arg_names(run_enrichment),
                        assign_to = "enrich_compare")
    lines <- c(lines, section("Enrichment across comparisons"), call$lines)
    notes <- c(notes, call$notes)
  }
  if (!is.null(bundles$integration)) {
    lines <- c(lines, section("Multi-omics integration"))
    # Built as its own block rather than by rewriting `lines` after the
    # fact: an earlier version ran a sub() over every line to append the
    # comma, which would have edited any other line that happened to
    # start with a `name` argument.
    ip <- bundles$integration$params %||% list()
    # The sample pairing correlation relies on, when the project has one.
    link <- project$sample_link
    link_code <- if (!is.null(link) && nrow(link) > 0L) render_data_frame(link)
    if (!is.null(link_code) && is.na(link_code)) {
      notes <- c(notes, "the project's sample_link could not be written out; set it before run_integration().")
      link_code <- NULL
    }
    lines <- c(
      lines,
      "project <- omics_project(",
      sprintf("  name        = %s,",
              render_value(project$name %||% "project")),
      sprintf("  experiments = list(%s)%s",
              paste(sprintf("%s = %s", vapply(names(input_vars), render_name, ""),
                            unname(input_vars)), collapse = ", "),
              if (!is.null(link_code)) "," else ""),
      if (!is.null(link_code)) sprintf("  sample_link = %s", link_code),
      ")"
    )
    # The feature link the run matched features with, read from the file
    # it was uploaded as (archived under raw/), so the pairs are the same.
    fl <- feature_link_script_lines(ip$feature_link_source, ip$experiments)
    lines <- c(lines, fl$lines)
    notes <- c(notes, fl$notes)
    ip$feature_link <- fl$arg
    # Each layer's differential result, made the way it was made --
    # including the partner layer's, which the app computes itself.
    if (length(ip$diff_params)) {
      dvars <- character(0)
      for (tag in names(ip$diff_params)) {
        var <- make.names(paste0("diff_", tag))
        dc <- diff_call_lines(ip$diff_params[[tag]],
                              input_vars[[tag]] %||% "input", var)
        lines <- c(lines, sprintf("# The differential result for '%s':", comment_text(tag)),
                   dc$lines)
        notes <- c(notes, dc$notes)
        dvars[[tag]] <- var
      }
      ip$diff_bundles <- script_code(sprintf("list(%s)", paste(
        sprintf("%s = %s", vapply(names(dvars), render_name, ""), dvars),
        collapse = ", ")))
    } else if (!identical(ip$method, "correlation")) {
      notes <- c(notes, paste(
        "run_integration() was given differential results that this project",
        "does not record; supply `diff_bundles =` for both layers."))
    }
    # The method's own settings travel through `...`, so they are named
    # here; omitting them ran every method at its defaults.
    method_args <- setdiff(names(ip), c(script_arg_names(run_integration), "method_info",
                                        "diff_params", "feature_link_source"))
    call <- render_call("run_integration", "project", ip,
                        c(script_arg_names(run_integration), method_args),
                        assign_to = "integration")
    lines <- c(lines, call$lines)
    notes <- c(notes, call$notes)
  }

  if (isTRUE(include_plots)) {
    lines <- c(lines, section("Figures"))
    if (!is.null(bundles$qc)) {
      lines <- c(lines,
                 'plot_qc(qc, view = "pca")',
                 'plot_qc(qc, view = "missing")')
    }
    if (!is.null(bundles$diff)) {
      lines <- c(lines,
        "# Drawn at plot_volcano()'s default cut, which is the figure the",
        "# app shows: its threshold sliders filter the hit table, not this.",
        "# The cut is printed in the plot's caption.",
        # One comparison per figure: a bundle holding several is refused.
        if (length(bundles$diff$params$comparison) > 1L) {
          if (identical(diff_var_for_enrich, "diff_shown")) "plot_volcano(diff_shown)"
          else sprintf("plot_volcano(select_comparison(diff, %s))", render_value(
            (bundles$diff$params$shown_comparison %||% bundles$diff$params$comparison)[[1L]]))
        } else "plot_volcano(diff)")
    }
    if (!is.null(bundles$enrich)) {
      lines <- c(lines, 'plot_enrichment(enrich, view = "dot", top_n = 12L)')
    }
    if (!is.null(bundles$integration)) {
      views <- INTEGRATION_FIGURES[[bundles$integration$params$method %||% "concordance"]]
      lines <- c(lines, sprintf('plot_integration(integration, view = "%s")', views))
    }
  }

  lines <- c(lines, section("Session"), "sessionInfo()")

  if (length(notes) > 0L) {
    # Notes go at the top: a reader must meet the caveats before the
    # code, not after they have already run it.
    header_end <- which(lines == "library(omicsCore)")[[1L]]
    lines <- append(lines, c("", as_notes(notes)), after = header_end - 2L)
  }

  if (!is.null(path)) {
    writeLines(lines, path)
    return(invisible(lines))
  }
  lines
}
