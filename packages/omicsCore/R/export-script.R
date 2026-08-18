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
    rendered <- render_value(value)
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

section <- function(title) {
  c("", sprintf("# ---- %s %s", title,
                strrep("-", max(0L, 62L - nchar(title)))))
}

as_notes <- function(notes) {
  if (length(notes) == 0L) return(character(0))
  paste0("# NOTE: ", notes)
}

# Find the experiment a bundle was computed on. Bundles record the
# omics_type they ran against, which is also the tag app_server() files
# experiments under, so this resolves in the common case and degrades to
# the first experiment with a note otherwise.
resolve_tag <- function(project, bundle) {
  want <- bundle$input_info$omics_type
  tags <- names(project$experiments)
  if (length(tags) == 0L) return(NULL)
  if (is.null(want)) return(tags[[1L]])
  types <- vapply(project$experiments,
                  function(e) e$omics_type %||% "", character(1))
  hit <- tags[types %in% want]
  if (length(hit) == 0L) tags[[1L]] else hit[[1L]]
}

script_input_var <- function(tag, n_experiments) {
  if (n_experiments == 1L) "input" else paste0("input_", tag)
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

  notes <- character(0)
  lines <- c(
    "# Reproducibility script generated by omicsCore::export_script().",
    "#",
    sprintf("# Project : %s", project$name %||% "(unnamed)"),
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
    for (tag in names(experiments)) {
      exp <- experiments[[tag]]
      var <- script_input_var(tag, n_exp)
      input_vars[[tag]] <- var
      src <- exp$source_path %||% NA_character_
      if (is.na(src)) {
        notes <- c(notes, sprintf(
          "the file behind '%s' was not archived; fill in the path below.",
          tag))
        src <- sprintf("<path-to-%s-file>", tag)
      } else {
        src <- file.path("raw", basename(src))
      }
      call <- render_call(
        "read_omics", render_value(src),
        params = list(omics_type = exp$omics_type,
                      assay_type = exp$assay_type),
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
  if (!is.null(bundles$diff)) {
    emit("diff", "Differential analysis", "run_diff", run_diff,
         input_for("diff"), "diff")
  }
  if (!is.null(bundles$gsva)) {
    emit("gsva", "Gene-set variation", "run_gsva", run_gsva,
         input_for("gsva"), "gsva")
  }
  if (!is.null(bundles$enrich)) {
    emit("enrich", "Pathway enrichment", "run_enrichment", run_enrichment,
         "diff", "enrich")
    if (is.null(bundles$diff)) {
      notes <- c(notes,
                 "enrichment ran on a differential result the project no longer holds.")
    }
  }
  if (!is.null(bundles$integration)) {
    lines <- c(lines, section("Multi-omics integration"))
    # Built as its own block rather than by rewriting `lines` after the
    # fact: an earlier version ran a sub() over every line to append the
    # comma, which would have edited any other line that happened to
    # start with a `name` argument.
    lines <- c(
      lines,
      "project <- omics_project(",
      sprintf("  name        = %s,",
              render_value(project$name %||% "project")),
      sprintf("  experiments = list(%s)",
              paste(sprintf("%s = %s", names(input_vars),
                            unname(input_vars)), collapse = ", ")),
      ")"
    )
    call <- render_call("run_integration", "project",
                        bundles$integration$params %||% list(),
                        script_arg_names(run_integration),
                        assign_to = "integration")
    lines <- c(lines, call$lines)
    notes <- c(notes, call$notes)
    if (!"diff_bundles" %in% names(bundles$integration$params %||% list())) {
      notes <- c(notes, paste(
        "run_integration() was given differential bundles that are not",
        "recorded in params; add `diff_bundles = list(diff)` if the",
        "method needs them."))
    }
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
        "# The app's threshold sliders re-colour the volcano for reading;",
        "# they do not change the analysis. This draws significance as",
        "# `run_diff()` determined it, which is the version worth citing.",
        "plot_volcano(diff)")
    }
    if (!is.null(bundles$enrich)) {
      lines <- c(lines, 'plot_enrichment(enrich, view = "dot", top_n = 12L)')
    }
    if (!is.null(bundles$integration)) {
      lines <- c(lines,
                 'plot_integration(integration, view = "dual_volcano")',
                 'plot_integration(integration, view = "effect_pair")')
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
