# Render a multi-panel HTML (or PDF) report for an omics_project. The
# heavy lifting lives in `inst/rmd/default-report.Rmd`; this file just
# resolves the template path and calls `rmarkdown::render()` with the
# right params. `rmarkdown` is Suggests-gated.

ensure_rmarkdown <- function() {
  if (!is_installed("rmarkdown")) {
    stop(
      "Package 'rmarkdown' is required for export_report(). ",
      "Install with: install.packages('rmarkdown').",
      call. = FALSE
    )
  }
}

#' Render an HTML report for an omics_project
#'
#' Renders a self-contained HTML (or PDF) report summarising every
#' experiment and any `analysis_bundle` objects attached to a project.
#' The bundled `inst/rmd/default-report.Rmd` template is used unless a
#' custom path is supplied.
#'
#' Bundles can be attached to the project in a `bundles` slot
#' (e.g. `project$bundles$diff <- run_diff(...)`); the template walks
#' over `names(project$bundles)` and emits per-bundle parameter and
#' top-rows tables.
#'
#' Requires the `rmarkdown` package.
#'
#' @param project An [`omics_project`][is_omics_project()].
#' @param path Output path (file extension is inferred from `format`).
#' @param format `"html"` (default) or `"pdf"`. PDF requires a working
#'   LaTeX engine in the user environment.
#' @param template Optional path to a custom Rmd template. Defaults to
#'   the package's `inst/rmd/default-report.Rmd`.
#' @param overwrite If `FALSE` (default) and `path` exists, raise an
#'   error.
#'
#' @return Invisibly returns `path`.
#' @export
#' @family persistence
#' @examples
#' \dontrun{
#'   p <- omics_project("demo", experiments = list(proteo = my_input))
#'   p$bundles <- list(diff = run_diff(my_input, ...))
#'   export_report(p, "demo.html")
#' }
export_report <- function(
  project,
  path,
  format = c("html", "pdf"),
  template = NULL,
  overwrite = FALSE
) {
  assert_string(template, "template", allow_null = TRUE)
  assert_flag(overwrite, "overwrite")
  if (!is_omics_project(project)) {
    stop("`project` must be an `omics_project`.")
  }
  format <- match.arg(format)
  if (!is.character(path) || length(path) != 1L || !nzchar(path)) {
    stop("`path` must be a non-empty single string.")
  }
  if (file.exists(path) && !isTRUE(overwrite)) {
    stop("Path already exists: ", path,
         " (pass `overwrite = TRUE` to replace).")
  }
  ensure_rmarkdown()
  if (identical(format, "pdf")) ensure_pdf_toolchain()

  rmd_path <- template %||% system.file(
    "rmd", "default-report.Rmd", package = "omicsCore"
  )
  if (!nzchar(rmd_path) || !file.exists(rmd_path)) {
    stop("Report template not found: ", rmd_path)
  }

  out_dir <- dirname(path)
  if (!nzchar(out_dir)) out_dir <- "."
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # Knitting writes its figure directory beside the input by default,
  # and the input is the template inside the installed package -- a
  # location that is read-only on a shared server and shared by every
  # session on any other. A private scratch directory instead, gone
  # when the render is.
  #
  # Created before the output format, because the PDF format writes a
  # LaTeX preamble and that belongs in here too -- left in tempdir() it
  # accumulates one file per report for the life of the session.
  scratch <- tempfile("omicsCore-report-")
  dir.create(scratch)
  on.exit(unlink(scratch, recursive = TRUE), add = TRUE)

  # rmarkdown's self-contained path copies each figure through a bare
  # `tempfile()` and does not remove it, so a report with plots leaves one
  # PNG per figure in tempdir() for the life of the session. On a server
  # where a container runs all day and every user renders reports, that is
  # the container's disk.
  #
  # It is rmarkdown's file, but this call is what caused it, and the window
  # is this function -- the render is synchronous, so anything that appears
  # at the top of tempdir() during it came from here.
  #
  # Never the report itself, though. `path` is frequently a bare
  # `tempfile()` -- that is exactly what Shiny's downloadHandler passes its
  # `content` function -- which puts the report in tempdir() under a name
  # matching the very pattern below. This cleanup deleted the file it had
  # just been asked to write, and did it on the way out, so the render
  # succeeded, the return value was the path, and the path was empty.
  before_tmp <- list.files(tempdir(), all.files = TRUE, no.. = TRUE)
  protected <- basename(path)
  on.exit({
    left <- setdiff(list.files(tempdir(), all.files = TRUE, no.. = TRUE),
                    before_tmp)
    left <- setdiff(left[grepl("^file[0-9a-f]+$", left)], protected)
    for (f in left) {
      p_ <- file.path(tempdir(), f)
      if (!dir.exists(p_)) unlink(p_)
    }
  }, add = TRUE)

  output_format <- switch(format,
    html = rmarkdown::html_document(self_contained = TRUE),
    pdf  = pdf_output_format(scratch)
  )

  # Render into the scratch directory under a name whose extension matches
  # the format, then move the result to `path`. Two reasons, and the first
  # one is why the download buttons produced nothing.
  #
  # rmarkdown appends the format's extension to `output_file` when it does
  # not already carry one, so asking it for `fileab12cd` yielded
  # `fileab12cd.html` -- and returned without an error, because from its
  # point of view nothing went wrong. Every test here passes
  # `tempfile(fileext = ".html")` and so never saw it; Shiny's
  # downloadHandler passes an extension-less temp path and saw nothing else.
  #
  # The second reason is XeLaTeX, and it is why the PDF button failed on
  # the server while working here. `intermediates_dir` moves knitr's files
  # but not LaTeX's: the engine still wrote `report.aux` and `report.log`
  # beside the input, and the input is the template inside the installed
  # package. On a development checkout that directory is writable, so the
  # only trace was stray .log files appearing in inst/rmd -- two of which
  # were committed. On the server it is /usr/local/lib/R/site-library,
  # owned by root and read-only to the account the container runs as, and
  # LaTeX stopped with `I can't write on file`.
  #
  # So the template is rendered as a copy inside the scratch directory. It
  # references nothing outside itself, so the copy is complete.
  render_name <- paste0("report", switch(format, html = ".html", pdf = ".pdf"))
  local_rmd <- file.path(scratch, "report.Rmd")
  if (!file.copy(rmd_path, local_rmd)) {
    stop("Could not stage the report template for rendering.", call. = FALSE)
  }

  # tinytex removes its working files on failure, the log among them. That
  # log holds the actual LaTeX error; without it the message is "LaTeX
  # failed to compile <file>" plus a link to a web page, addressed to
  # someone who has neither the file nor a shell on the server.
  old_opts <- options(tinytex.clean = FALSE)
  on.exit(options(old_opts), add = TRUE)

  tryCatch(
    rmarkdown::render(
      input = local_rmd,
      output_file = render_name,
      output_dir = scratch,
      intermediates_dir = scratch,
      # knitr writes figures to `fig.path`, which is relative -- and
      # relative to the knit root, not to intermediates_dir. Unset, the
      # report's plots landed in tempdir() and stayed there: one PNG per
      # figure per report, for the life of the session. Rooting the knit in
      # the scratch directory means they go where the cleanup already is.
      knit_root_dir = scratch,
      output_format = output_format,
      params = list(
        project = project,
        project_name = project$name %||% "omics_project"
      ),
      envir = new.env(parent = globalenv()),
      quiet = TRUE
    ),
    error = function(e) {
      stop(render_failure_message(e, scratch), call. = FALSE)
    }
  )

  produced <- file.path(scratch, render_name)
  if (!file.exists(produced)) {
    stop("The report render finished but produced no file.", call. = FALSE)
  }
  # file.copy() rather than file.rename(): scratch is under tempdir() and
  # `path` is wherever the caller asked, which on a server is a different
  # filesystem -- rename() cannot cross one, and it says so by returning
  # FALSE rather than raising. The scratch copy goes with the directory
  # on exit.
  if (!file.copy(produced, path, overwrite = TRUE)) {
    stop("Could not write the report to ", path, call. = FALSE)
  }

  invisible(path)
}

# The last error LaTeX reported, for a message a user will actually read.
#
# A .log line starting with `!` is TeX's error; the `l.<n>` line after it is
# the source line that caused it. Both are worth carrying. Returns the
# original message unchanged when there is no log to read -- an HTML render
# has none, and neither does a failure that happened before LaTeX started.
render_failure_message <- function(e, scratch) {
  msg <- conditionMessage(e)
  logs <- list.files(scratch, pattern = "\\.log$", full.names = TRUE,
                     recursive = TRUE)
  if (length(logs) == 0L) return(msg)
  logs <- logs[order(file.info(logs)$mtime, decreasing = TRUE)]
  txt <- tryCatch(readLines(logs[[1L]], warn = FALSE),
                  error = function(e) character(0))
  hits <- grep("^!", txt)
  if (length(hits) == 0L) return(msg)
  detail <- txt[seq(hits[[1L]], min(length(txt), hits[[1L]] + 3L))]
  paste0(msg, "\nLaTeX said: ", paste(trimws(detail), collapse = " "))
}


# ---- PDF ---------------------------------------------------------------

# The font to typeset CJK with, first one present wins. Named rather than
# discovered because a fallback the user cannot see is how a report of
# empty boxes gets produced: XeLaTeX does not stop when a glyph is missing.
CJK_FONT_CANDIDATES <- c("Noto Sans CJK SC", "Noto Serif CJK SC",
                         "Source Han Sans SC", "Songti SC",
                         "PingFang SC", "SimSun")

#' Which CJK font this machine can typeset with
#'
#' @return A font name, or `NA_character_` when none is installed.
#' @keywords internal
#' @noRd
available_cjk_font <- function(candidates = CJK_FONT_CANDIDATES) {
  listed <- tryCatch(
    suppressWarnings(system2("fc-list", ":family", stdout = TRUE, stderr = FALSE)),
    error = function(e) character(0))
  if (length(listed) == 0L) return(NA_character_)
  hit <- candidates[vapply(candidates, function(f) {
    any(grepl(f, listed, fixed = TRUE))
  }, logical(1))]
  if (length(hit) == 0L) NA_character_ else hit[[1L]]
}

#' Is there a working PDF toolchain here?
#'
#' rmarkdown being installed says nothing about it: the engine and the
#' fonts are system packages, and the Report view used to offer the button
#' on the strength of the R package alone.
#'
#' @return `list(ok = , message = , engine = , cjk_font = )`.
#' @keywords internal
#' @export
pdf_toolchain_status <- function() {
  if (!is_installed("rmarkdown")) {
    return(list(ok = FALSE, engine = NA_character_, cjk_font = NA_character_,
                message = "rmarkdown is not installed."))
  }
  if (!rmarkdown::pandoc_available()) {
    return(list(ok = FALSE, engine = NA_character_, cjk_font = NA_character_,
                message = "pandoc is not available."))
  }
  engine <- if (nzchar(Sys.which("xelatex"))) "xelatex" else NA_character_
  if (is.na(engine)) {
    return(list(ok = FALSE, engine = NA_character_, cjk_font = NA_character_,
                message = paste0(
                  "No XeLaTeX engine found. PDF export needs one: pdfLaTeX ",
                  "stops on any non-Latin character, which includes most ",
                  "project names this application is used for.")))
  }
  font <- available_cjk_font()
  if (is.na(font)) {
    # Deliberately a refusal rather than a warning. XeLaTeX without a CJK
    # font does not fail -- it typesets empty boxes -- and a report of
    # blank rectangles is worse than no report.
    return(list(ok = FALSE, engine = engine, cjk_font = NA_character_,
                message = paste0(
                  "XeLaTeX is present but no CJK font is installed, so any ",
                  "non-Latin text would be typeset as empty boxes. Install ",
                  "fonts-noto-cjk.")))
  }
  list(ok = TRUE, engine = engine, cjk_font = font, message = "")
}

ensure_pdf_toolchain <- function() {
  st <- pdf_toolchain_status()
  if (isTRUE(st$ok)) return(invisible(TRUE))
  stop("Cannot render a PDF here: ", st$message, call. = FALSE)
}

pdf_output_format <- function(scratch = tempdir()) {
  st <- pdf_toolchain_status()
  rmarkdown::pdf_document(
    latex_engine = st$engine %||% "xelatex",
    includes = rmarkdown::includes(
      in_header = cjk_header_file(st$cjk_font, scratch))
  )
}

# xeCJK selects a font per script, so Latin text keeps its own face and
# only the CJK runs switch. Setting the main font to a CJK face instead
# would typeset the whole report in it.
cjk_header_file <- function(font, scratch = tempdir()) {
  path <- tempfile("omicsCore-cjk-", tmpdir = scratch, fileext = ".tex")
  writeLines(c(
    "\\usepackage{xeCJK}",
    sprintf("\\setCJKmainfont{%s}", font),
    sprintf("\\setCJKsansfont{%s}", font),
    sprintf("\\setCJKmonofont{%s}", font),
    "\\XeTeXlinebreaklocale \"zh\"",
    "\\XeTeXlinebreakskip = 0pt plus 1pt"
  ), path)
  path
}
