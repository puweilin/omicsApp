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

  rmd_path <- template %||% system.file(
    "rmd", "default-report.Rmd", package = "omicsCore"
  )
  if (!nzchar(rmd_path) || !file.exists(rmd_path)) {
    stop("Report template not found: ", rmd_path)
  }

  output_format <- switch(format,
    html = rmarkdown::html_document(self_contained = TRUE),
    pdf  = rmarkdown::pdf_document()
  )

  out_dir <- dirname(path)
  if (!nzchar(out_dir)) out_dir <- "."
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # Knitting writes its figure directory beside the input by default,
  # and the input is the template inside the installed package -- a
  # location that is read-only on a shared server and shared by every
  # session on any other. A private scratch directory instead, gone
  # when the render is.
  scratch <- tempfile("omicsCore-report-")
  dir.create(scratch)
  on.exit(unlink(scratch, recursive = TRUE), add = TRUE)

  rmarkdown::render(
    input = rmd_path,
    output_file = basename(path),
    output_dir = out_dir,
    intermediates_dir = scratch,
    output_format = output_format,
    params = list(
      project = project,
      project_name = project$name %||% "omics_project"
    ),
    envir = new.env(parent = globalenv()),
    quiet = TRUE
  )

  invisible(path)
}
