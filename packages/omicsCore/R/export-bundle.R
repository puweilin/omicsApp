# Write the data side of an analysis_bundle (tables, matrices, plots,
# provenance) to disk. The contract from `docs/export-manifest.md` is:
# every disk-writing function name starts with `export_` or `save_`, and
# the return value is an artifact registry describing what landed where.

#' Export an analysis_bundle to disk
#'
#' Writes the data-bearing fields of a bundle (`results$*_df` /
#' `results$*_matrix`), provenance (`params`), and -- if `plots` is
#' non-empty -- the supplied ggplot/heatmap objects, in a target
#' directory.
#'
#' Output layout under `dir`:
#'
#' ```
#' <dir>/
#'   <prefix><bundle_name>_<table>.xlsx
#'   <prefix><bundle_name>_<table>.tsv
#'   <prefix><bundle_name>_<matrix>.tsv
#'   <prefix><bundle_name>_<plot>.pdf      (or .png)
#'   <prefix><bundle_name>_params.json
#' ```
#'
#' `prefix` defaults to `""` (no prefix). The bundle's existing artifact
#' registry is merged with the rows added by this call so downstream
#' code can keep tracking every file from one entry point.
#'
#' @param bundle An [`analysis_bundle`][is_analysis_bundle()].
#' @param dir Output directory. Created if it does not exist.
#' @param formats Character vector of table/plot formats. Tables support
#'   `"xlsx"` and `"tsv"`; plots support `"pdf"` and `"png"`. Defaults to
#'   `c("xlsx", "tsv", "pdf")`.
#' @param prefix Optional file-name prefix. Useful when exporting
#'   multiple bundles to the same directory.
#' @param plots Optional named list of plot objects to write alongside the
#'   tables (e.g. `list(volcano = plot_volcano(bundle))`). Both `ggplot`
#'   and `ComplexHeatmap` objects are supported.
#' @param width,height Plot dimensions in inches. Defaults `7 x 5`.
#'
#' @return A `data.frame` artifact registry (one row per written file).
#' @export
#' @family persistence
#' @examples
#' \dontrun{
#'   d <- run_diff(input, method = "ttest", analysis_type = "group",
#'                 group_col = "g", control_group = "ctrl", case_group = "case")
#'   export_bundle(d, dir = "out/", plots = list(volcano = plot_volcano(d)))
#' }
export_bundle <- function(
  bundle,
  dir,
  formats = c("xlsx", "tsv", "pdf"),
  prefix = NULL,
  plots = NULL,
  width = 7,
  height = 5
) {
  if (!is_analysis_bundle(bundle)) {
    stop("`bundle` must be an analysis_bundle.")
  }
  if (!is.character(dir) || length(dir) != 1L || !nzchar(dir)) {
    stop("`dir` must be a non-empty single string.")
  }
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)

  formats <- match.arg(formats,
                       choices = c("xlsx", "tsv", "pdf", "png"),
                       several.ok = TRUE)
  prefix <- prefix %||% ""
  base <- file.path(dir, paste0(prefix, bundle$analysis_name))

  registry <- bundle$artifacts %||% new_artifact_registry()

  # ---- tables / matrices --------------------------------------------
  for (nm in names(bundle$results)) {
    obj <- bundle$results[[nm]]
    if (is.data.frame(obj)) {
      registry <- write_table(obj, paste0(base, "_", nm), formats, registry, label = nm)
    } else if (is.matrix(obj) && is.numeric(obj)) {
      registry <- write_matrix(obj, paste0(base, "_", nm), registry, label = nm)
    }
    # Lists / model objects are intentionally skipped here; they round-trip
    # via save_project() instead.
  }

  # ---- provenance ----------------------------------------------------
  params_path <- paste0(base, "_params.json")
  json_payload <- list(
    analysis_name = bundle$analysis_name,
    input_info = bundle$input_info,
    params = bundle$params,
    messages = bundle$messages,
    warnings = bundle$warnings,
    exported_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  )
  writeLines(
    jsonlite::toJSON(json_payload, auto_unbox = TRUE, pretty = TRUE,
                     force = TRUE, null = "null", na = "null"),
    con = params_path
  )
  registry <- register_artifact(registry, "params", "bundle_params", params_path)

  # ---- plots ---------------------------------------------------------
  if (!is.null(plots)) {
    if (!is.list(plots) || is.null(names(plots))) {
      stop("`plots` must be a named list of plot objects.")
    }
    plot_fmts <- intersect(c("pdf", "png"), formats)
    if (length(plot_fmts) == 0L) {
      plot_fmts <- "pdf"  # default to PDF if caller asked for tables only
    }
    for (nm in names(plots)) {
      for (fmt in plot_fmts) {
        path <- paste0(base, "_", nm, ".", fmt)
        write_plot(plots[[nm]], path, fmt, width = width, height = height)
        registry <- register_artifact(registry, "plot", nm, path)
      }
    }
  }

  rownames(registry) <- NULL
  registry
}

# ---- internal helpers --------------------------------------------------

write_table <- function(df, base, formats, registry, label) {
  if ("xlsx" %in% formats) {
    path <- paste0(base, ".xlsx")
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "data")
    openxlsx::writeData(wb, "data", df)
    openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
    registry <- register_artifact(registry, "table", label, path)
  }
  if ("tsv" %in% formats) {
    path <- paste0(base, ".tsv")
    utils::write.table(df, file = path, sep = "\t",
                       quote = FALSE, row.names = FALSE)
    registry <- register_artifact(registry, "table", label, path)
  }
  registry
}

write_matrix <- function(mat, base, registry, label) {
  path <- paste0(base, ".tsv")
  # Matrices stay TSV-only -- xlsx column limits hurt wide pathway/sample tables.
  df <- data.frame(row_id = rownames(mat), mat, check.names = FALSE,
                   stringsAsFactors = FALSE)
  utils::write.table(df, file = path, sep = "\t",
                     quote = FALSE, row.names = FALSE)
  register_artifact(registry, "matrix", label, path)
}

write_plot <- function(plot, path, fmt, width, height) {
  if (inherits(plot, "ggplot")) {
    ggplot2::ggsave(filename = path, plot = plot,
                    width = width, height = height, units = "in",
                    device = fmt)
    return(invisible(path))
  }
  if (inherits(plot, c("Heatmap", "HeatmapList"))) {
    if (!is_installed("ComplexHeatmap")) {
      stop("ComplexHeatmap is required to export heatmap plots.")
    }
    open_device(path, fmt, width, height)
    on.exit(grDevices::dev.off(), add = TRUE)
    ComplexHeatmap::draw(plot)
    return(invisible(path))
  }
  if (inherits(plot, "recordedplot")) {
    open_device(path, fmt, width, height)
    on.exit(grDevices::dev.off(), add = TRUE)
    grDevices::replayPlot(plot)
    return(invisible(path))
  }
  stop("Unsupported plot type for export: ", paste(class(plot), collapse = "/"))
}

open_device <- function(path, fmt, width, height) {
  switch(fmt,
    pdf = grDevices::pdf(file = path, width = width, height = height),
    png = grDevices::png(filename = path, width = width, height = height,
                         units = "in", res = 150),
    stop("Unsupported plot format: ", fmt)
  )
}
