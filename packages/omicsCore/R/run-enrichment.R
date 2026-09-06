# Public entry point for ORA / GSEA. Dispatches to enrich-ora.R or
# enrich-gsea.R, then wraps results in an analysis_bundle. GSVA has its own
# entry point (run_gsva) since its inputs are an omics_input + gene-set list
# rather than a diff bundle.

SUPPORTED_ENRICH_TYPES <- c("ora", "gsea")

#' Run pathway enrichment from a differential bundle
#'
#' Single entry point for over-representation analysis (`type = "ora"`) and
#' rank-based gene-set enrichment analysis (`type = "gsea"`). Both backends
#' run against the requested MSigDB database in symbol space via
#' `clusterProfiler::enricher()` / `GSEA()`, so no `org.*` annotation
#' package is required.
#'
#' For ORA, features are split into up/down sets unless
#' `direction = "both"`; for GSEA, the full ranked vector is always passed
#' to the backend and `direction` filters the standardized output.
#'
#' For GSVA-style sample-level scoring, see [run_gsva()].
#'
#' @param diff_bundle An [`analysis_bundle`][is_analysis_bundle()] produced
#'   by [run_diff()].
#' @param type One of `"ora"` or `"gsea"`.
#' @param database One of `"hallmark"`, `"kegg"`, `"reactome"`,
#'   `"wikipathways"`, `"go_bp"`, `"go_mf"`, `"go_cc"`. Pass a character
#'   vector to query multiple databases.
#' @param organism Organism shorthand (e.g. `"Hs"`).
#' @param direction One of `"both"` (default), `"up"`, or `"down"`.
#' @param p_cutoff Significance cutoff for selecting diff features (ORA
#'   only; GSEA ranks the whole list).
#' @param output_p_cutoff Bound on the returned table. Defaults to
#'   `p_cutoff`. Pass `1` to keep every pathway, so a caller can choose
#'   raw or adjusted p at display time without re-running.
#' @param p_preference For ORA feature selection: `"adjusted"` (default) or
#'   `"raw"`.
#' @param effect_cutoff Optional |effect| cutoff for ORA feature selection.
#' @param p_adjust_method Multiple-testing correction method.
#' @param min_size,max_size GSEA min/max gene-set sizes.
#' @param ... Reserved for backend-specific extensions.
#'
#' @return An [`analysis_bundle`][is_analysis_bundle()] with
#'   `results$enrich_result_df` (standardized schema) and
#'   `results$enrich_object` (named list of clusterProfiler objects keyed by
#'   database, or for ORA by `<direction>__<database>`).
#' @export
#' @family enrich
#' @examples
#' \dontrun{
#'   diff <- run_diff(input, method = "limma", analysis_type = "group",
#'                    group_col = "treatment",
#'                    control_group = "DMSO", case_group = "Drug")
#'   enr <- run_enrichment(diff, type = "ora", database = "hallmark")
#'   head(enr$results$enrich_result_df)
#' }
run_enrichment <- function(
  diff_bundle,
  type = c("ora", "gsea"),
  database = c("hallmark", "kegg", "reactome", "go_bp", "go_mf", "go_cc",
               "wikipathways"),
  organism = "Hs",
  direction = c("both", "up", "down"),
  p_cutoff = 0.05,
  output_p_cutoff = NULL,
  p_preference = c("adjusted", "raw"),
  effect_cutoff = NULL,
  p_adjust_method = "BH",
  min_size = 10L,
  max_size = 500L,
  ...
) {
  if (!is_analysis_bundle(diff_bundle) ||
      !identical(diff_bundle$analysis_name, "run_diff")) {
    stop("`diff_bundle` must be an analysis_bundle from run_diff().")
  }

  type <- match.arg(type)
  direction <- match.arg(direction)
  p_preference <- match.arg(p_preference)

  # Accept a single value or vector and coerce through normalization.
  if (missing(database)) database <- "hallmark"
  assert_names(database, "database")
  assert_number(p_cutoff, "p_cutoff", lower = 0, upper = 1)
  assert_number(output_p_cutoff, "output_p_cutoff", lower = 0, upper = 1,
                allow_null = TRUE)
  assert_number(effect_cutoff, "effect_cutoff", lower = 0, allow_null = TRUE)
  assert_choice(p_adjust_method, "p_adjust_method", stats::p.adjust.methods)
  assert_count(min_size, "min_size", lower = 1L)
  assert_count(max_size, "max_size", lower = 1L)
  if (min_size > max_size) {
    stop("`min_size` must not exceed `max_size`.", call. = FALSE)
  }
  databases <- vapply(database, normalize_enrich_database, character(1L))
  databases <- unique(databases)
  organism <- normalize_organism(organism)

  per_db <- lapply(databases, function(db) {
    if (type == "ora") {
      run_ora_from_bundle(
        diff_bundle = diff_bundle,
        database = db,
        organism = organism,
        direction = direction,
        p_cutoff = p_cutoff,
        output_p_cutoff = output_p_cutoff,
        p_preference = p_preference,
        effect_cutoff = effect_cutoff,
        p_adjust_method = p_adjust_method
      )
    } else {
      run_gsea_from_bundle(
        diff_bundle = diff_bundle,
        database = db,
        organism = organism,
        direction = direction,
        p_cutoff = p_cutoff,
        output_p_cutoff = output_p_cutoff,
        p_adjust_method = p_adjust_method,
        min_size = min_size,
        max_size = max_size
      )
    }
  })
  names(per_db) <- databases

  enrich_object <- if (type == "ora") {
    flat <- list()
    for (db in databases) {
      objs <- per_db[[db]]$objects
      for (dir_label in names(objs)) {
        flat[[paste0(dir_label, "__", db)]] <- objs[[dir_label]]
      }
    }
    flat
  } else {
    lapply(per_db, `[[`, "object")
  }

  enrich_result_df <- dplyr::bind_rows(lapply(per_db, `[[`, "std"))
  check_enrich_result_schema(enrich_result_df)

  # Which pathway definitions produced this result. Matters for `kegg`,
  # where a refreshed cache holds current KEGG REST pathways while the
  # fallback is the 2011 MSigDB KEGG_LEGACY snapshot.
  geneset_sources <- stats::setNames(
    vapply(databases, geneset_table_source, character(1L), organism = organism),
    databases
  )

  new_analysis_bundle(
    analysis_name = "run_enrichment",
    input_info = diff_bundle$input_info,
    params = list(
      type = type,
      database = databases,
      organism = organism,
      direction = direction,
      p_cutoff = p_cutoff,
      p_preference = p_preference,
      effect_cutoff = effect_cutoff,
      p_adjust_method = p_adjust_method,
      min_size = min_size,
      max_size = max_size,
      geneset_sources = geneset_sources,
      comparison = diff_bundle$params$comparison
    ),
    results = list(
      enrich_result_df = enrich_result_df,
      enrich_object = enrich_object
    )
  )
}
