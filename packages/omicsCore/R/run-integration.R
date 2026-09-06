# Public entry point for cross-omics integration. Dispatches to one of the
# three backends in `R/integration-*.R` and wraps the result in an
# analysis_bundle. Backends differ in what they need:
#
#   * `correlation`     -- needs the expression matrices on disk in the
#                          project; no diff bundles required.
#   * `concordance`     -- needs a `diff_bundles` list keyed by the two
#                          experiment tags.
#   * `active_pathways` -- needs `diff_bundles` plus a MSigDB database and
#                          the `ActivePathways` Suggests package.

SUPPORTED_INTEGRATION_METHODS <- c("correlation", "concordance", "active_pathways")

#' Run a cross-omics integration analysis
#'
#' Single entry point for dual-omics integration. Three methods are
#' supported:
#'
#' * `"correlation"` -- per-feature Pearson or Spearman correlation across
#'   paired samples from two experiments. Useful for RNA / protein layers
#'   that share donors.
#' * `"concordance"` -- per-feature agreement between two differential
#'   analyses, classified into the four sign quadrants and combined via
#'   Fisher's method.
#' * `"active_pathways"` -- pathway-level combined-p enrichment via the
#'   `ActivePathways` package, fed by the p-values of two differential
#'   analyses.
#'
#' All methods return an `analysis_bundle` whose `results$integration_df`
#' follows the schema documented in [check_integration_result_schema()].
#'
#' @param project An [`omics_project`][is_omics_project()] containing the
#'   two experiments to integrate.
#' @param method One of `"correlation"`, `"concordance"`, `"active_pathways"`.
#' @param experiments Length-2 character vector of experiment tags. If
#'   `NULL` and the project has exactly two layers, both are used.
#' @param diff_bundles Named list of [`run_diff()`] bundles keyed by
#'   experiment tag. Required for `"concordance"` and `"active_pathways"`.
#' @param by Feature key used to join layers (default `"feature_symbol"`).
#' @param ... Method-specific arguments forwarded to the backend. See
#'   the per-method sections below.
#'
#' @section Correlation arguments:
#' * `cor_method` -- `"spearman"` (default) or `"pearson"`.
#' * `p_adjust_method` -- defaults to `"BH"`.
#' * `min_samples` -- minimum paired samples (default `4`).
#' * `p_cutoff` -- significance cutoff (default `0.05`).
#'
#' @section Concordance arguments:
#' * `p_preference` -- `"adjusted"` (default) or `"raw"`.
#' * `p_cutoff` -- significance cutoff (default `0.05`).
#' * `p_adjust_method` -- defaults to `"BH"`.
#'
#' @section ActivePathways arguments:
#' * `database` -- MSigDB shorthand (default `"hallmark"`).
#' * `organism` -- defaults to `"Hs"`.
#' * `p_preference` -- `"raw"` (default) or `"adjusted"`. ActivePathways
#'   prefers raw p-values since it applies its own multiple-testing
#'   correction across pathways.
#' * `significant` -- pathway-level cutoff (default `0.05`).
#' * `geneset_filter` -- length-2 integer vector of min/max gene-set sizes
#'   (default `c(5L, 1000L)`).
#' * `merge_method` -- `"Brown"` (default) or another method accepted by
#'   `ActivePathways::ActivePathways()`.
#'
#' @return An [`analysis_bundle`][is_analysis_bundle()] with
#'   `results$integration_df` (standardized schema) and, for
#'   `"active_pathways"`, `results$integration_raw` carrying the raw
#'   `ActivePathways` table.
#' @export
#' @family integration
#' @examples
#' \dontrun{
#'   # Correlation across paired RNA / protein samples
#'   cor_bundle <- run_integration(project, method = "correlation",
#'                                 experiments = c("rna", "prot"))
#'
#'   # Concordance from two diff bundles
#'   diff_rna <- run_diff(project$experiments$rna, ...)
#'   diff_prot <- run_diff(project$experiments$prot, ...)
#'   con_bundle <- run_integration(
#'     project, method = "concordance",
#'     experiments = c("rna", "prot"),
#'     diff_bundles = list(rna = diff_rna, prot = diff_prot)
#'   )
#' }
run_integration <- function(
  project,
  method = c("correlation", "concordance", "active_pathways"),
  experiments = NULL,
  diff_bundles = NULL,
  by = "feature_symbol",
  ...
) {
  method <- match.arg(method)
  assert_list(diff_bundles, "diff_bundles", allow_null = TRUE)
  assert_string(by, "by")
  experiments <- resolve_experiment_pair(project, experiments)
  tag_a <- experiments[[1L]]
  tag_b <- experiments[[2L]]
  dots <- list(...)

  if (method %in% c("concordance", "active_pathways") && is.null(diff_bundles)) {
    stop("`diff_bundles` is required for method = '", method, "'.")
  }

  if (method == "correlation") {
    cor_method <- dots$cor_method %||% "spearman"
    p_adjust_method <- dots$p_adjust_method %||% "BH"
    min_samples <- dots$min_samples %||% 4L
    p_cutoff <- dots$p_cutoff %||% 0.05

    backend <- run_integration_correlation(
      project = project,
      experiments = experiments,
      method = cor_method,
      by = by,
      p_adjust_method = p_adjust_method,
      min_samples = min_samples,
      p_cutoff = p_cutoff
    )
    integration_df <- backend$std
    integration_raw <- NULL
    method_info <- backend$info
    method_params <- list(
      cor_method = cor_method,
      p_adjust_method = p_adjust_method,
      min_samples = min_samples,
      p_cutoff = p_cutoff
    )
  } else if (method == "concordance") {
    p_preference <- dots$p_preference %||% "adjusted"
    p_cutoff <- dots$p_cutoff %||% 0.05
    p_adjust_method <- dots$p_adjust_method %||% "BH"

    backend <- run_integration_concordance(
      project = project,
      experiments = experiments,
      diff_bundles = diff_bundles,
      by = by,
      p_preference = p_preference,
      p_cutoff = p_cutoff,
      p_adjust_method = p_adjust_method
    )
    integration_df <- backend$std
    integration_raw <- NULL
    method_info <- backend$info
    method_params <- list(
      p_preference = p_preference,
      p_cutoff = p_cutoff,
      p_adjust_method = p_adjust_method
    )
  } else {
    database <- dots$database %||% "hallmark"
    organism <- dots$organism %||% "Hs"
    p_preference <- dots$p_preference %||% "raw"
    significant <- dots$significant %||% 0.05
    geneset_filter <- dots$geneset_filter %||% c(5L, 1000L)
    merge_method <- dots$merge_method %||% "Brown"

    backend <- run_integration_active_pathways(
      project = project,
      experiments = experiments,
      diff_bundles = diff_bundles,
      database = database,
      organism = organism,
      by = by,
      p_preference = p_preference,
      significant = significant,
      geneset_filter = geneset_filter,
      merge_method = merge_method
    )
    integration_df <- backend$std
    integration_raw <- backend$raw
    method_info <- backend$info
    method_params <- list(
      database = normalize_enrich_database(database),
      organism = normalize_organism(organism),
      p_preference = p_preference,
      significant = significant,
      geneset_filter = geneset_filter,
      merge_method = merge_method
    )
  }

  check_integration_result_schema(integration_df)

  input_info <- list(
    experiments = experiments,
    omics_type = c(
      project$experiments[[tag_a]]$omics_type,
      project$experiments[[tag_b]]$omics_type
    ),
    assay_type = c(
      project$experiments[[tag_a]]$assay_type,
      project$experiments[[tag_b]]$assay_type
    ),
    n_features = c(
      nrow(project$experiments[[tag_a]]$expr_mat),
      nrow(project$experiments[[tag_b]]$expr_mat)
    ),
    n_samples = c(
      ncol(project$experiments[[tag_a]]$expr_mat),
      ncol(project$experiments[[tag_b]]$expr_mat)
    )
  )

  params <- c(
    list(
      method = method,
      experiments = experiments,
      by = by
    ),
    method_params,
    list(method_info = method_info)
  )

  results <- list(integration_df = integration_df)
  if (!is.null(integration_raw)) {
    results$integration_raw <- integration_raw
  }

  new_analysis_bundle(
    analysis_name = "run_integration",
    input_info = input_info,
    params = params,
    results = results
  )
}
