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
#'   `ActivePathways` package, fed by the p-values and fold-change signs
#'   of two differential analyses.
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
#' * `p_preference` -- `"adjusted"` (default) or `"raw"`: which p-value of
#'   each layer is held to `p_cutoff` when deciding whether a feature is a
#'   hit in that layer. The combined (Fisher) p-value is always built from
#'   the raw p-values and corrected once across features.
#' * `p_cutoff` -- significance cutoff (default `0.05`).
#' * `effect_cutoff` -- minimum absolute effect in each layer for a
#'   feature to count as a hit there (default `0`, no bound).
#' * `p_adjust_method` -- defaults to `"BH"`.
#'
#' Besides the schema columns, the concordance table keeps what each layer
#' said on its own: `effect_a`, `effect_b`, `p_value_a`, `p_value_b`,
#' `adj_p_value_a`, `adj_p_value_b`, `significant_a`, `significant_b`,
#' `feature_id_a`, `feature_id_b`. Features are matched on `by`
#' ignoring case and surrounding whitespace; where several features of
#' one layer share a symbol the most abundant (`base_mean`) is kept.
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
#' * `merge_method` -- `"DPM"` (default), the directional extension of
#'   Brown's method, or another method accepted by
#'   `ActivePathways::ActivePathways()` (`"Brown"` ignores direction, as
#'   results made before the directional method did).
#' * `constraints_vector` -- for a directional `merge_method`, the expected
#'   sign relation of the two layers' fold changes. Default `c(1, 1)`: the
#'   layers are expected to move the same way (RNA and protein of a gene
#'   rise together), so genes whose layers disagree are penalised and
#'   contribute less to any pathway. `c(1, -1)` would expect them to move
#'   oppositely.
#'
#' With a directional merge, the ActivePathways table also says which way
#' each pathway went: `direction` is `"up"` or `"down"` when the genes that
#' drove it (ActivePathways' overlap) moved that way in both layers, and
#' `"mixed"` when the layers disagree; `direction_a` / `direction_b` give
#' each layer's own reading (the sign most driving genes share in it),
#' `layers_agree` whether those match, and `n_genes_agree` /
#' `n_genes_disagree` count the driving genes measured in both layers
#' whose signs agree or not. `evidence` says which layers found the
#' pathway on their own: `"shared"` (both), `"unique"` (one), or
#' `"combined"` (only once merged). Results made before this held the
#' evidence in `direction`. The direction is reported for a non-directional
#' merge too, but there it describes the result without having shaped it.
#' When the installed ActivePathways predates the directional method (it
#' needs 2.0 or later), Brown's method is used instead, the bundle's
#' `params$merge_method` records `"Brown"`, and a warning says so.
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
    effect_cutoff <- dots$effect_cutoff %||% 0
    assert_number(effect_cutoff, "effect_cutoff", lower = 0)

    backend <- run_integration_concordance(
      project = project,
      experiments = experiments,
      diff_bundles = diff_bundles,
      by = by,
      p_preference = p_preference,
      p_cutoff = p_cutoff,
      p_adjust_method = p_adjust_method,
      effect_cutoff = effect_cutoff
    )
    integration_df <- backend$std
    integration_raw <- NULL
    method_info <- backend$info
    method_params <- list(
      p_preference = p_preference,
      p_cutoff = p_cutoff,
      effect_cutoff = effect_cutoff,
      p_adjust_method = p_adjust_method
    )
  } else {
    database <- dots$database %||% "hallmark"
    organism <- dots$organism %||% "Hs"
    p_preference <- dots$p_preference %||% "raw"
    significant <- dots$significant %||% 0.05
    geneset_filter <- dots$geneset_filter %||% c(5L, 1000L)
    merge_method <- dots$merge_method %||% "DPM"
    constraints_vector <- dots$constraints_vector %||% c(1, 1)

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
      merge_method = merge_method,
      constraints_vector = constraints_vector
    )
    integration_df <- backend$std
    integration_raw <- backend$raw
    method_info <- backend$info
    # The method that actually ran: a directional request on an older
    # ActivePathways falls back to Brown's, and the script must repeat
    # what ran.
    method_params <- list(
      database = normalize_enrich_database(database),
      organism = normalize_organism(organism),
      p_preference = p_preference,
      significant = significant,
      geneset_filter = geneset_filter,
      merge_method = method_info$merge_method %||% merge_method,
      constraints_vector = method_info$constraints_vector
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
    list(method_info = method_info),
    # How each layer's differential result was made, so export_script()
    # can make them again: the partner layer's run happens inside the
    # app and exists nowhere else.
    if (!is.null(diff_bundles)) {
      list(diff_params = lapply(diff_bundles[intersect(experiments, names(diff_bundles))],
                                function(b) b$params))
    }
  )

  results <- list(integration_df = integration_df)
  if (!is.null(integration_raw)) {
    results$integration_raw <- integration_raw
  }

  warns <- as.character(method_info$notes %||% character(0))
  n_amb <- method_info$n_ambiguous_samples %||% 0L
  if (n_amb > 0L) {
    warns <- c(warns, sprintf(
      "%d sample(s) share a donor with another sample of the same layer and were left out of the pairing; only one sample per donor and layer is correlated.",
      n_amb))
  }

  new_analysis_bundle(
    analysis_name = "run_integration",
    input_info = input_info,
    params = params,
    results = results,
    warnings = warns
  )
}
