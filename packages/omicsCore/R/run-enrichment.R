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
#' For ORA, `direction` decides which gene list is tested. The default,
#' `"separate"`, tests the up-regulated and the down-regulated hits as two
#' lists, each with its own multiple-testing correction, and returns both
#' in one table whose `direction` column says which list a pathway came
#' from -- the same rows two runs with `"up"` and `"down"` would give.
#' `"both"` pools up and down into one list (the behaviour before
#' `"separate"` existed, and what saved results and scripts with
#' `direction = "both"` still mean); pooling mixes opposite biology, so a
#' pathway half up and half down can look enriched.
#'
#' For GSEA the full ranked vector is always passed to the backend, so up
#' and down are always separate; `"up"` or `"down"` keep one sign only,
#' and `"separate"` and `"both"` keep both.
#'
#' A result with no direction (an ANOVA or spline fit) can only be
#' enriched as one pooled list: ORA then defaults to `"both"`, and asking
#' for a direction is an error.
#'
#' When several databases are queried, `p_adjust_scope` decides the
#' family each pathway is corrected in: `"database"` (default) adjusts
#' within each database, as running them one at a time would; `"all"`
#' adjusts across every database's pathways together (within each ORA
#' gene list), the stricter choice when the databases are read as one
#' search. Under `"all"`, `q_value` is set to `NA` (it is a per-database
#' quantity), and the clusterProfiler objects in `enrich_object` keep
#' their per-database adjustment -- the table is the result.
#'
#' Gene sets come in the chosen species' own symbols (see
#' [enrichment_species()]). When the data's symbols match them poorly as
#' written but well ignoring case (a mouse table in upper case, say), the
#' match ignoring case is used, the result's gene lists keep the data's
#' spelling, and the bundle carries a warning saying so
#' (`params$symbol_case` is then `"ignored"`).
#'
#' For GSVA-style sample-level scoring, see [run_gsva()].
#'
#' @param diff_bundle An [`analysis_bundle`][is_analysis_bundle()] produced
#'   by [run_diff()].
#' @param type One of `"ora"` or `"gsea"`.
#' @param database One of `"hallmark"`, `"kegg"`, `"reactome"`,
#'   `"wikipathways"`, `"go_bp"`, `"go_mf"`, `"go_cc"`. Pass a character
#'   vector to query multiple databases.
#' @param organism Species of the gene sets: a code or name from
#'   [enrichment_species()] (`"Hs"`, `"Mm"`, `"Rn"`, `"Dr"`, `"Dm"`,
#'   `"Sc"`, `"Ce"`, or e.g. `"mouse"`), or any full species name
#'   `msigdbr::msigdbr_species()` lists.
#' @param direction One of `"separate"` (default: up and down tested as
#'   two lists), `"up"`, `"down"`, or `"both"` (up and down pooled into one
#'   list). See Details.
#' @param p_cutoff Significance cutoff for selecting diff features (ORA
#'   only; GSEA ranks the whole list).
#' @param output_p_cutoff Bound on the returned table. Defaults to
#'   `p_cutoff`. Pass `1` to keep every pathway, so a caller can choose
#'   raw or adjusted p at display time without re-running.
#' @param p_preference For ORA feature selection: `"adjusted"` (default) or
#'   `"raw"`.
#' @param effect_cutoff Optional |effect| cutoff for ORA feature selection.
#' @param p_adjust_method Multiple-testing correction method.
#' @param p_adjust_scope `"database"` (default) to correct within each
#'   database, or `"all"` to correct across all queried databases
#'   together. Makes a difference only with more than one database.
#' @param min_size,max_size Min/max gene-set sizes, for ORA and GSEA alike.
#' @param ... Reserved for backend-specific extensions.
#'
#' @return An [`analysis_bundle`][is_analysis_bundle()] with
#'   `results$enrich_result_df` (standardized schema) and
#'   `results$enrich_object` (named list of clusterProfiler objects keyed by
#'   database, or for ORA by `<list>__<database>` where `<list>` is `up`,
#'   `down` or `both`).
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
  direction = c("separate", "up", "down", "both"),
  p_cutoff = 0.05,
  output_p_cutoff = NULL,
  p_preference = c("adjusted", "raw"),
  effect_cutoff = NULL,
  p_adjust_method = "BH",
  p_adjust_scope = c("database", "all"),
  min_size = 10L,
  max_size = 500L,
  ...
) {
  if (!is_analysis_bundle(diff_bundle) ||
      !identical(diff_bundle$analysis_name, "run_diff")) {
    stop("`diff_bundle` must be an analysis_bundle from run_diff().")
  }

  # A run with several comparisons stacks them in one table; enriching
  # that table would pool the genes of every comparison (ORA) or keep
  # each gene's largest effect across them (GSEA) -- a gene list that
  # belongs to no comparison at all.
  n_cmp <- length(unique(stats::na.omit(diff_bundle$results$diff_result_df$comparison)))
  if (n_cmp > 1L) {
    stop("`diff_bundle` holds ", n_cmp, " comparisons. Take one out with ",
         "select_comparison(), or enrich them all side by side with ",
         "compare_enrichment().", call. = FALSE)
  }

  type <- match.arg(type)
  direction_given <- !missing(direction)
  direction <- match.arg(direction)
  p_preference <- match.arg(p_preference)
  p_adjust_scope <- match.arg(p_adjust_scope)

  # A global test (ANOVA, LRT) or a spline fit says *whether* a feature
  # changes, not which way. GSEA on it ranked by an unsigned score, so
  # "down" meant "least variable"; ORA "up"/"down" split nothing.
  at <- unique(stats::na.omit(diff_bundle$results$diff_result_df$analysis_type))
  undirected <- any(at %in% c("anova", "continuous_spline"))
  # There is only one list to test, so the default is to test it; only a
  # direction actually asked for is refused.
  if (undirected && type == "ora" && !direction_given) direction <- "both"
  if (undirected && (type == "gsea" || direction != "both")) {
    stop("This differential result (", paste(at, collapse = ", "), ") has no ",
         "direction: it says which features change, not which way. Use ORA ",
         "with direction = \"both\", or enrich a pairwise comparison.",
         call. = FALSE)
  }

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

  # Corrected across databases: every pathway the backends tested is
  # kept, so the correction sees the whole family, and the bound is
  # applied after it.
  cross_db <- identical(p_adjust_scope, "all") && length(databases) > 1L
  final_cutoff <- output_p_cutoff %||% p_cutoff
  backend_cutoff <- if (cross_db) 1 else output_p_cutoff

  warns <- character(0)
  case_fix <- match_symbol_case(diff_bundle, databases, organism)
  if (!is.null(case_fix)) {
    diff_bundle <- case_fix$bundle
    warns <- c(warns, case_fix$note)
  }

  per_db <- lapply(databases, function(db) {
    if (type == "ora") {
      run_ora_from_bundle(
        diff_bundle = diff_bundle,
        database = db,
        organism = organism,
        direction = direction,
        p_cutoff = p_cutoff,
        output_p_cutoff = backend_cutoff,
        p_preference = p_preference,
        effect_cutoff = effect_cutoff,
        p_adjust_method = p_adjust_method,
        min_size = min_size,
        max_size = max_size
      )
    } else {
      run_gsea_from_bundle(
        diff_bundle = diff_bundle,
        database = db,
        organism = organism,
        # GSEA always keeps the two signs apart, so "separate" and the
        # older "both" ask the same of it.
        direction = if (direction == "separate") "both" else direction,
        p_cutoff = p_cutoff,
        output_p_cutoff = backend_cutoff,
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
  if (cross_db) {
    enrich_result_df <- adjust_enrich_across_databases(
      enrich_result_df, type = type, p_adjust_method = p_adjust_method,
      cutoff = final_cutoff)
  }
  if (!is.null(case_fix) && nrow(enrich_result_df)) {
    enrich_result_df$overlap_features <-
      restore_symbol_case(enrich_result_df$overlap_features, case_fix$back)
    enrich_result_df$leading_features <-
      restore_symbol_case(enrich_result_df$leading_features, case_fix$back)
  }
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
      output_p_cutoff = output_p_cutoff,
      p_preference = p_preference,
      effect_cutoff = effect_cutoff,
      p_adjust_method = p_adjust_method,
      p_adjust_scope = p_adjust_scope,
      min_size = min_size,
      max_size = max_size,
      symbol_case = if (is.null(case_fix)) "exact" else "ignored",
      rank_metric = if (type == "gsea") gsea_rank_metric(diff_bundle),
      geneset_sources = geneset_sources,
      comparison = diff_bundle$params$comparison
    ),
    results = list(
      enrich_result_df = enrich_result_df,
      enrich_object = enrich_object
    ),
    warnings = warns
  )
}

# The case-insensitive fallback of `run_enrichment()` (see the notes on
# symbol_case_plan()). Returns NULL when the symbols match as written --
# the common case -- or when the gene sets cannot be read here (the
# backend then reports that itself); otherwise the bundle with its
# symbols respelled, the map back to the data's spelling, and the note
# for the user.
match_symbol_case <- function(diff_bundle, databases, organism) {
  df <- diff_bundle$results$diff_result_df
  if (!"feature_symbol" %in% names(df)) return(NULL)
  if (!is_installed("clusterProfiler") || !is_installed("msigdbr")) return(NULL)
  reference <- tryCatch(
    unique(unlist(lapply(databases, function(db) {
      unique(build_term_tables(database = db, organism = organism)$term2gene$gene)
    }), use.names = FALSE)),
    error = function(e) NULL
  )
  if (!length(reference)) return(NULL)
  plan <- symbol_case_plan(df$feature_symbol, reference)
  if (!isTRUE(plan$apply)) return(NULL)
  respelled <- apply_symbol_case(df$feature_symbol, plan)
  changed <- !is.na(respelled) & respelled != df$feature_symbol
  back <- df$feature_symbol[changed]
  names(back) <- respelled[changed]
  back <- back[!duplicated(names(back))]
  example <- if (any(changed)) {
    sprintf(" (for example %s was read as %s)",
            df$feature_symbol[changed][[1L]], respelled[changed][[1L]])
  } else ""
  df$feature_symbol <- respelled
  diff_bundle$results$diff_result_df <- df
  list(
    bundle = diff_bundle,
    back = back,
    note = sprintf(paste(
      "The gene names are written differently from the %s gene sets: only %d",
      "of %d matched as written, %d when upper and lower case are ignored.",
      "They were matched ignoring case%s; check that %s is the right species."),
      organism, plan$n_exact, plan$n_symbols, plan$n_ci, example, organism)
  )
}
