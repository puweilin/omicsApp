# GSVA (Gene Set Variation Analysis) entry point. Returns a pathways-by-
# samples score matrix; sample-level association testing (gsva_diff) is
# deferred to a later slice to avoid a circular dependency on run_diff.

#' Run Gene Set Variation Analysis
#'
#' Computes sample-level pathway activity scores via the `GSVA` package. The
#' caller may either pass a named list of gene sets directly (`gene_sets`)
#' or pick one of the supported MSigDB databases (`database`). RNA-seq raw
#' counts are log2-transformed before scoring.
#'
#' @param input A validated `omics_input`.
#' @param gene_sets Optional named list of character vectors. When supplied,
#'   `database` is ignored.
#' @param database Database key when `gene_sets` is `NULL`. See
#'   [list_gene_sets()] for supported keys.
#' @param organism Organism shorthand (e.g. `"Hs"`).
#' @param method One of `"gsva"` (default) or `"ssgsea"`.
#' @param min_size,max_size Min / max gene-set size after intersecting with
#'   the assay's feature set.
#' @param kcdf Kernel CDF used by `GSVA::gsvaParam()` — typically
#'   `"Gaussian"` for log-scale data, `"Poisson"` for raw counts.
#' @param ... Reserved for future arguments.
#'
#' @return An [`analysis_bundle`][is_analysis_bundle()] with
#'   `results$gsva_matrix` (pathways × samples) and
#'   `results$gsva_gene_sets` (the gene-set list actually used).
#' @export
#' @family enrich
#' @examples
#' \dontrun{
#'   gsva_bundle <- run_gsva(input, database = "hallmark")
#'   head(gsva_bundle$results$gsva_matrix[, 1:4])
#' }
run_gsva <- function(
  input,
  gene_sets = NULL,
  database = "hallmark",
  organism = "Hs",
  method = c("gsva", "ssgsea"),
  min_size = 10L,
  max_size = 500L,
  kcdf = NULL,
  ...
) {
  validate_omics_input(input)
  method <- match.arg(method)
  assert_list(gene_sets, "gene_sets", allow_null = TRUE)
  assert_character(database, "database", allow_null = TRUE)
  assert_count(min_size, "min_size", lower = 1L)
  assert_count(max_size, "max_size", lower = 1L)
  assert_choice(kcdf, "kcdf", c("Gaussian", "Poisson", "none"), allow_null = TRUE)

  if (!is_installed("GSVA")) {
    stop(
      "Package 'GSVA' is required for run_gsva(). ",
      "Install with: omicsCore::install_optional('enrichment').",
      call. = FALSE
    )
  }

  expr_mat <- as.matrix(input$expr_mat)
  if (identical(input$assay_type, "raw_count")) {
    expr_mat <- log2(expr_mat + 1)
  }

  # Replace expression rownames with feature_symbol so the gene-set lookup
  # works against HGNC symbols (the space msigdbr returns).
  #
  # Matched on the feature_id column, not on rownames(feature_df): the
  # constructor does not set those, so a feature_df built by hand has
  # 1..n for row names, and indexing it by feature id returned NA for
  # every feature. Every row was then dropped and GSVA refused the empty
  # matrix -- and the tests skipped over the refusal as an environment
  # problem.
  if (!is.null(input$feature_df) && "feature_symbol" %in% colnames(input$feature_df)) {
    idx <- match(rownames(expr_mat), input$feature_df$feature_id)
    syms <- as.character(input$feature_df$feature_symbol[idx])
    valid <- !is.na(syms) & nzchar(syms)
    expr_mat <- expr_mat[valid, , drop = FALSE]
    syms <- syms[valid]
    expr_mat <- expr_mat[!duplicated(syms), , drop = FALSE]
    rownames(expr_mat) <- syms[!duplicated(syms)]
  }

  resolved_database <- NA_character_
  if (is.null(gene_sets)) {
    ensure_enrichment_deps()
    resolved_database <- normalize_enrich_database(database)
    gene_sets <- get_gene_set_list(
      database = resolved_database,
      organism = organism,
      min_size = min_size,
      max_size = max_size
    )
  } else {
    if (!is.list(gene_sets) || is.null(names(gene_sets))) {
      stop("`gene_sets` must be a named list of character vectors.")
    }
  }

  if (length(gene_sets) == 0L) {
    stop("No gene sets available after filtering (min_size=", min_size,
         ", max_size=", max_size, ").")
  }

  if (is.null(kcdf)) {
    kcdf <- "Gaussian"
  }

  param_obj <- switch(
    method,
    gsva = GSVA::gsvaParam(
      exprData = expr_mat,
      geneSets = gene_sets,
      minSize = min_size,
      maxSize = max_size,
      kcdf = kcdf,
      maxDiff = TRUE
    ),
    ssgsea = GSVA::ssgseaParam(
      exprData = expr_mat,
      geneSets = gene_sets,
      minSize = min_size,
      maxSize = max_size
    )
  )

  gsva_mat <- GSVA::gsva(param_obj, verbose = FALSE)

  new_analysis_bundle(
    analysis_name = "run_gsva",
    input_info = list(
      omics_type = input$omics_type,
      assay_type = input$assay_type,
      n_samples = ncol(input$expr_mat),
      n_features = nrow(input$expr_mat)
    ),
    params = list(
      method = method,
      database = resolved_database,
      organism = normalize_organism(organism),
      min_size = min_size,
      max_size = max_size,
      kcdf = kcdf,
      gene_sets_supplied = !is.null(database) && is.na(resolved_database) == FALSE
    ),
    results = list(
      gsva_matrix = gsva_mat,
      gsva_gene_sets = gene_sets
    )
  )
}
