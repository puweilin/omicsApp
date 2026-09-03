#' Built-in example data fixtures
#'
#' These builders produce small but realistic `omicsCore` objects so the
#' Phase 2 UI views have something to render before a user has uploaded
#' any data. All builders are deterministic (seeded), package-internal,
#' and rebuilt on demand — we keep them as functions rather than
#' shipping `.rda` data to avoid the `data-raw/` maintenance burden for
#' fixtures this small.
#'
#' Slices 2C–2F call these from inside view modules wrapped in
#' [shiny::reactive()] so each session pays the build cost once.
#'
#' @keywords internal
#' @name example-data
#' @noRd
NULL

#' Build a demo `omics_input`
#'
#' Proteomics: 50 features x 12 samples (6 G1 + 6 G2), log2-scale
#' Gaussian noise with signal injected on the first 10 features so
#' downstream diff / volcano views have something to point at.
#' RNA-seq: 60 features x 12 samples, Poisson counts with 3x
#' fold-change on the first 10 features.
#'
#' Metadata carries `group`, `age`, `sex`, `donor_id`, `batch` —
#' richer than the omicsCore test fixtures so the project / QC views
#' can demo group-by, covariate adjustment, and batch indicators.
#'
#' @param omics_type One of `"proteomics"`, `"rnaseq"`.
#'
#' @return An `omics_input` object.
#'
#' @keywords internal
#' @noRd
example_input <- function(omics_type = c("proteomics", "rnaseq")) {
  omics_type <- match.arg(omics_type)

  n_per_group <- 6L
  n_samples   <- 2L * n_per_group
  samp_ids    <- sprintf("S%02d", seq_len(n_samples))
  groups      <- rep(c("G1", "G2"), each = n_per_group)

  # Deterministic age stratified by group: G1 younger, G2 older with overlap.
  set.seed(2026L)
  ages <- c(round(stats::rnorm(n_per_group, mean = 32, sd = 6)),
            round(stats::rnorm(n_per_group, mean = 48, sd = 7)))
  sex     <- rep(c("F", "M"), length.out = n_samples)
  batch   <- rep(c("A", "B"), times = n_per_group)
  donor   <- sprintf("D%03d", seq_len(n_samples))

  meta <- data.frame(
    group     = groups,
    age       = ages,
    sex       = sex,
    donor_id  = donor,
    batch     = batch,
    row.names = samp_ids,
    stringsAsFactors = FALSE
  )

  if (omics_type == "proteomics") {
    n_features <- 50L
    feat_ids   <- sprintf("P%03d", seq_len(n_features))
    symbols    <- example_protein_symbols(n_features)
    descriptions <- paste0("Demo protein ", symbols)
    assay_type <- "normalized_intensity"

    set.seed(2026L + 1L)
    expr <- matrix(
      stats::rnorm(n_features * n_samples, mean = 18, sd = 0.7),
      nrow = n_features,
      dimnames = list(feat_ids, samp_ids)
    )
    # Inject signal: features 1-5 up in G2, 6-10 down in G2.
    # Effect sizes are deliberately wider than within-group noise so
    # the demo volcano in slice 2D clearly highlights the seeded hits.
    case <- groups == "G2"
    expr[1:5,  case] <- expr[1:5,  case] + 3.5
    expr[6:10, case] <- expr[6:10, case] - 3.0
  } else {
    n_features <- 60L
    feat_ids   <- sprintf("ENSG%011d", seq_len(n_features))
    symbols    <- example_gene_symbols(n_features)
    descriptions <- paste0("Demo gene ", symbols)
    assay_type <- "raw_count"

    set.seed(2026L + 2L)
    expr <- matrix(
      stats::rpois(n_features * n_samples, lambda = 200),
      nrow = n_features,
      dimnames = list(feat_ids, samp_ids)
    )
    case <- groups == "G2"
    expr[1:5,  case] <- expr[1:5,  case] * 4L
    expr[6:10, case] <- pmax(round(expr[6:10, case] / 4L), 1L)
  }

  feat <- data.frame(
    feature_id      = feat_ids,
    feature_symbol  = symbols,
    description     = descriptions,
    row.names       = feat_ids,
    stringsAsFactors = FALSE
  )

  omicsCore::omics_input(
    expr_mat   = expr,
    meta_df    = meta,
    feature_df = feat,
    omics_type = omics_type,
    assay_type = assay_type
  )
}

# Private cache for the demo fixtures, shared by every demo view.
.example_cache <- new.env(parent = emptyenv())

# The demo views draw through the same `omicsCore::plot_*()` dispatchers
# as the live ones, which take bundles rather than data frames. Wrapping
# the fixtures keeps both paths on one set of drawing code -- when the
# demo has its own, the two drift and the demo stops resembling what
# users get.
example_enrich_bundle <- function() {
  if (!is.null(.example_cache$enrich)) return(.example_cache$enrich)
  .example_cache$enrich <- omicsCore::new_analysis_bundle(
    analysis_name = "run_enrichment",
    input_info = list(omics_type = "proteomics"),
    params = list(type = "gsea", database = "hallmark", organism = "Hs",
                  direction = "both", comparison = "G2_vs_G1"),
    results = list(enrich_result_df = example_enrich_table())
  )
}

example_integration_bundle <- function() {
  if (!is.null(.example_cache$integration)) return(.example_cache$integration)
  .example_cache$integration <- omicsCore::new_analysis_bundle(
    analysis_name = "run_integration",
    input_info = list(omics_type = c("rnaseq", "proteomics")),
    params = list(method = "concordance",
                  experiments = c("rnaseq", "proteomics"),
                  by = "feature_symbol"),
    results = list(
      integration_df = example_integration_tables()$concordance_df)
  )
}

# Cached accessor for the demo proteomics input. Used by both
# `example_diff_bundle()` and `example_qc_bundle()` (and the diff
# view's contrast UI when no project is loaded) so all demo views
# operate on the same synthetic dataset — otherwise each call to
# `example_input("proteomics")` reseeded and produced different
# matrices, leaving the QC and Diff demos visibly inconsistent.
example_proteomics_input <- function() {
  if (!is.null(.example_cache$proteomics_input)) {
    return(.example_cache$proteomics_input)
  }
  .example_cache$proteomics_input <- example_input("proteomics")
}

#' Build a demo `analysis_bundle` (limma diff, G2 vs G1, age-adjusted)
#'
#' Runs the same contrast that the Phase 1 validation script exercises
#' (`scripts/validation/cheek_g2_vs_g1_omicscore.R`) so any regression
#' in the limma path that escapes omicsCore's own tests will also surface
#' here.
#'
#' @return An `analysis_bundle` from [omicsCore::run_diff()].
#'
#' @keywords internal
#' @noRd
example_diff_bundle <- function() {
  if (!is.null(.example_cache$diff)) return(.example_cache$diff)
  .example_cache$diff <- omicsCore::run_diff(
    example_proteomics_input(),
    method         = "limma",
    analysis_type  = "group",
    group_col      = "group",
    control_group  = "G1",
    case_group     = "G2",
    covariates     = "age"
  )
}

#' Build a demo `omics_project` with both proteomics and RNA-seq
#'
#' Used by the Project view (slice 2C) to demonstrate the dual-omics
#' container without requiring a real user upload.
#'
#' @return An `omics_project` object.
#'
#' @keywords internal
#' @noRd
example_project <- function() {
  omicsCore::omics_project(
    name = "Cheek \u00B7 G2 vs G1",
    experiments = list(
      proteomics = example_proteomics_input(),
      rnaseq     = example_input("rnaseq")
    )
  )
}

#' Build a demo enrichment-result `data.frame`
#'
#' Hand-rolled 15-row GSEA-style enrichment table over MSigDB Hallmark
#' names. Matches the standard `enrich_result_df` schema used by
#' [omicsCore::run_enrichment()] so the Enrich view (slice 2E) can
#' render a dotplot and table without dragging Bioconductor deps
#' (`clusterProfiler`, `msigdbr`) into the demo path.
#'
#' Effect (NES) is a mix of positive and negative; p-values are
#' deliberately tiny so every row clears default thresholds.
#'
#' @return A `data.frame` carrying the standard enrichment columns.
#'
#' @keywords internal
#' @noRd
example_enrich_table <- function() {
  pathways <- c(
    "INFLAMMATORY RESPONSE",
    "TNFA SIGNALING VIA NFKB",
    "IL6 JAK STAT3 SIGNALING",
    "EPITHELIAL MESENCHYMAL TRANSITION",
    "COMPLEMENT",
    "OXIDATIVE PHOSPHORYLATION",
    "MYOGENESIS",
    "APICAL JUNCTION",
    "HYPOXIA",
    "APOPTOSIS",
    "INTERFERON GAMMA RESPONSE",
    "KRAS SIGNALING UP",
    "ADIPOGENESIS",
    "Fatty acid metabolism",
    "MYC targets v1"
  )
  ids <- c(
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_IL6_JAK_STAT3_SIGNALING",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_COMPLEMENT",
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    "HALLMARK_MYOGENESIS",
    "HALLMARK_APICAL_JUNCTION",
    "HALLMARK_HYPOXIA",
    "HALLMARK_APOPTOSIS",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_KRAS_SIGNALING_UP",
    "HALLMARK_ADIPOGENESIS",
    "HALLMARK_FATTY_ACID_METABOLISM",
    "HALLMARK_MYC_TARGETS_V1"
  )
  set.seed(2026L + 4L)
  n <- length(pathways)
  nes <- c(2.31, 2.18, 2.04, 1.92, 1.81, -1.74, -1.61, -1.52,
           1.45, 1.38, 1.92, 1.55, -1.42, -1.30, 1.20)
  p_value     <- 10 ^ -seq(from = 8, by = -0.35, length.out = n)
  adj_p_value <- p_value * 1.5
  gene_set_size <- as.integer(c(200, 200, 87, 200, 200, 200, 200, 200,
                                200, 161, 200, 200, 200, 158, 200))
  overlap_size  <- as.integer(c(52, 48, 31, 46, 42, 38, 35, 33,
                                39, 36, 41, 28, 30, 22, 26))
  data.frame(
    database      = "hallmark",
    result_type   = "gsea",
    comparison    = "G2_vs_G1",
    pathway_id    = ids,
    pathway_name  = pathways,
    effect        = nes,
    effect_type   = "nes",
    direction     = ifelse(nes >= 0, "up", "down"),
    p_value       = p_value,
    adj_p_value   = adj_p_value,
    q_value       = adj_p_value,
    gene_set_size = gene_set_size,
    overlap_size  = overlap_size,
    overlap_features = NA_character_,
    leading_features = NA_character_,
    source_label  = "gsea_hallmark",
    stringsAsFactors = FALSE
  )
}

#' Build demo integration tables (concordance + ActivePathways)
#'
#' Hand-rolled fixtures for the Integration view (slice 2E). Avoids
#' the live `omicsCore::run_integration()` calls so the demo path
#' doesn't depend on `ActivePathways` (Bioconductor Suggests).
#'
#' Returns a list:
#'
#' * `concordance_df` — feature-level data frame with effect sizes from
#'   both omics layers (`effect_a` = RNA log2FC, `effect_b` = protein
#'   log2FC), their difference, a p-value pair, and the four-quadrant
#'   concordance label (`up_up`, `down_down`, `up_down`, `down_up`,
#'   `ns`). ~60 rows so the dual-volcano and scatter cards have
#'   something to render.
#' * `active_pathways_df` — pathway-level data frame with per-omics
#'   p-values (`p_a` = protein, `p_b` = RNA) and a combined p-value,
#'   matching the six pathways shown in the mockup. All rows are
#'   `direction = "shared"` since the mockup only displays shared
#'   findings.
#'
#' @return Named list with `concordance_df` and `active_pathways_df`.
#'
#' @keywords internal
#' @noRd
example_integration_tables <- function() {
  set.seed(2026L + 5L)
  n_feat <- 60L
  symbols <- example_protein_symbols(n_feat)
  feat_ids <- sprintf("F%03d", seq_len(n_feat))

  # Concordant signal on the first 30 features (split up/down), the
  # remaining 30 are noise around zero so the scatter shows a clear
  # diagonal cloud plus an "ns" centre.
  rna   <- c(stats::rnorm(15, mean =  1.8, sd = 0.5),
             stats::rnorm(15, mean = -1.6, sd = 0.5),
             stats::rnorm(30, mean =  0.0, sd = 0.45))
  prot  <- c(rna[1:15]  + stats::rnorm(15, mean =  0.2, sd = 0.4),
             rna[16:30] + stats::rnorm(15, mean = -0.2, sd = 0.4),
             stats::rnorm(30, mean =  0.0, sd = 0.5))

  effect_diff <- prot - rna
  p_value     <- 10 ^ -stats::runif(n_feat, min = 1.0, max = 8.0)
  # Strengthen the seeded signal so significance follows the design.
  p_value[1:30] <- p_value[1:30] / 1e3
  adj_p_value   <- pmin(p_value * 1.4, 1)

  thr <- 0.5  # log2FC threshold for quadrant labels
  quadrant <- ifelse(rna >  thr & prot >  thr, "up_up",
              ifelse(rna < -thr & prot < -thr, "down_down",
              ifelse(rna >  thr & prot < -thr, "up_down",
              ifelse(rna < -thr & prot >  thr, "down_up", "ns"))))

  concordance_df <- data.frame(
    feature_id     = feat_ids,
    feature_symbol = symbols,
    effect_a       = rna,
    effect_b       = prot,
    effect_diff    = effect_diff,
    p_value        = p_value,
    adj_p_value    = adj_p_value,
    quadrant       = quadrant,
    # The columns below carry no extra information -- they exist so the
    # fixture satisfies INTEGRATION_RESULT_REQUIRED_COLS and can be
    # plotted by `omicsCore::plot_integration()`. Without them the demo
    # would need its own drawing code, and a demo drawn by different
    # code from the live view is a demo that can quietly stop matching
    # what users actually get.
    result_type    = "concordance",
    experiments    = "rnaseq vs proteomics",
    comparison     = "G2_vs_G1",
    # The schema defines `effect` as effect_a - effect_b.
    effect         = rna - prot,
    effect_type    = "log2fc_difference",
    statistic      = rna - prot,
    statistic_type = "difference",
    direction      = ifelse(rna - prot >= 0, "up", "down"),
    is_significant = adj_p_value < 0.05,
    source_label   = "demo fixture",
    stringsAsFactors = FALSE
  )

  active_pathways_df <- data.frame(
    pathway_id = c(
      "R-HSA-IL6_SIGNALING",
      "R-HSA-SASP",
      "R-HSA-ECM_ORGANIZATION",
      "R-HSA-COLLAGEN_BIOSYNTHESIS",
      "R-HSA-INNATE_IMMUNE_SYSTEM",
      "R-HSA-OX_PHOS"
    ),
    pathway_name = c(
      "Interleukin-6 signaling",
      "Senescence-associated secretory phenotype",
      "Extracellular matrix organization",
      "Collagen biosynthesis",
      "Innate immune system",
      "OXIDATIVE PHOSPHORYLATION"
    ),
    p_a        = c(1.2e-12, 8.6e-09, 2.4e-08, 5.8e-08, 3.1e-06, 3.5e-06),
    p_b        = c(4.1e-08, 2.0e-07, 7.8e-06, 1.2e-05, 8.4e-05, 2.1e-04),
    p_combined = c(3.4e-18, 2.1e-14, 9.5e-12, 3.0e-11, 1.6e-08, 4.2e-08),
    direction  = "shared",
    stringsAsFactors = FALSE
  )

  list(
    concordance_df     = concordance_df,
    active_pathways_df = active_pathways_df
  )
}

#' Build a demo QC `analysis_bundle`
#'
#' Returns `omicsCore::run_qc()` on the proteomics demo input with
#' ~5% of cells set to NA so the QC view (slice 2D) has something
#' visible in the missingness panel. Outlier detection is forced to
#' `"iqr"` so we don't depend on the optional WGCNA suggest.
#'
#' @return An `analysis_bundle` from [omicsCore::run_qc()].
#'
#' @keywords internal
#' @noRd
example_qc_bundle <- function() {
  if (!is.null(.example_cache$qc)) return(.example_cache$qc)
  .example_cache$qc <- omicsCore::run_qc(
    example_qc_input(),
    missing_threshold = 0.5,
    outlier_method    = "iqr"
  )
}

#' The input behind the demo QC bundle
#'
#' Split out from [example_qc_bundle()] so the QC view can run the demo
#' through `run_qc()` with the *live* slider and radio values. Handing
#' back a bundle computed at fixed settings left both controls visibly
#' enabled and doing nothing, which reads as a broken app rather than as
#' a demo.
#'
#' 50 x 12 proteomics, so a re-run costs milliseconds and needs no cache
#' of its own.
#'
#' @return An `omics_input` with ~5% of cells set to `NA`.
#' @keywords internal
#' @noRd
example_qc_input <- function() {
  if (!is.null(.example_cache$qc_input)) return(.example_cache$qc_input)
  input <- example_proteomics_input()
  expr  <- input$expr_mat
  set.seed(2026L + 3L)
  n_cells <- length(expr)
  na_idx  <- sample.int(n_cells, size = ceiling(0.05 * n_cells))
  expr[na_idx] <- NA_real_
  input$expr_mat <- expr
  .example_cache$qc_input <- input
}

# ---- internal symbol pools -------------------------------------------

# Recognisable-looking symbols for the demo. Recycled if `n` exceeds the
# pool. Picked from common skin-aging proteomics / RNA-seq markers so
# the demo looks plausible without being a real dataset.
example_protein_symbols <- function(n) {
  pool <- c("COL1A1","COL1A2","COL3A1","FBN1","ELN","LOX","MMP1","MMP2",
            "MMP9","TIMP1","TIMP2","DCN","BGN","LUM","DPT","ACTA2",
            "VIM","KRT5","KRT14","KRT10","KRT1","FLG","LOR","IVL",
            "DSG1","DSG3","DSC1","CDSN","CASP14","SPRR2A","S100A8",
            "S100A9","DEFB1","DEFB4A","LCE3D","TGM1","TGM3","CDH1",
            "ITGB1","ITGA6","PECAM1","CD31","VWF","COL4A1","COL4A2",
            "LAMB1","LAMA5","NID1","HSPG2","PERLECAN")
  rep_len(pool, n)
}

example_gene_symbols <- function(n) {
  pool <- c("IFNG","IL6","IL10","IL1B","TNF","NFKB1","STAT1","STAT3",
            "JAK1","JAK2","CXCL10","CXCL9","CCL2","CCL5","FOXP3",
            "TGFB1","BMP2","WNT3A","CTNNB1","MYC","TP53","CDKN1A",
            "MKI67","PCNA","ESR1","AR","NR3C1","PPARG","PPARGC1A",
            "SIRT1","FOXO1","FOXO3","MTOR","AKT1","PIK3CA","RPS6",
            "EIF4E","HIF1A","VEGFA","ANGPT1","ANGPT2","KDR","FLT1",
            "ICAM1","VCAM1","SELP","SELE","CDH5","TJP1","CLDN1",
            "OCLN","ZO1","CDKN2A","RB1","BRCA1","CCND1","CCNE1",
            "CDK4","CDK6","E2F1")
  rep_len(pool, n)
}
