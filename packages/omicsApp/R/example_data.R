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
  omicsCore::run_diff(
    example_input("proteomics"),
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
    name = "CHISSS demo \u00B7 Cheek \u00B7 G2 vs G1",
    experiments = list(
      proteomics = example_input("proteomics"),
      rnaseq     = example_input("rnaseq")
    )
  )
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
