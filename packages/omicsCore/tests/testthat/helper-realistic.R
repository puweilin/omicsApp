# A fixture with real gene symbols.
#
# Most fixtures in this suite name their features gene_1, GENE07, P0004.
# That is fine for testing shapes and errors, and useless for testing
# anything that looks a gene up: no pathway database contains gene_1, so
# every enrichment over such a fixture is an enrichment of nothing, and
# a test that only checks the class of what comes back cannot tell a
# working GSEA from one that failed inside a tryCatch. It happened.
#
# The symbols below are members of four MSigDB Hallmark sets, chosen
# because those sets are old, large and stable. Signal is injected into
# one of them, so a test can ask the question a user would: did the
# analysis find the pathway that was actually perturbed?

REAL_GENE_SETS <- list(
  G2M = c(
    "ABL1", "AMD1", "ARID4A", "ATF5", "ATRX", "AURKA", "AURKB", "BARD1",
    "BCL3", "BIRC5", "BRCA2", "BUB1", "BUB3", "CASP8AP2", "CBX1", "CCNA2",
    "CCNB2", "CCND1", "CCNF", "CCNT1", "CDC20", "CDC25A", "CDC25B", "CDC27",
    "CDC45", "CDC6", "CDC7", "CDK1", "CDK4", "CDKN1B", "CDKN3", "CENPA",
    "CENPE", "CENPF", "CHAF1A", "CHEK1", "CHMP1A", "CKS1B", "CKS2", "CTCF"
  ),
  IFNG = c(
    "ADAR", "APOL6", "ARID5B", "AUTS2", "B2M", "BANK1", "BATF2", "BPGM",
    "BST2", "BTG1", "C1R", "C1S", "CASP1", "CASP3", "CASP4", "CASP8",
    "CCL2", "CCL5", "CCL7", "CD274", "CD38", "CD40", "CD69", "CD74",
    "CD86", "CDKN1A", "CFB", "CFH", "CIITA", "CMKLR1", "CMPK2", "CMTR1",
    "CSF2RB", "CXCL10", "CXCL11", "CXCL9", "DDX60", "DHX58", "EIF2AK2", "FAS"
  ),
  OXPHOS = c(
    "ABCB7", "ACAA1", "ACADSB", "ACADVL", "ACAT1", "AFG3L2", "ALAS1",
    "ALDH6A1", "ATP1B1", "ATP5F1A", "ATP5F1B", "ATP5F1C", "ATP5F1D",
    "ATP5F1E", "ATP5MC1", "ATP5MC2", "ATP5MC3", "ATP5ME", "ATP5MF", "ATP5MG",
    "ATP5PB", "ATP5PD", "ATP5PF", "ATP6AP1", "ATP6V0B", "ATP6V0C", "ATP6V0E1",
    "ATP6V1C1", "ATP6V1D", "ATP6V1E1", "ATP6V1F", "ATP6V1G1", "ATP6V1H",
    "BAX", "BDH2", "COX10", "COX11", "COX15", "ACADM", "ACO2"
  ),
  ADIPO = c(
    "ABCA1", "ABCB8", "ACADL", "ACADS", "ACLY", "ACOX1", "ADCY6", "ADIG",
    "ADIPOQ", "ADIPOR2", "AGPAT3", "AK2", "ALDH2", "ALDOA", "ANGPT1",
    "ANGPTL4", "APLP2", "APOE", "ARAF", "ATL2", "ATP1B3", "BAZ2A", "BCL2L13",
    "BCL6", "C3", "CAT", "CAVIN1", "CAVIN2", "CCNG2", "CD151", "CD302",
    "CD36", "CHCHD10", "CHUK", "CIDEA", "CMBL", "CMPK1", "ACAA2", "AIFM1",
    "BCKDHA"
  )
)

# Which Hallmark set each block above comes from, as run_enrichment()
# labels it: the HALLMARK_ prefix stripped and underscores replaced.
REAL_GENE_SET_LABELS <- c(
  G2M    = "G2M CHECKPOINT",
  IFNG   = "INTERFERON GAMMA RESPONSE",
  OXPHOS = "OXIDATIVE PHOSPHORYLATION",
  ADIPO  = "ADIPOGENESIS"
)

realistic_symbols <- function() {
  unique(unlist(REAL_GENE_SETS, use.names = FALSE))
}

#' An omics_input whose features are real genes
#'
#' @param omics_type "proteomics" (log2 intensities) or "rnaseq" (counts).
#' @param n_per_group Samples per group; groups are G1 (control) and G2.
#' @param signal Name of the block in REAL_GENE_SETS that is shifted up in
#'   G2, or NULL for a null dataset.
#' @param effect Shift on the log2 scale.
#' @param seed Seed, so two calls with the same arguments agree exactly.
realistic_input <- function(omics_type = c("proteomics", "rnaseq"),
                            n_per_group = 6L,
                            signal = "G2M",
                            effect = 1.5,
                            seed = 2027L) {
  omics_type <- match.arg(omics_type)
  set.seed(seed)
  symbols <- realistic_symbols()
  n_feat <- length(symbols)
  n_samp <- 2L * n_per_group
  sample_ids <- sprintf("S%02d", seq_len(n_samp))
  group <- rep(c("G1", "G2"), each = n_per_group)
  is_case <- group == "G2"

  shift <- matrix(0, n_feat, n_samp)
  if (!is.null(signal)) {
    hit <- symbols %in% REAL_GENE_SETS[[signal]]
    shift[hit, is_case] <- effect
  }

  if (omics_type == "proteomics") {
    feature_ids <- paste0("PROT_", symbols)
    base <- stats::rnorm(n_feat, mean = 20, sd = 2)
    mat <- matrix(stats::rnorm(n_feat * n_samp, sd = 0.6), n_feat, n_samp) +
      base + shift
    assay_type <- "normalized_intensity"
  } else {
    feature_ids <- paste0("ENSG", formatC(seq_len(n_feat), width = 11,
                                          flag = "0"))
    mu <- exp(stats::rnorm(n_feat, mean = log(400), sd = 1))
    mu_mat <- mu * 2^shift
    mat <- matrix(stats::rnbinom(n_feat * n_samp, mu = mu_mat, size = 15),
                  n_feat, n_samp)
    storage.mode(mat) <- "integer"
    assay_type <- "raw_count"
  }
  dimnames(mat) <- list(feature_ids, sample_ids)

  meta <- data.frame(
    group = group,
    age = round(stats::runif(n_samp, 25, 70)),
    sex = rep(c("F", "M"), length.out = n_samp),
    batch = rep(c("A", "B"), times = n_per_group),
    row.names = sample_ids,
    stringsAsFactors = FALSE
  )
  feat <- data.frame(
    feature_id = feature_ids,
    feature_symbol = symbols,
    stringsAsFactors = FALSE
  )
  omics_input(mat, meta, feat, omics_type = omics_type,
              assay_type = assay_type)
}

realistic_diff_bundle <- function(method = "limma", ...) {
  run_diff(realistic_input(...), method = method, analysis_type = "group",
           group_col = "group", control_group = "G1", case_group = "G2")
}
