# The tutorial: an example project a first-time user can take through
# every step of the app, and the page furniture that walks them through
# it.
#
# The demo that each view falls back on when nothing is loaded is a
# picture, not a project: it cannot be saved, integrated or reported,
# and its two layers share no genes. This one is a real `omics_project`
# built so that each step has something to find --
#
#   * two layers, proteomics (log intensities) and RNA-seq (raw counts),
#     measuring the same genes, so Integration has features to join;
#   * three groups -- Control and two treatments -- so the Differential
#     view shows the several-groups-against-one-control design;
#   * different sample ids in the two layers but a shared `donor`
#     column, so the pairing card has a real pairing to show;
#   * genes from real MSigDB Hallmark sets carrying the signal, so
#     Enrichment finds the pathways that were put there:
#       TreatA raises INFLAMMATORY_RESPONSE genes in both layers,
#       TreatB lowers OXIDATIVE_PHOSPHORYLATION genes in both layers,
#       TreatB raises E2F_TARGETS in the RNA only (a discordant slice).

TUTORIAL_GENES <- list(
  inflammatory = c(
    "ABCA1", "ABI1", "ACVR1B", "ACVR2A", "ADGRE1", "ADM", "ADORA2B",
    "ADRM1", "AHR", "APLNR", "AQP9", "ATP2A2", "ATP2B1", "ATP2C1", "AXL",
    "BDKRB1", "BEST1", "BST2", "BTG2", "C3AR1", "C5AR1", "CALCRL",
    "CCL17", "CCL2", "CCL20", "CCL22", "CCL24", "CCL5", "CCL7", "CCR7",
    "CCRL2", "CD14", "CD40", "CD48", "CD55", "CD69" 
  ),
  oxphos = c(
    "ABCB7", "ACAA1", "ACAA2", "ACADM", "ACADSB", "ACADVL", "ACAT1",
    "ACO2", "AFG3L2", "AIFM1", "ALAS1", "ALDH6A1", "ATP1B1", "ATP5F1A",
    "ATP5F1B", "ATP5F1C", "ATP5F1D", "ATP5F1E", "ATP5MC1", "ATP5MC2",
    "ATP5MC3", "ATP5ME", "ATP5MF", "ATP5MG", "ATP5PB", "ATP5PD",
    "ATP5PF", "ATP5PO", "ATP6AP1", "ATP6V0B", "ATP6V0C", "ATP6V0E1",
    "ATP6V1C1", "ATP6V1D", "ATP6V1E1", "ATP6V1F" 
  ),
  e2f = c(
    "AK2", "ANP32E", "ASF1A", "ASF1B", "ATAD2", "AURKA", "AURKB",
    "BARD1", "BIRC5", "BRCA1", "BRCA2", "BRMS1L", "BUB1B", "CBX5",
    "CCNB2", "CCNE1", "CCP110", "CDC20", "CDC25A", "CDC25B", "CDCA3",
    "CDCA8", "CDK1", "CDK4", "CDKN1A", "CDKN1B", "CDKN2A", "CDKN2C",
    "CDKN3", "CENPE" 
  ),
  background = c(
    "ABCB8", "ABCF2", "ACACA", "ACADL", "ACADS", "ACLY", "ACOX1",
    "ACSL3", "ACTA1", "ACTB", "ACTC1", "ACTG1", "ACTG2", "ACTN1",
    "ACTN2", "ACTN3", "ACTN4", "ACTR2", "ACTR3", "ADAM10", "ADAM15",
    "ADAM23", "ADAM9", "ADAMTS5", "ADCY6", "ADD3", "ADIG", "ADIPOQ",
    "ADIPOR2", "ADRA1B", "AGPAT3", "AK4", "AKT2", "AKT3", "ALDH2",
    "ALDOA", "ALOX15B", "AMH", "AMIGO1", "AMIGO2", "ANGPT1", "ANGPTL4",
    "AP1G1", "AP2B1", "AP2M1", "AP2S1", "AP3B1", "AP3S1", "APLP2",
    "APOE", "ARAF", "ARCN1", "ARF1", "ARFGAP3", "ARFGEF1", "ARFGEF2",
    "ARFIP1", "ARHGEF6", "ARL4A", "ARPC2", "ARPC5L", "ASNS", "ATL2",
    "ATP1A1", "ATP1A3", "ATP1B3", "ATP6V1B1", "ATP6V1H", "ATP7A",
    "B4GALT1", "BAIAP2", "BAZ2A", "BCAT1", "BCKDHA", "BCL2L13", "BCL6",
    "BET1", "BHLHE40", "BMP1", "BNIP3", "BUB1", "C3", "CACYBP", "CADM2",
    "CADM3", "CALB2", "CALR", "CANX", "CAP1", "CAT", "CAV2", "CAVIN1",
    "CAVIN2", "CCNF", "CCNG1", "CCNG2", "CCT6A", "CD151", "CD209",
    "CD274", "CD276", "CD302", "CD34", "CD36", "CD63", "CD86", "CD9",
    "CD99", "CDH1", "CDH11", "CDH15", "CDH3", "CDH4", "CDH6", "CDH8",
    "CDK8", "CDSN", "CERCAM", "CFP", "CHCHD10", "CHUK", "CIDEA", "CLCN3",
    "CLDN11", "CLDN14", "CLDN15", "CLDN18", "CLDN19", "CLDN4", "CLDN5",
    "CLDN6", "CLDN7", "CLDN8", "CLDN9", "CLN5", "CLTA", "CLTC", "CMBL",
    "CMPK1", "CNN2", "CNTN1", "COG2", "COL15A1", "COL16A1", "COL17A1",
    "COL4A1", "COL9A1", "COPB1", "COPB2", "COPE" 
  )
)

#' Build the tutorial project
#'
#' Deterministic: seeded through [local_seed()], so it never touches the
#' session's random stream and every user sees the same numbers.
#'
#' @return An `omics_project` with layers `proteomics` and `rnaseq`.
#' @keywords internal
#' @noRd
tutorial_project <- function() {
  if (!is.null(.example_cache$tutorial)) return(.example_cache$tutorial)
  local_seed(20260930L)

  genes <- unlist(TUTORIAL_GENES, use.names = FALSE)
  set_of <- rep(names(TUTORIAL_GENES), lengths(TUTORIAL_GENES))
  n_genes <- length(genes)

  groups <- rep(c("Control", "TreatA", "TreatB"), each = 4L)
  n <- length(groups)
  donors <- sprintf("D%02d", seq_len(n))
  meta <- function(prefix) {
    ids <- sprintf("%s%02d", prefix, seq_len(n))
    data.frame(
      group = groups,
      donor = donors,
      age   = c(34, 41, 29, 52, 38, 45, 31, 49, 36, 43, 30, 55),
      sex   = rep(c("F", "M"), length.out = n),
      batch = rep(c("B1", "B2"), times = n / 2L),
      row.names = ids,
      stringsAsFactors = FALSE
    )
  }

  # Log2 effect of each group on each gene, shared by both layers except
  # where noted. A per-donor offset makes the layers of one person
  # resemble each other, which is what sample-level correlation reads.
  effect <- matrix(0, n_genes, n)
  effect[set_of == "inflammatory", groups == "TreatA"] <- 1.6
  effect[set_of == "oxphos", groups == "TreatB"] <- -1.3
  donor_offset <- matrix(stats::rnorm(n_genes * n, sd = 0.35), n_genes, n)

  # Proteomics: every protein, log2 intensity.
  base_p <- stats::rnorm(n_genes, mean = 24, sd = 1.5)
  prot <- base_p + effect + donor_offset +
    matrix(stats::rnorm(n_genes * n, sd = 0.3), n_genes, n)
  pid <- sprintf("PROT%04d", seq_len(n_genes))
  pmeta <- meta("P")
  dimnames(prot) <- list(pid, rownames(pmeta))
  # A few proteins not detected in every sample, as in real DIA data.
  miss <- sample.int(length(prot), size = round(0.02 * length(prot)))
  prot[miss] <- NA_real_

  # RNA-seq: the same genes as counts, with E2F targets up in TreatB on
  # the RNA side only.
  effect_r <- effect
  effect_r[set_of == "e2f", groups == "TreatB"] <- 1.4
  base_r <- stats::runif(n_genes, log2(80), log2(4000))
  lib <- stats::runif(n, 0.8, 1.25)
  mu <- 2^(base_r + effect_r + donor_offset)
  mu <- sweep(mu, 2L, lib, "*")
  counts <- matrix(stats::rnbinom(length(mu), mu = mu, size = 25),
                   n_genes, n)
  gid <- sprintf("ENSG%011d", 100000L + seq_len(n_genes))
  rmeta <- meta("R")
  dimnames(counts) <- list(gid, rownames(rmeta))

  proteomics <- omicsCore::omics_input(
    expr_mat   = prot,
    meta_df    = pmeta,
    feature_df = data.frame(feature_id = pid, feature_symbol = genes,
                            row.names = pid, stringsAsFactors = FALSE),
    omics_type = "proteomics",
    assay_type = "normalized_intensity"
  )
  rnaseq <- omicsCore::omics_input(
    expr_mat   = counts,
    meta_df    = rmeta,
    feature_df = data.frame(feature_id = gid, feature_symbol = genes,
                            row.names = gid, stringsAsFactors = FALSE),
    omics_type = "rnaseq",
    assay_type = "raw_count"
  )
  # As the import view would record it: the groups, and the control.
  proteomics <- omicsCore::set_study_design(proteomics, "group", "Control")
  rnaseq <- omicsCore::set_study_design(rnaseq, "group", "Control")
  proj <- omicsCore::omics_project(
    name = "Tutorial \u00B7 two treatments vs control",
    experiments = list(proteomics = proteomics, rnaseq = rnaseq)
  )
  proj$bundles <- list()
  # Nothing visited yet, so the Workflow card walks through every step.
  proj$visited_steps <- character(0)
  .example_cache$tutorial <- proj
}
