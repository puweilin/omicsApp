# A two-layer study with a known answer, for checking Integration.
#
#   Rscript make_demo_study.R [outdir] [n_donors]
#
# Writes three files:
#
#   demo_rnaseq.xlsx       raw counts, one sheet per the RNA-seq template
#   demo_proteomics.xlsx   intensities, per the proteomics template
#   demo_truth.csv         what each shared gene was built to do
#
# The point is the third file. A demo where every gene is random tells you
# the application runs; it cannot tell you the concordance is right,
# because there is nothing to be right about. Here each shared symbol is
# assigned a class before any numbers exist, and the classes are chosen to
# separate the ways this can be wrong:
#
#   concordant_up     up in both layers
#   concordant_down   down in both
#   discordant        up in RNA, down in protein -- the case a sign error
#                     or a swapped contrast turns into a concordant one
#   rna_only          moves in RNA, flat in protein
#   protein_only      flat in RNA, moves in protein
#   null              flat in both
#
# The two layers deliberately use different sample ids and the same
# donors, as the templates do and as a real study does: one person, two
# tissues, two sample names. Integration has to pair them through the
# `donor` column, and the ids are shuffled so pairing by position would
# give the wrong answer rather than accidentally the right one.

suppressPackageStartupMessages({
  library(omicsCore)
})
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  stop("openxlsx is required to write the workbooks.")
}

args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1L) args[[1L]] else "."
n_donors <- if (length(args) >= 2L) as.integer(args[[2L]]) else 60L
if (is.na(n_donors) || n_donors < 6L || n_donors %% 2L != 0L) {
  stop("n_donors must be an even number, at least 6.")
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260908)

# ---- who -------------------------------------------------------------
donors <- sprintf("D%03d", seq_len(n_donors))
condition <- rep(c("G1", "G2"), each = n_donors / 2L)
age <- round(c(stats::rnorm(n_donors / 2L, 27, 4),
               stats::rnorm(n_donors / 2L, 58, 5)))
sex <- sample(c("F", "M"), n_donors, replace = TRUE)

# ---- which genes -----------------------------------------------------
#
# Real Ensembl ids, because the RNA-seq layer has no feature sheet: the
# symbol is derived by mapping the id against the bundled HGNC table, and
# a made-up ENSG maps to nothing. Integration then joins on the symbol,
# so an invented id is a gene that silently cannot be paired.
hgnc <- readRDS(system.file("extdata", "hgnc_ensembl.rds",
                            package = "omicsCore"))
hgnc <- hgnc[hgnc$locus_group == "protein-coding gene", ]
hgnc <- hgnc[!duplicated(hgnc$symbol) & !duplicated(hgnc$ensembl_gene_id), ]

n_shared <- 900L
n_rna_only_genes <- 3000L        # measured in RNA and not in protein
n_prot_only_genes <- 300L
stopifnot(nrow(hgnc) >= n_shared + n_rna_only_genes + n_prot_only_genes)

pick <- hgnc[sample(nrow(hgnc), n_shared + n_rna_only_genes + n_prot_only_genes), ]
shared <- pick[seq_len(n_shared), ]
rna_extra <- pick[n_shared + seq_len(n_rna_only_genes), ]
prot_extra <- pick[n_shared + n_rna_only_genes + seq_len(n_prot_only_genes), ]

# ---- the answer, decided first ---------------------------------------
classes <- c(concordant_up = 90L, concordant_down = 90L, discordant = 60L,
             rna_only = 80L, protein_only = 80L)
classes <- c(classes, null = n_shared - sum(classes))
truth <- data.frame(
  gene_symbol = shared$symbol,
  ensembl_gene_id = shared$ensembl_gene_id,
  class = rep(names(classes), times = classes),
  stringsAsFactors = FALSE)
truth <- truth[sample(nrow(truth)), ]
rownames(truth) <- NULL

# log2 fold changes the generator will impose. Protein follows RNA at
# about half the size, which is what a real pair looks like -- the
# concordance is in the sign and the ranking, not in the magnitude.
lfc <- function(cls, layer) {
  up <- 1.4; down <- -1.4
  prot <- function(x) x * 0.75
  switch(cls,
    concordant_up   = if (layer == "rna") up else prot(up),
    concordant_down = if (layer == "rna") down else prot(down),
    discordant      = if (layer == "rna") up else prot(down),
    rna_only        = if (layer == "rna") up else 0,
    protein_only    = if (layer == "rna") 0 else prot(up),
    0)
}
truth$rna_log2fc <- vapply(truth$class, lfc, numeric(1), layer = "rna")
truth$protein_log2fc <- vapply(truth$class, lfc, numeric(1), layer = "protein")
truth$expected_concordance <- ifelse(
  truth$rna_log2fc == 0 | truth$protein_log2fc == 0, "not comparable",
  ifelse(sign(truth$rna_log2fc) == sign(truth$protein_log2fc),
         "concordant", "discordant"))

# ---- RNA-seq: counts -------------------------------------------------
rna_genes <- rbind(shared[, c("ensembl_gene_id", "symbol")],
                   rna_extra[, c("ensembl_gene_id", "symbol")])
rna_lfc <- c(truth$rna_log2fc[match(shared$symbol, truth$gene_symbol)],
             rep(0, nrow(rna_extra)))

# A spread of baselines, and a tail of genes nobody detected -- the QC
# low-expression filter exists for those, and a demo without them makes
# it look like a control that does nothing.
base_mu <- c(stats::rlnorm(nrow(rna_genes) - 600L, log(60), 1.4),
             rep(0, 400L), stats::runif(200L, 0, 3))
base_mu <- sample(base_mu)

rna_counts <- matrix(0L, nrow(rna_genes), n_donors,
                     dimnames = list(rna_genes$ensembl_gene_id, NULL))
# `if`, not `ifelse`: ifelse() returns a value shaped like its *test*,
# and the test here is one string. `ifelse(TRUE, 2^rna_lfc, 1)` is
# `2^rna_lfc[1]` -- one number applied to every gene, which after
# library-size normalisation is not a fold change at all but an invisible
# constant. The first version of this file did that, and the demo it
# produced showed no signal in RNA and the same shift on every protein.
for (j in seq_len(n_donors)) {
  mu <- base_mu * if (condition[[j]] == "G2") 2^rna_lfc else 1
  rna_counts[, j] <- stats::rnbinom(nrow(rna_genes), mu = mu, size = 8)
}
rna_counts[is.na(rna_counts)] <- 0L

# ---- proteomics: intensities ----------------------------------------
prot_genes <- rbind(shared[, c("ensembl_gene_id", "symbol")],
                    prot_extra[, c("ensembl_gene_id", "symbol")])
prot_lfc <- c(truth$protein_log2fc[match(shared$symbol, truth$gene_symbol)],
              rep(0, nrow(prot_extra)))
prot_base <- stats::rnorm(nrow(prot_genes), 22, 1.6)

prot_mat <- matrix(NA_real_, nrow(prot_genes), n_donors)
for (j in seq_len(n_donors)) {
  shift <- if (condition[[j]] == "G2") prot_lfc else 0
  prot_mat[, j] <- prot_base + shift + stats::rnorm(nrow(prot_genes), 0, 0.35)
}
# Left-censored missingness, which is what the imputation control is for:
# the lowest-abundance proteins are the ones that go missing. Kept modest
# on purpose -- MinProb fills a gap near the minimum, so on an up-regulated
# protein it pulls the estimate back down, and enough of it turns a signal
# this file put in into one nothing can find. The effect is meant to be
# recoverable; the point of the demo is checking the concordance, not
# defeating it.
p_missing <- stats::plogis((stats::median(prot_base) - prot_base) * 1.1) * 0.22
for (i in seq_len(nrow(prot_mat))) {
  drop <- stats::runif(n_donors) < p_missing[[i]]
  prot_mat[i, drop] <- NA_real_
}
# A UniProt-shaped accession, since that is what this layer's ids are.
prot_ids <- sprintf("Q%05d", seq_len(nrow(prot_genes)))
rownames(prot_mat) <- prot_ids

# ---- ids: same people, different sample names, shuffled --------------
#
# Shuffled independently per layer so the two files do not line up
# column-by-column. Pairing by position would then be visibly wrong,
# rather than right by luck -- which is the failure this demo is for.
rna_order <- sample(n_donors)
prot_order <- sample(n_donors)
rna_ids <- sprintf("RNA%03d", seq_len(n_donors))
prot_ids_s <- sprintf("PRT%03d", seq_len(n_donors))

rna_expr <- data.frame(feature_id = rownames(rna_counts),
                       rna_counts[, rna_order, drop = FALSE],
                       check.names = FALSE, stringsAsFactors = FALSE)
names(rna_expr)[-1L] <- rna_ids
rna_info <- data.frame(sample_id = rna_ids,
                       donor = donors[rna_order],
                       condition = condition[rna_order],
                       age = age[rna_order], sex = sex[rna_order],
                       stringsAsFactors = FALSE)

prot_expr <- data.frame(feature_id = rownames(prot_mat),
                        round(prot_mat[, prot_order, drop = FALSE], 3),
                        check.names = FALSE, stringsAsFactors = FALSE)
names(prot_expr)[-1L] <- prot_ids_s
prot_info <- data.frame(sample_id = prot_ids_s,
                        donor = donors[prot_order],
                        condition = condition[prot_order],
                        age = age[prot_order], sex = sex[prot_order],
                        stringsAsFactors = FALSE)
prot_feat <- data.frame(feature_id = prot_ids,
                        gene_symbol = prot_genes$symbol,
                        protein_description = paste(prot_genes$symbol,
                                                    "protein"),
                        stringsAsFactors = FALSE)

# ---- write ------------------------------------------------------------
rna_path <- file.path(outdir, "demo_rnaseq.xlsx")
prot_path <- file.path(outdir, "demo_proteomics.xlsx")
truth_path <- file.path(outdir, "demo_truth.csv")

openxlsx::write.xlsx(
  list(expression = rna_expr, sample_info = rna_info),
  file = rna_path, overwrite = TRUE)
openxlsx::write.xlsx(
  list(expression = prot_expr, sample_info = prot_info,
       feature_info = prot_feat),
  file = prot_path, overwrite = TRUE)
utils::write.csv(truth, truth_path, row.names = FALSE)

cat(sprintf("\n%d donors, %d per group\n", n_donors, n_donors / 2L))
cat(sprintf("  %-24s %5d genes x %d samples\n", basename(rna_path),
            nrow(rna_expr), n_donors))
cat(sprintf("  %-24s %5d proteins x %d samples, %.0f%% missing\n",
            basename(prot_path), nrow(prot_expr), n_donors,
            100 * mean(is.na(prot_mat))))
cat(sprintf("  %-24s %5d shared symbols\n", basename(truth_path), nrow(truth)))
cat("\nbuilt to be:\n")
print(table(truth$class))
cat("\nof the shared genes that move in both layers:\n")
print(table(truth$expected_concordance[
  truth$expected_concordance != "not comparable"]))
cat("\n")
