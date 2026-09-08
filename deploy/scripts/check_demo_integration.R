# Does Integration recover what make_demo_study.R put in?
#
#   Rscript check_demo_integration.R [dir]
#
# Runs the pipeline the application runs -- import, QC, a differential per
# layer, concordance -- and cross-tabulates the answer against the classes
# the generator assigned before any numbers existed. A correct run puts
# every designed class on its own row and nothing off the diagonal.
#
# `is_significant` is not the column to filter on, and the name is why
# this is worth saying: in the concordance table it means "significant in
# both layers AND in the same direction", so every discordant pair is
# excluded by definition. Filtering on it makes a working integration look
# like one that finds no discordance at all. The per-layer p-values are
# what "significant in both" means here.

suppressPackageStartupMessages({
  library(omicsCore)
})

args <- commandArgs(trailingOnly = TRUE)
dir <- if (length(args) >= 1L) args[[1L]] else "."
rna_path <- file.path(dir, "demo_rnaseq.xlsx")
prot_path <- file.path(dir, "demo_proteomics.xlsx")
truth_path <- file.path(dir, "demo_truth.csv")
for (f in c(rna_path, prot_path, truth_path)) {
  if (!file.exists(f)) stop("Missing ", f, ". Run make_demo_study.R first.")
}

rna <- read_omics(rna_path, omics_type = "rnaseq",
                  assay_type = "raw_count")$input
prt <- read_omics(prot_path, omics_type = "proteomics",
                  assay_type = "normalized_intensity")$input
cat(sprintf("\nimported   RNA %d x %d   protein %d x %d\n",
            nrow(rna$expr_mat), ncol(rna$expr_mat),
            nrow(prt$expr_mat), ncol(prt$expr_mat)))
cat(sprintf("symbols    RNA %d mapped from Ensembl   protein %d from the file\n",
            sum(!is.na(rna$feature_df$feature_symbol)),
            sum(!is.na(prt$feature_df$feature_symbol))))

qr <- run_qc(rna)$results$cleaned_input
qp <- run_qc(prt)$results$cleaned_input
cat(sprintf("after QC   RNA %d genes (low-expression filter)   protein %d\n",
            nrow(qr$expr_mat), nrow(qp$expr_mat)))

dr <- run_diff(qr, method = "deseq2", analysis_type = "group",
               group_col = "condition", control_group = "G1",
               case_group = "G2")
dp <- run_diff(qp, method = "limma", analysis_type = "group",
               group_col = "condition", control_group = "G1",
               case_group = "G2")

project <- omics_project("demo", experiments = list(rnaseq = qr,
                                                    proteomics = qp))
pairing <- sample_pairing_preview(project, "rnaseq", "proteomics")
cat(sprintf("pairing    %d pairs, by %s\n", nrow(pairing$pairs),
            pairing$source))
if (!identical(pairing$source, "donor")) {
  cat("  ! expected pairing by donor; the sample ids are deliberately\n",
      "    different between layers, so anything else is a fallback\n")
}

ig <- run_integration(project, method = "concordance",
                      experiments = c("rnaseq", "proteomics"),
                      diff_bundles = list(rnaseq = dr, proteomics = dp))
df <- ig$results$integration_df
truth <- utils::read.csv(truth_path, stringsAsFactors = FALSE)
m <- merge(df, truth, by.x = "feature_symbol", by.y = "gene_symbol")

# Significant in both layers -- not `is_significant`, see the header.
both <- m[!is.na(m$p_a) & !is.na(m$p_b) & m$p_a < 0.05 & m$p_b < 0.05, ]
cat(sprintf("paired     %d shared symbols, %d significant in both layers\n\n",
            nrow(m), nrow(both)))

tab <- table(designed = both$class, reported = both$direction)
print(tab)

designed <- c(concordant_up = "concordant", concordant_down = "concordant",
              discordant = "discordant")
wrong <- 0L
for (cls in names(designed)) {
  if (!cls %in% rownames(tab)) next
  right <- designed[[cls]]
  for (rep_ in colnames(tab)) {
    if (!identical(rep_, right)) wrong <- wrong + tab[cls, rep_]
  }
}
n_designed <- sum(tab[intersect(names(designed), rownames(tab)), ])
cat(sprintf("\n%d of %d designed pairs reported the wrong direction.\n",
            wrong, n_designed))
cat("Rows for `null`, `rna_only` and `protein_only` are noise reaching\n",
    "significance in both layers by chance; a few are expected.\n")

# A non-zero exit, so this can gate a deployment rather than only inform
# one. Printing a mismatch count and returning 0 is a check that passes
# while reporting failure, which is the shape of check nobody notices has
# broken.
#
# Three ways this run is not evidence that integration works:
failures <- character(0)
if (!identical(pairing$source, "donor")) {
  failures <- c(failures, sprintf(
    "pairing fell back to '%s'; the sample ids differ between layers by design, so donor is the only correct answer",
    pairing$source))
}
if (wrong > 0L) {
  failures <- c(failures, sprintf("%d designed pairs reported the wrong direction",
                                  wrong))
}
# Coverage: a run where almost nothing reached significance in both
# layers agrees with the truth table trivially, having tested almost
# nothing. Two thirds of the designed changes is well inside what this
# generator produces and well outside what a broken pipeline would.
expected_designed <- sum(truth$expected_concordance != "not comparable")
if (n_designed < 0.5 * expected_designed) {
  failures <- c(failures, sprintf(
    "only %d of %d designed changes reached significance in both layers; too few to conclude anything",
    n_designed, expected_designed))
}

if (length(failures) > 0L) {
  cat("\nFAILED:\n")
  cat(paste0("  - ", failures, collapse = "\n"), "\n\n")
  quit(status = 1L)
}
cat("\nPASSED\n\n")
