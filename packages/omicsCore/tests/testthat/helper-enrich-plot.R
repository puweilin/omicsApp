# Small standardized enrichment tables for the plot tests
# (test-plot-enrichment-dot.R, test-plot-enrichment-bar.R). The ORA
# default mirrors the realistic project: one pathway at -log10 p = 248,
# the rest between 1.5 and 15.

ora_df <- function(p = c(1e-248, 4e-16, 5e-12, 2e-9, 2e-9, 4e-3, 1e-2, 3e-2),
                   gene_set_size = c(152, 156, 64, 151, 143, 74, 150, 111),
                   overlap_size = c(151, 37, 21, 28, 27, 11, 16, 12),
                   direction = "up") {
  n <- length(p)
  data.frame(
    database = "hallmark", result_type = "ora", comparison = "B_vs_A",
    pathway_id = paste0("P", seq_len(n)), pathway_name = paste("pathway", seq_len(n)),
    effect = NA_real_, effect_type = NA_character_,
    direction = rep_len(direction, n), p_value = p, adj_p_value = p,
    q_value = NA_real_, gene_set_size = gene_set_size, overlap_size = overlap_size,
    overlap_features = "A", leading_features = NA_character_, source_label = "ora",
    stringsAsFactors = FALSE)
}

gsea_df <- function() {
  data.frame(
    database = "hallmark", result_type = "gsea", comparison = "B_vs_A",
    pathway_id = paste0("P", 1:4), pathway_name = paste("pathway", 1:4),
    effect = c(2.2, 1.7, -1.5, -1.9), effect_type = "nes",
    direction = c("up", "up", "down", "down"),
    p_value = c(1e-9, 1e-4, 1e-2, 1e-3), adj_p_value = c(1e-8, 1e-3, 3e-2, 5e-3),
    q_value = NA_real_, gene_set_size = c(40, 20, 10, 60), overlap_size = NA_real_,
    overlap_features = NA_character_, leading_features = "A", source_label = "gsea",
    stringsAsFactors = FALSE)
}

enrich_bundle <- function(df, type = "ora") {
  new_analysis_bundle("run_enrichment",
                      params = list(type = type, direction = "separate"),
                      results = list(enrich_result_df = df))
}
