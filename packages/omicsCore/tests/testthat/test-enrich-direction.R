# ORA's gene lists, the multiple-testing family across databases, and the
# species the gene sets come in.
#
# ORA used to pool the up- and down-regulated hits into one list by
# default, which mixes opposite biology: a pathway half up and half down
# looked enriched. The default now tests the two lists separately and
# returns both, with "both" kept as the pooled analysis it always meant.

skip_if_no_enrichment <- function() {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("msigdbr")
}

# Signal in two blocks pushed opposite ways, so each list has a pathway.
two_way_bundle <- function() {
  inp <- realistic_input(signal = "G2M", effect = 1.5)
  hit <- inp$feature_df$feature_symbol %in% REAL_GENE_SETS$OXPHOS
  case <- inp$meta_df$group == "G2"
  inp$expr_mat[hit, case] <- inp$expr_mat[hit, case] - 1.5
  run_diff(inp, method = "limma", analysis_type = "group",
           group_col = "group", control_group = "G1", case_group = "G2")
}

std_cols <- function(df) {
  df <- df[order(df$direction, df$pathway_id), , drop = FALSE]
  rownames(df) <- NULL
  df
}

# ---- up and down separately -------------------------------------------

test_that("ORA tests up and down separately by default", {
  skip_if_no_enrichment()
  b <- two_way_bundle()
  res <- run_enrichment(b, type = "ora", database = "hallmark")
  expect_identical(res$params$direction, "separate")
  df <- res$results$enrich_result_df
  expect_setequal(unique(df$direction), c("up", "down"))
  expect_true("G2M CHECKPOINT" %in% df$pathway_name[df$direction == "up"])
  expect_true("OXIDATIVE PHOSPHORYLATION" %in% df$pathway_name[df$direction == "down"])
  expect_false("OXIDATIVE PHOSPHORYLATION" %in% df$pathway_name[df$direction == "up"])
  expect_setequal(names(res$results$enrich_object), c("up__hallmark", "down__hallmark"))
})

test_that("separate gives exactly the rows of an up run and a down run", {
  skip_if_no_enrichment()
  b <- two_way_bundle()
  args <- list(diff_bundle = b, type = "ora", database = "hallmark",
               p_cutoff = 0.05, output_p_cutoff = 1)
  sep <- do.call(run_enrichment, c(args, direction = "separate"))
  up <- do.call(run_enrichment, c(args, direction = "up"))
  down <- do.call(run_enrichment, c(args, direction = "down"))
  expect_equal(std_cols(sep$results$enrich_result_df),
               std_cols(rbind(up$results$enrich_result_df,
                              down$results$enrich_result_df)))
  # Each list carries its own correction: the adjusted p of an up row is
  # the up run's, not one corrected over both lists.
  expect_equal(
    sep$results$enrich_result_df$adj_p_value[sep$results$enrich_result_df$direction == "up"],
    up$results$enrich_result_df$adj_p_value)
})

test_that("'both' still pools up and down into one list", {
  skip_if_no_enrichment()
  b <- two_way_bundle()
  res <- run_enrichment(b, type = "ora", database = "hallmark",
                        direction = "both", output_p_cutoff = 1)
  df <- res$results$enrich_result_df
  expect_true(all(is.na(df$direction)))
  expect_identical(names(res$results$enrich_object), "both__hallmark")
  # The same as handing clusterProfiler the pooled hit list by hand.
  rd <- b$results$diff_result_df
  hits <- filter_diff_results(rd, p_cutoff = 0.05, p_preference = "adjusted")
  universe <- unique(stats::na.omit(rd$feature_symbol[!is.na(rd$p_value)]))
  direct <- run_ora_database(unique(hits$feature_symbol), universe, "hallmark",
                             p_cutoff = 1)
  expect_setequal(res$results$enrich_object$both__hallmark@gene,
                  unique(hits$feature_symbol))
  expect_equal(df$p_value, as.data.frame(direct)$pvalue)
  expect_equal(df$adj_p_value, as.data.frame(direct)$p.adjust)
})

test_that("the comparison enrichment follows the separate default and shows the list", {
  skip_if_no_enrichment()
  inp <- realistic_input(signal = "G2M")
  meta <- inp$meta_df
  # A third group whose OXPHOS block went down.
  m3 <- inp$expr_mat[, meta$group == "G1"] +
    matrix(stats::rnorm(sum(meta$group == "G1") * nrow(inp$expr_mat), sd = 0.3),
           nrow(inp$expr_mat))
  hit <- inp$feature_df$feature_symbol %in% REAL_GENE_SETS$OXPHOS
  m3[hit, ] <- m3[hit, ] - 1.5
  colnames(m3) <- paste0("T", seq_len(ncol(m3)))
  meta3 <- meta[meta$group == "G1", , drop = FALSE]
  rownames(meta3) <- colnames(m3)
  meta3$group <- "G3"
  x <- omics_input(cbind(inp$expr_mat, m3), rbind(meta, meta3), inp$feature_df,
                   omics_type = "proteomics", assay_type = "normalized_intensity")
  b <- run_diff(x, method = "limma", group_col = "group", control_group = "G1",
                case_group = c("G2", "G3"))
  ce <- compare_enrichment(b, type = "ora", database = "hallmark")
  df <- ce$results$enrich_result_df
  expect_identical(ce$params$direction, "separate")
  expect_true(any(df$comparison == "G2_vs_G1" & df$direction == "up"))
  expect_true(any(df$comparison == "G3_vs_G1" & df$direction == "down"))
  p <- plot_enrichment_comparison(ce)
  built <- ggplot2::ggplot_build(p)
  expect_identical(p$scales$get_scales("colour")$name, "found among")
  expect_match(p$labels$subtitle, "tested separately", fixed = TRUE)
  expect_true(all(c(omics_colors$up, omics_colors$down) %in% built$data[[1L]]$colour))
})

test_that("a pathway found in both lists of one comparison is drawn as 'Up and down'", {
  df <- data.frame(
    database = "hallmark", result_type = "ora", comparison = "B_vs_A",
    pathway_id = c("P1", "P1", "P2"), pathway_name = c("one", "one", "two"),
    effect = NA_real_, effect_type = NA_character_,
    direction = c("up", "down", "down"), p_value = c(1e-4, 1e-3, 1e-3),
    adj_p_value = c(1e-3, 1e-2, 1e-2), q_value = NA_real_,
    gene_set_size = 50, overlap_size = c(10, 8, 7), overlap_features = "A",
    leading_features = "A", source_label = "ora", stringsAsFactors = FALSE)
  b <- new_analysis_bundle("compare_enrichment",
                           params = list(type = "ora", direction = "separate",
                                         database = "hallmark", comparison = "B_vs_A"),
                           results = list(enrich_result_df = df))
  p <- plot_enrichment_comparison(b)
  expect_setequal(as.character(p$data$.list), c("Up and down", "Down-regulated genes"))
})

test_that("the dot plot splits ORA's two lists into labelled panels", {
  df <- data.frame(
    database = "hallmark", result_type = "ora", comparison = "B_vs_A",
    pathway_id = paste0("P", 1:6), pathway_name = paste("pathway", 1:6),
    effect = NA_real_, effect_type = NA_character_,
    direction = c("up", "up", "up", "up", "down", "down"),
    p_value = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 2e-2),
    adj_p_value = c(1e-5, 1e-4, 1e-3, 1e-2, 3e-2, 4e-2), q_value = NA_real_,
    gene_set_size = 50, overlap_size = 5, overlap_features = "A",
    leading_features = "A", source_label = "ora", stringsAsFactors = FALSE)
  b <- new_analysis_bundle("run_enrichment",
                           params = list(type = "ora", direction = "separate"),
                           results = list(enrich_result_df = df))
  p <- plot_enrichment(b, top_n = 4L)
  # top_n is shared between the lists, alternating: the down list is not
  # crowded out by four stronger up pathways.
  expect_setequal(p$data$pathway_id, c("P1", "P2", "P5", "P6"))
  expect_setequal(as.character(p$data$.list), c("Up-regulated genes", "Down-regulated genes"))
  expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
  expect_s3_class(plot_enrichment(b, view = "bar"), "ggplot")
})

test_that("a result with no direction is enriched pooled by default, and 'separate' is refused", {
  skip_if_no_enrichment()
  inp <- realistic_input(n_per_group = 6L, signal = "G2M")
  inp$meta_df$group <- rep(c("G1", "G2", "G3"), each = 4L)
  hit <- inp$feature_df$feature_symbol %in% REAL_GENE_SETS$G2M
  inp$expr_mat[hit, inp$meta_df$group == "G3"] <-
    inp$expr_mat[hit, inp$meta_df$group == "G3"] + 2
  an <- run_diff(inp, method = "limma", analysis_type = "anova", group_col = "group")
  res <- run_enrichment(an, type = "ora", database = "hallmark")
  expect_identical(res$params$direction, "both")
  expect_true(all(is.na(res$results$enrich_result_df$direction)))
  expect_error(run_enrichment(an, type = "ora", direction = "separate"), "no direction")
})

test_that("GSEA treats 'separate' as both signs, as before", {
  skip_if_no_enrichment()
  skip_if_not_installed("fgsea")
  b <- two_way_bundle()
  sep <- suppressWarnings(run_enrichment(b, type = "gsea", database = "hallmark"))
  both <- suppressWarnings(run_enrichment(b, type = "gsea", database = "hallmark",
                                          direction = "both"))
  expect_equal(sep$results$enrich_result_df, both$results$enrich_result_df)
})

# ---- multiple testing across databases --------------------------------

test_that("adjusting across databases treats each gene list as its own family", {
  df <- data.frame(
    database = c("a", "a", "b", "b", "a", "b"),
    direction = c("up", "up", "up", "up", "down", "down"),
    p_value = c(0.001, 0.02, 0.01, 0.04, 0.003, 0.03),
    adj_p_value = 0, q_value = 0.5, stringsAsFactors = FALSE)
  out <- adjust_enrich_across_databases(df, "ora", "BH", cutoff = 1)
  up <- out$direction == "up"
  expect_equal(out$adj_p_value[up], stats::p.adjust(df$p_value[1:4], "BH"))
  expect_equal(out$adj_p_value[!up], stats::p.adjust(df$p_value[5:6], "BH"))
  expect_true(all(is.na(out$q_value)))
  # The bound is applied after the correction.
  cut <- adjust_enrich_across_databases(df, "ora", "BH", cutoff = 0.03)
  expect_true(all(cut$adj_p_value <= 0.03 & cut$p_value <= 0.03))
  # GSEA has one ranking, so one family.
  g <- adjust_enrich_across_databases(df, "gsea", "BH", cutoff = 1)
  expect_equal(g$adj_p_value, stats::p.adjust(df$p_value, "BH"))
})

test_that("p_adjust_scope = 'all' corrects across databases; the default does not", {
  skip_if_no_enrichment()
  b <- two_way_bundle()
  dbs <- c("hallmark", "kegg")
  per_db <- run_enrichment(b, type = "ora", database = dbs, output_p_cutoff = 1)
  expect_identical(per_db$params$p_adjust_scope, "database")
  alone <- run_enrichment(b, type = "ora", database = "hallmark", output_p_cutoff = 1)
  expect_equal(std_cols(per_db$results$enrich_result_df[
    per_db$results$enrich_result_df$database == "hallmark", ]),
    std_cols(alone$results$enrich_result_df))

  all_db <- run_enrichment(b, type = "ora", database = dbs, output_p_cutoff = 1,
                           p_adjust_scope = "all")
  expect_identical(all_db$params$p_adjust_scope, "all")
  df <- all_db$results$enrich_result_df
  expect_gt(length(unique(df$database)), 1L)
  for (d in c("up", "down")) {
    i <- df$direction == d
    expect_equal(df$adj_p_value[i], stats::p.adjust(df$p_value[i], "BH"))
  }
  # Raw p-values are untouched; only the correction changed.
  m <- merge(per_db$results$enrich_result_df, df,
             by = c("database", "pathway_id", "direction"))
  expect_equal(m$p_value.x, m$p_value.y)

  # With the default bound, what is kept passes it after the correction.
  bounded <- run_enrichment(b, type = "ora", database = dbs, p_adjust_scope = "all")
  kept <- bounded$results$enrich_result_df
  expect_true(all(kept$adj_p_value <= 0.05))
  expect_setequal(paste(kept$database, kept$pathway_id, kept$direction),
                  with(df[df$p_value <= 0.05 & df$adj_p_value <= 0.05, ],
                       paste(database, pathway_id, direction)))
})

test_that("p_adjust_scope makes no difference with one database", {
  skip_if_no_enrichment()
  b <- two_way_bundle()
  a <- run_enrichment(b, type = "ora", database = "hallmark")
  z <- run_enrichment(b, type = "ora", database = "hallmark", p_adjust_scope = "all")
  expect_equal(a$results$enrich_result_df, z$results$enrich_result_df)
  expect_error(run_enrichment(b, p_adjust_scope = "everything"))
})

test_that("the script repeats the direction and the correction scope", {
  inp <- realistic_input()
  b <- realistic_diff_bundle()
  eb <- new_analysis_bundle("run_enrichment", input_info = b$input_info,
                            params = list(type = "ora", database = c("hallmark", "kegg"),
                                          organism = "Rattus norvegicus",
                                          direction = "separate",
                                          p_adjust_scope = "all",
                                          symbol_case = "ignored"))
  proj <- omics_project("p", list(proteomics = inp))
  proj$bundles <- list(diff = b, enrich = eb)
  txt <- paste(export_script(proj, include_plots = FALSE), collapse = "\n")
  expect_match(txt, 'direction      = "separate"', fixed = TRUE)
  expect_match(txt, 'p_adjust_scope = "all"', fixed = TRUE)
  expect_match(txt, 'organism       = "Rattus norvegicus"', fixed = TRUE)
  # Not an argument: the run decides it again from the data.
  expect_false(grepl("symbol_case", txt, fixed = TRUE))
  # An old bundle's "both" is written as it was: pooled.
  eb$params <- list(type = "ora", database = "hallmark", direction = "both")
  proj$bundles$enrich <- eb
  txt <- paste(export_script(proj, include_plots = FALSE), collapse = "\n")
  expect_match(txt, 'direction = "both"', fixed = TRUE)
})

# ---- species ------------------------------------------------------------

test_that("the species offered resolve to msigdbr's names", {
  sp <- enrichment_species()
  expect_true(all(c("Hs", "Mm", "Rn", "Dr", "Dm", "Sc") %in% sp$code))
  for (i in seq_len(nrow(sp))) {
    expect_identical(normalize_organism(sp$code[[i]]), sp$species[[i]])
    expect_identical(normalize_organism(sp$species[[i]]), sp$species[[i]])
  }
  expect_identical(normalize_organism("rat"), "Rattus norvegicus")
  expect_identical(normalize_organism("Zebrafish"), "Danio rerio")
  expect_identical(normalize_organism("fly"), "Drosophila melanogaster")
  expect_identical(normalize_organism("yeast"), "Saccharomyces cerevisiae")
  expect_error(normalize_organism("Xx"), "Unsupported organism")
  expect_identical(organism_code("Rattus norvegicus"), "Rn")
})

test_that("other species msigdbr knows are accepted by full name", {
  skip_if_not_installed("msigdbr")
  extra <- msigdbr_species_names()
  skip_if(!"Bos taurus" %in% extra, "msigdbr does not list Bos taurus")
  expect_identical(normalize_organism("Bos taurus"), "Bos taurus")
})

# The msigdbr call, or a skip when its data are not available here (msigdbr
# 10 keeps them in msigdbdf, which may be missing offline).
gene_sets_or_skip <- function(organism) {
  skip_if_not_installed("msigdbr")
  out <- tryCatch(list_gene_sets("hallmark", organism = organism),
                  error = function(e) NULL)
  if (is.null(out) || !nrow(out)) skip(paste("no msigdbr data for", organism))
  out
}

test_that("gene sets come in each species' own symbols", {
  rat <- gene_sets_or_skip("Rn")
  fish <- gene_sets_or_skip("Dr")
  human <- gene_sets_or_skip("Hs")
  expect_true("Tp53" %in% rat$gene_symbol)
  expect_false("TP53" %in% rat$gene_symbol)
  expect_true("tp53" %in% fish$gene_symbol)
  expect_true("TP53" %in% human$gene_symbol)
  # The same pathways, carried across by orthology.
  expect_true(length(intersect(unique(rat$pathway_name), unique(human$pathway_name))) > 40L)
  yeast <- gene_sets_or_skip("Sc")
  fly <- gene_sets_or_skip("Dm")
  expect_gt(nrow(yeast), 0L)
  expect_gt(nrow(fly), 0L)
})

test_that("the species is part of every cache key", {
  skip_if_no_enrichment()
  gene_sets_or_skip("Rn")
  rat <- build_term_tables("hallmark", "Rn")
  human <- build_term_tables("hallmark", "Hs")
  expect_false(identical(rat$term2gene, human$term2gene))
  expect_true("Tp53" %in% rat$term2gene$gene)
  expect_false(identical(
    term_table_key("hallmark", "Rattus norvegicus", "s"),
    term_table_key("hallmark", "Homo sapiens", "s")))
  expect_false(identical(
    geneset_cache_file("d", "hallmark", "Rattus norvegicus"),
    geneset_cache_file("d", "hallmark", "Homo sapiens")))
})

test_that("symbols are matched ignoring case only when the exact match is poor", {
  ref <- c("Tp53", "Egfr", "Myc", paste0("Gene", 1:30))
  upper <- c("TP53", "EGFR", "MYC", toupper(paste0("Gene", 1:30)))
  plan <- symbol_case_plan(upper, ref)
  expect_true(plan$apply)
  expect_identical(plan$n_exact, 0L)
  expect_identical(apply_symbol_case(c("TP53", "unknown", NA), plan),
                   c("Tp53", "unknown", NA))
  # Already matching: nothing to do.
  expect_false(symbol_case_plan(ref, ref)$apply)
  # Too few genes gained to be worth second-guessing the data.
  expect_false(symbol_case_plan(c("TP53", "EGFR"), ref)$apply)
  # A target that differs from another set symbol only by case is never used.
  amb <- symbol_case_plan(upper, c(ref, "MYC"))
  expect_identical(apply_symbol_case("MYC", amb), "MYC")
  expect_identical(restore_symbol_case(c("Tp53/Egfr/X", NA),
                                       c(Tp53 = "TP53", Egfr = "EGFR")),
                   c("TP53/EGFR/X", NA))
})

test_that("human-cased genes against the mouse sets are matched ignoring case, with a note", {
  skip_if_no_enrichment()
  gene_sets_or_skip("Mm")
  b <- realistic_diff_bundle()
  res <- run_enrichment(b, type = "ora", database = "hallmark", organism = "Mm")
  expect_identical(res$params$symbol_case, "ignored")
  expect_match(res$warnings, "ignoring case", all = FALSE)
  df <- res$results$enrich_result_df
  expect_identical(df$pathway_name[which.min(df$adj_p_value)], "G2M CHECKPOINT")
  # The genes are named as the user's table names them.
  genes <- unlist(strsplit(df$overlap_features[which.min(df$adj_p_value)], "/"))
  expect_true(all(genes %in% b$results$diff_result_df$feature_symbol))

  same <- run_enrichment(b, type = "ora", database = "hallmark", organism = "Hs")
  expect_identical(same$params$symbol_case, "exact")
  expect_length(same$warnings, 0L)
})
