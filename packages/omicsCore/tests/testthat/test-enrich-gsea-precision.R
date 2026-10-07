# fgsea's precision settings, `eps` and `nPermSimple`, reach GSEA through
# run_enrichment(eps =, n_perm_simple =). Left alone they must give the
# result GSEA gave before they could be set; set, they must take effect,
# be recorded, and be repeated by the exported script.

skip_if_no_gsea <- function() {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("msigdbr")
  skip_if_not_installed("fgsea")
}

# G2M shifted far enough that its true p-value is well below fgsea's
# default floor of 1e-10.
strong_gsea_bundle <- function() {
  run_diff(realistic_input("proteomics", effect = 4), method = "limma",
           analysis_type = "group", group_col = "group",
           control_group = "G1", case_group = "G2")
}

quiet_gsea <- function(...) {
  suppressMessages(suppressWarnings(run_enrichment(..., type = "gsea",
                                                    database = "hallmark",
                                                    output_p_cutoff = 1)))
}

test_that("GSEA at the default precision is the call made before it could be set", {
  skip_if_no_gsea()
  b <- strong_gsea_bundle()
  got <- quiet_gsea(b)$results$enrich_result_df

  # The call as it was written before eps / nPermSimple were passed.
  df <- b$results$diff_result_df
  ranked <- gsea_rank_vector(df, "feature_symbol")
  ranked <- sort(ranked[!duplicated(names(ranked))], decreasing = TRUE)
  terms <- build_term_tables(database = "hallmark", organism = "Hs")
  old <- with_fixed_seed(123L, suppressMessages(suppressWarnings(
    clusterProfiler::GSEA(geneList = ranked, TERM2GENE = terms$term2gene,
                          TERM2NAME = terms$term2name, pvalueCutoff = 1,
                          pAdjustMethod = "BH", minGSSize = 10L,
                          maxGSSize = 500L, seed = TRUE))))
  old <- as.data.frame(old)
  expect_setequal(got$pathway_id, old$ID)
  expect_equal(got$p_value[match(old$ID, got$pathway_id)], old$pvalue)
  expect_equal(got$effect[match(old$ID, got$pathway_id)], old$NES)
  # The floor at work: the perturbed pathway is reported at 1e-10.
  expect_equal(min(got$p_value), 1e-10)
})

test_that("eps = 0 estimates p-values below the default floor", {
  skip_if_no_gsea()
  b <- strong_gsea_bundle()
  g2m <- "HALLMARK_G2M_CHECKPOINT"
  floor <- quiet_gsea(b)$results$enrich_result_df
  free <- quiet_gsea(b, eps = 0)$results$enrich_result_df
  expect_equal(floor$p_value[floor$pathway_id == g2m], 1e-10)
  expect_lt(free$p_value[free$pathway_id == g2m], 1e-20)
  # Pathways nowhere near the floor are untouched by it.
  weak <- floor$pathway_id[floor$p_value > 1e-3]
  expect_equal(free$p_value[match(weak, free$pathway_id)],
               floor$p_value[match(weak, floor$pathway_id)])
})

test_that("eps and n_perm_simple are recorded for GSEA, not for ORA, and checked", {
  skip_if_no_gsea()
  b <- strong_gsea_bundle()
  def <- quiet_gsea(b)
  expect_identical(def$params$eps, 1e-10)
  expect_identical(def$params$n_perm_simple, 1000L)
  set <- quiet_gsea(b, eps = 0, n_perm_simple = 5000)
  expect_identical(set$params$eps, 0)
  expect_identical(set$params$n_perm_simple, 5000L)
  # A larger first stage changes the estimates of the middling pathways.
  expect_false(isTRUE(all.equal(
    def$results$enrich_result_df$p_value,
    set$results$enrich_result_df$p_value[match(def$results$enrich_result_df$pathway_id,
                                               set$results$enrich_result_df$pathway_id)])))

  ora <- suppressMessages(suppressWarnings(
    run_enrichment(b, type = "ora", database = "hallmark")))
  expect_null(ora$params$eps)
  expect_null(ora$params$n_perm_simple)

  expect_error(run_enrichment(b, type = "gsea", eps = -1), "`eps`")
  expect_error(run_enrichment(b, type = "gsea", n_perm_simple = 10.5), "`n_perm_simple`")
})

test_that("the exported script repeats eps and n_perm_simple, and reproduces the run", {
  skip_if_no_gsea()
  inp <- realistic_input("proteomics", effect = 4)
  b <- run_diff(inp, method = "limma", analysis_type = "group",
                group_col = "group", control_group = "G1", case_group = "G2")
  enr <- quiet_gsea(b, eps = 0, n_perm_simple = 2000L)
  proj <- omics_project("gsea", experiments = list(prot = inp))
  proj$bundles <- list(diff = b, enrich = enr)
  lines <- export_script(proj, include_plots = FALSE)
  call <- paste(lines, collapse = "\n")
  expect_match(call, "eps\\s+= 0\\b")
  expect_match(call, "n_perm_simple\\s+= 2000L")

  # Run the enrichment call as written against the same differential result.
  start <- grep("^enrich <- run_enrichment\\(", lines)
  end <- start + which(lines[start:length(lines)] == ")")[[1L]] - 1L
  env <- new.env(parent = asNamespace("omicsCore"))
  env$diff <- b
  again <- suppressMessages(suppressWarnings(
    eval(parse(text = lines[start:end]), envir = env)))
  expect_equal(env$enrich$results$enrich_result_df$p_value,
               enr$results$enrich_result_df$p_value)

  # A GSEA result saved before the settings existed carries neither, so
  # its script calls GSEA at the defaults -- the run that made it.
  old <- quiet_gsea(b)
  old$params$eps <- NULL
  old$params$n_perm_simple <- NULL
  proj$bundles$enrich <- old
  old_lines <- export_script(proj, include_plots = FALSE)
  expect_false(any(grepl("eps|n_perm_simple", old_lines)))
})
