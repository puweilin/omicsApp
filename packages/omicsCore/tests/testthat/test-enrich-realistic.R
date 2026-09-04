# Enrichment on data where the answer is known.
#
# The fixture in helper-realistic.R shifts one Hallmark set (G2M
# checkpoint) up in the case group and leaves three others alone. Every
# test here asks the question a user would ask of the result -- did it
# find the pathway that was perturbed? -- rather than whether an object
# of the right class came back. The latter was all the suite asked, and
# it passed while GSEA was failing inside a tryCatch on every run.

collect_warnings <- function(expr) {
  msgs <- character(0)
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = msgs)
}

skip_if_no_enrichment <- function() {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("msigdbr")
}

g2m_label <- REAL_GENE_SET_LABELS[["G2M"]]

# ---- GSEA ---------------------------------------------------------------

test_that("GSEA finds the pathway that was perturbed, and does not fail quietly", {
  skip_if_no_enrichment()
  skip_if_not_installed("fgsea")
  b <- realistic_diff_bundle()
  set.seed(1)
  got <- collect_warnings(run_enrichment(b, type = "gsea", database = "hallmark"))

  # fgsea warns about its own estimates on a 160-gene list; those are
  # fine. What must not appear is the backend's own "GSEA failed", which
  # is what a fixture with made-up gene names produced.
  expect_false(any(grepl("GSEA failed", got$warnings, fixed = TRUE)),
               info = paste(got$warnings, collapse = " | "))

  df <- got$value$results$enrich_result_df
  expect_gt(nrow(df), 0L)
  top <- df$pathway_name[which.min(df$adj_p_value)]
  expect_identical(top, g2m_label)
  # NES is positive: the set went up in the case group
  expect_gt(df$effect[df$pathway_name == top], 0)
  expect_identical(df$direction[df$pathway_name == top], "up")
})

test_that("GSEA on counts, through DESeq2, finds the same pathway", {
  skip_if_no_enrichment()
  skip_if_not_installed("fgsea")
  skip_if_not_installed("DESeq2")
  b <- realistic_diff_bundle(method = "deseq2", omics_type = "rnaseq")
  set.seed(1)
  res <- suppressWarnings(run_enrichment(b, type = "gsea", database = "hallmark"))
  df <- res$results$enrich_result_df
  expect_gt(nrow(df), 0L)
  expect_identical(df$pathway_name[which.min(df$adj_p_value)], g2m_label)
})

test_that("GSEA is reproducible under a seed", {
  skip_if_no_enrichment()
  skip_if_not_installed("fgsea")
  b <- realistic_diff_bundle()
  set.seed(7)
  a <- suppressWarnings(run_enrichment(b, type = "gsea", database = "hallmark"))
  set.seed(7)
  z <- suppressWarnings(run_enrichment(b, type = "gsea", database = "hallmark"))
  expect_equal(a$results$enrich_result_df$effect, z$results$enrich_result_df$effect)
  expect_equal(a$results$enrich_result_df$p_value, z$results$enrich_result_df$p_value)
})

# ---- ORA ----------------------------------------------------------------

test_that("ORA finds the perturbed pathway", {
  skip_if_no_enrichment()
  b <- realistic_diff_bundle()
  res <- run_enrichment(b, type = "ora", database = "hallmark")
  df <- res$results$enrich_result_df
  expect_gt(nrow(df), 0L)
  expect_identical(df$pathway_name[which.min(df$adj_p_value)], g2m_label)
})

test_that("ORA enriches exactly the features the thresholds admit", {
  # The Enrichment view once enriched a different gene list from the one
  # in the hit table, silently. This pins the contract: what clusterProfiler
  # was handed is what filter_diff_results() returns for the same cutoffs.
  skip_if_no_enrichment()
  b <- realistic_diff_bundle()
  df <- b$results$diff_result_df
  for (pref in c("adjusted", "raw")) {
    res <- run_enrichment(b, type = "ora", database = "hallmark",
                          p_cutoff = 0.01, p_preference = pref,
                          effect_cutoff = 0.5)
    handed <- res$results$enrich_object[["both__hallmark"]]@gene
    expected <- unique(stats::na.omit(filter_diff_results(
      df, p_cutoff = 0.01, p_preference = pref, effect_cutoff = 0.5
    )$feature_symbol))
    expect_setequal(handed, expected)
  }
})

test_that("ORA by direction splits the list by the sign of the effect", {
  skip_if_no_enrichment()
  b <- realistic_diff_bundle()
  up <- run_enrichment(b, type = "ora", database = "hallmark", direction = "up")
  down <- run_enrichment(b, type = "ora", database = "hallmark", direction = "down")
  expect_identical(up$results$enrich_result_df$pathway_name[
    which.min(up$results$enrich_result_df$adj_p_value)], g2m_label)
  # Nothing was pushed down, so the down list carries no G2M signal
  expect_false(g2m_label %in% down$results$enrich_result_df$pathway_name)
})

test_that("a null dataset produces no confident pathway", {
  skip_if_no_enrichment()
  b <- realistic_diff_bundle(signal = NULL)
  res <- run_enrichment(b, type = "ora", database = "hallmark")
  df <- res$results$enrich_result_df
  expect_true(nrow(df) == 0L || all(df$adj_p_value > 0.05))
})

# ---- GSVA ---------------------------------------------------------------

test_that("GSVA scores the perturbed set higher in the perturbed group", {
  skip_if_not_installed("GSVA")
  skip_if_not_installed("msigdbr")
  inp <- realistic_input()
  g <- run_gsva(inp, database = "hallmark")
  m <- g$results$gsva_matrix
  expect_true(is.matrix(m))
  expect_identical(colnames(m), colnames(inp$expr_mat))

  row <- grep("G2M_CHECKPOINT", rownames(m), value = TRUE)
  expect_length(row, 1L)
  grp <- inp$meta_df[colnames(m), "group"]
  shifted <- stats::t.test(m[row, grp == "G2"], m[row, grp == "G1"],
                           alternative = "greater")
  expect_lt(shifted$p.value, 1e-3)

  # GSVA scores are relative within a sample, so pushing 40 genes up
  # pushes every other set's ranks down a little. The untouched sets
  # therefore do move -- the other way, and by less. What identifies the
  # perturbed set is that it moved up more than anything else.
  shift_by_set <- rowMeans(m[, grp == "G2", drop = FALSE]) -
    rowMeans(m[, grp == "G1", drop = FALSE])
  expect_identical(names(which.max(shift_by_set)), row)
  quiet <- grep("OXIDATIVE_PHOSPHORYLATION", rownames(m), value = TRUE)
  expect_length(quiet, 1L)
  expect_lt(shift_by_set[[quiet]], 0)
})

test_that("GSVA accepts the caller's own gene sets", {
  skip_if_not_installed("GSVA")
  g <- run_gsva(realistic_input(), gene_sets = REAL_GENE_SETS, min_size = 5L)
  m <- g$results$gsva_matrix
  expect_setequal(rownames(m), names(REAL_GENE_SETS))
  expect_identical(g$params$database, NA_character_)
})

test_that("GSVA on counts log-transforms before scoring", {
  skip_if_not_installed("GSVA")
  inp <- realistic_input("rnaseq")
  g <- run_gsva(inp, gene_sets = REAL_GENE_SETS, min_size = 5L)
  m <- g$results$gsva_matrix
  grp <- inp$meta_df[colnames(m), "group"]
  expect_lt(stats::t.test(m["G2M", grp == "G2"], m["G2M", grp == "G1"],
                          alternative = "greater")$p.value, 1e-3)
})
