# run_enrichment_grid(): every combination in one pass.
#
# Kept to one database. What is worth asserting here is the shape of the
# grid and what it does when a combination fails -- not the enrichment
# itself, which test-enrich.R covers, and running seven databases three
# ways would put a minute into every test run to say nothing new.

grid_diff_bundle <- function(n = 800L, ns = 8L) {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("msigdbr")
  set.seed(4)
  sets <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
  sym <- make.unique(rep_len(unique(sets$gene_symbol), n))
  m <- matrix(stats::rnorm(n * ns, 8, 1), nrow = n,
              dimnames = list(paste0("g", seq_len(n)), paste0("s", seq_len(ns))))
  m[1:120, (ns / 2 + 1):ns] <- m[1:120, (ns / 2 + 1):ns] + 2
  x <- omics_input(m,
                   data.frame(group = rep(c("a", "b"), each = ns / 2),
                              row.names = colnames(m)),
                   data.frame(feature_id = rownames(m), feature_symbol = sym,
                              row.names = rownames(m)),
                   omics_type = "proteomics",
                   assay_type = "normalized_intensity")
  run_diff(x, method = "ttest", analysis_type = "group", group_col = "group",
           control_group = "a", case_group = "b")
}

# GSEA's `direction` filters its output; it does not change the analysis,
# which reads every feature in rank order. Running it once per direction
# would be the same computation three times, and the rows say which way
# each pathway went anyway. ORA is the opposite: direction decides which
# genes are tested, so all three are real.
test_that("ORA runs per direction and GSEA runs once per database", {
  b <- grid_diff_bundle()
  g <- run_enrichment_grid(b, databases = "hallmark")
  lg <- g$results$grid_log

  expect_equal(sum(lg$type == "ora"), 3L)
  expect_setequal(lg$direction[lg$type == "ora"], c("both", "up", "down"))
  expect_equal(sum(lg$type == "gsea"), 1L)
  expect_identical(lg$direction[lg$type == "gsea"], "both")
  expect_equal(g$params$n_combinations, nrow(lg))
})

test_that("the stacked rows say which combination they came from", {
  b <- grid_diff_bundle()
  g <- run_enrichment_grid(b, databases = "hallmark")
  df <- g$results$enrich_result_df
  skip_if(nrow(df) == 0L, "no pathways reached significance in this fixture")

  expect_true(all(c("database", "direction", "result_type") %in% names(df)))
  expect_true(all(df$database == "hallmark"))
  expect_true(all(df$result_type %in% c("ora", "gsea")))
  # Stacking is only legitimate because the schema is shared.
  expect_true(all(ENRICH_RESULT_REQUIRED_COLS %in% names(df)))
})

# Seven databases against one differential will not all succeed: an
# organism a database does not cover, gene sets that cannot be fetched, a
# direction with nothing in it. One of those must not cost the other
# twenty-seven.
test_that("a combination that fails does not stop the grid", {
  b <- grid_diff_bundle()
  calls <- 0L
  testthat::local_mocked_bindings(
    run_enrichment = function(diff_bundle, type, database, direction, ...) {
      calls <<- calls + 1L
      if (identical(direction, "up")) stop("gene sets unavailable")
      new_analysis_bundle(
        analysis_name = "run_enrichment",
        results = list(enrich_result_df = empty_enrich_result_df()))
    },
    .package = "omicsCore"
  )
  g <- run_enrichment_grid(b, databases = "hallmark")

  expect_equal(calls, 4L)                       # all four still attempted
  lg <- g$results$grid_log
  expect_equal(sum(lg$status == "failed"), 1L)
  expect_match(lg$message[lg$status == "failed"], "gene sets unavailable")
  # and the failure is carried where a reader will see it
  expect_length(g$warnings, 1L)
  expect_match(g$warnings, "hallmark")
})

test_that("the log accounts for every combination, empty included", {
  b <- grid_diff_bundle()
  g <- run_enrichment_grid(b, databases = "hallmark")
  lg <- g$results$grid_log
  expect_true(all(lg$status %in% c("ok", "empty", "failed")))
  expect_equal(sum(lg$rows), nrow(g$results$enrich_result_df))
  expect_true(all(lg$seconds >= 0))
})

test_that("run_enrichment_grid rejects a bundle that is not a differential", {
  expect_error(run_enrichment_grid(list()), "run_diff")
  # A real bundle of the wrong kind, which is the case that would
  # otherwise reach run_enrichment() and fail somewhere less obvious.
  not_diff <- new_analysis_bundle(analysis_name = "run_qc",
                                  results = list(qc_summary = list()))
  expect_error(run_enrichment_grid(not_diff), "run_diff")
})

test_that("on_progress is told where the grid has got to", {
  b <- grid_diff_bundle()
  seen <- list()
  run_enrichment_grid(
    b, databases = "hallmark",
    on_progress = function(done, total, label) {
      seen[[length(seen) + 1L]] <<- list(done = done, total = total,
                                         label = label)
    })
  expect_gt(length(seen), 1L)
  expect_true(all(vapply(seen, function(s) s$total, numeric(1)) == 4L))
  # ends at the end, so a caller driving a progress bar can close it
  expect_equal(seen[[length(seen)]]$done, 4L)
})
