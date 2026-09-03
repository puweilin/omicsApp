# run_enrichment()'s p_cutoff did two jobs: for ORA it chose which
# features went in, and for both backends it bounded the table that came
# out. One number could not serve both, so a bundle only ever held what
# the input threshold admitted and a caller could not offer raw or
# adjusted p as a reading of the same result.

test_that("both backends take an output bound separate from p_cutoff", {
  for (fn in list(run_ora_from_bundle, run_gsea_from_bundle)) {
    expect_true("output_p_cutoff" %in% names(formals(fn)))
  }
  expect_true("output_p_cutoff" %in% names(formals(run_enrichment)))
})

test_that("it defaults to p_cutoff, so old callers are unaffected", {
  # NULL, not 0.05: the default has to follow whatever p_cutoff was
  # given, not a fixed number that would quietly widen a caller who
  # passed a stricter one.
  expect_null(eval(formals(run_ora_from_bundle)$output_p_cutoff))
  expect_null(eval(formals(run_gsea_from_bundle)$output_p_cutoff))
  expect_null(eval(formals(run_enrichment)$output_p_cutoff))
})

test_that("filter_enrich_results reads whichever p it is told to", {
  # This is what the separation is for: one stored result, two readings.
  set.seed(5L)
  n <- 200L
  df <- data.frame(
    database = "go_bp", result_type = "gsea", comparison = "c",
    pathway_id = sprintf("GO%04d", seq_len(n)),
    pathway_name = sprintf("P%d", seq_len(n)),
    effect = rnorm(n), effect_type = "nes", direction = "both",
    p_value     = c(runif(60L, 0, 0.049), runif(n - 60L, 0.05, 1)),
    adj_p_value = c(runif(9L,  0, 0.049), runif(n - 9L,  0.05, 1)),
    q_value = NA_real_, gene_set_size = 100L, overlap_size = 20L,
    overlap_features = NA_character_, leading_features = NA_character_,
    source_label = "s", stringsAsFactors = FALSE)

  expect_equal(nrow(filter_enrich_results(df, 0.05, "adjusted")), 9L)
  expect_equal(nrow(filter_enrich_results(df, 0.05, "raw")), 60L)
  # And the stored table still holds everything, which is the point.
  expect_equal(nrow(df), 200L)
})
