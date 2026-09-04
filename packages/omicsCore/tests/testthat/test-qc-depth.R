# The missingness panels answer the proteomics question: a peptide that
# was not detected is a hole in the matrix. A counts matrix has no holes
# -- every gene has a number for every sample and most are zero -- so
# the panel reported "63,241 features, all at 0%", which is true, carries
# no information, and occupied the space that should have been showing
# whether a library was under-sequenced.

depth_input <- function(lib_scale = rep(1, 6), n_feat = 40L) {
  ids <- paste0("S", seq_along(lib_scale))
  set.seed(7)
  base <- matrix(as.numeric(stats::rpois(n_feat * length(ids), 100)),
                 nrow = n_feat, dimnames = list(paste0("G", seq_len(n_feat)), ids))
  mat <- sweep(base, 2L, lib_scale, "*")
  omics_input(mat,
              data.frame(sample_id = ids, condition = "G1", row.names = ids,
                         stringsAsFactors = FALSE),
              data.frame(feature_id = rownames(mat), row.names = rownames(mat),
                         stringsAsFactors = FALSE),
              omics_type = "rnaseq", assay_type = "raw_count")
}

test_that("qc_depth reports library size and detection per sample", {
  d <- qc_depth(depth_input())
  expect_equal(nrow(d), 6L)
  expect_setequal(names(d), c("sample_id", "library_size", "n_detected",
                              "detection_rate", "library_size_ratio"))
  expect_true(all(d$library_size > 0))
  expect_true(all(d$detection_rate <= 1))
})

test_that("a shallow library is flagged against the median, not a constant", {
  # What counts as shallow depends entirely on the experiment; a fixed
  # count would be wrong for every study but one.
  d <- qc_depth(depth_input(c(1, 1, 1, 1, 1, 0.1)))
  expect_identical(qc_depth_outliers(d), "S6")

  # The same matrix scaled up is not suddenly healthy.
  d2 <- qc_depth(depth_input(c(1, 1, 1, 1, 1, 0.1) * 1000))
  expect_identical(qc_depth_outliers(d2), "S6")
})

test_that("an even cohort flags nothing", {
  expect_length(qc_depth_outliers(qc_depth(depth_input())), 0L)
})

test_that("detection counts signal, whether absence is written 0 or NA", {
  # A counts matrix says "nothing seen" with 0 and an intensity matrix
  # with NA. The question is the same one.
  inp <- depth_input()
  inp$expr_mat[1:10, "S1"] <- 0
  inp$expr_mat[1:10, "S2"] <- NA_real_
  d <- qc_depth(inp)
  n <- nrow(inp$expr_mat)
  expect_equal(d$n_detected[d$sample_id == "S1"], n - 10L)
  expect_equal(d$n_detected[d$sample_id == "S2"], n - 10L)
})

test_that("run_qc carries a depth summary for both modalities", {
  for (type in c("rnaseq", "proteomics")) {
    inp <- depth_input()
    inp$omics_type <- type
    inp$assay_type <- if (type == "rnaseq") "raw_count" else "raw_intensity"
    b <- run_qc(inp)
    expect_equal(nrow(b$results$qc_summary$depth), 6L, info = type)
  }
})

test_that("the depth view draws, and refuses a bundle that has none", {
  b <- run_qc(depth_input())
  expect_s3_class(plot_qc(b, view = "depth"), c("patchwork", "ggplot"))

  b$results$qc_summary$depth <- NULL
  expect_error(plot_qc(b, view = "depth"), "no depth summary")
})
