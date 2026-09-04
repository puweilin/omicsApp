# A scaled PCA divides each column by its standard deviation, so one
# feature that never varies stops the decomposition:
#
#   cannot rescale a constant/zero column to unit variance
#
# In RNA-seq that is the ordinary case. A counts matrix carries every
# annotated gene and a fifth of them are zero in every sample of a given
# tissue -- 14,259 of 63,241 in the follicle dataset. The message names
# neither the features nor the reason, so it reads as a broken file.

mat_with_constants <- function(n_var = 20L, n_const = 5L, n_samples = 6L) {
  set.seed(42)
  varying <- matrix(stats::rnorm(n_var * n_samples), nrow = n_var)
  constant <- matrix(0, nrow = n_const, ncol = n_samples)
  out <- rbind(varying, constant)
  dimnames(out) <- list(paste0("F", seq_len(nrow(out))),
                        paste0("S", seq_len(n_samples)))
  out
}

test_that("a scaled prcomp is what fails on constant features", {
  # The counter-test: without it there is nothing to say this guard earns
  # its place.
  m <- mat_with_constants()
  expect_error(stats::prcomp(t(m), scale. = TRUE), "constant/zero")
})

test_that("pca_over_samples runs anyway and says how many it dropped", {
  m <- mat_with_constants(n_var = 20L, n_const = 5L)
  pca <- pca_over_samples(m)
  expect_s3_class(pca, "prcomp")
  expect_identical(attr(pca, "n_dropped"), 5L)
  expect_equal(nrow(pca$x), ncol(m))
})

test_that("nothing is dropped when every feature varies", {
  m <- mat_with_constants(n_var = 20L, n_const = 0L)
  expect_identical(attr(pca_over_samples(m), "n_dropped"), 0L)
})

test_that("a matrix with nothing left to decompose says so", {
  m <- mat_with_constants(n_var = 1L, n_const = 10L)
  expect_error(pca_over_samples(m), "vary across samples")
})

test_that("row_variance agrees with var() on which rows are constant", {
  m <- mat_with_constants()
  m[3L, 2L] <- NA_real_          # NAs must not turn a varying row constant
  fast <- unname(row_variance(m) > 0)
  slow <- vapply(seq_len(nrow(m)), function(i) {
    v <- stats::var(m[i, ], na.rm = TRUE)
    !is.na(v) && v > 0
  }, logical(1L))
  expect_identical(fast, slow)
})

test_that("an all-NA row counts as constant rather than erroring", {
  m <- mat_with_constants()
  m[1L, ] <- NA_real_
  expect_identical(unname(row_variance(m)[1L]), 0)
  expect_identical(attr(pca_over_samples(m), "n_dropped"), 6L)
})
