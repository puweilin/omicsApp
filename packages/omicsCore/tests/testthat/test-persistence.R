# Round-trip + gating tests for save_project() / load_project().

make_persistence_project <- function() {
  set.seed(2027)
  expr <- matrix(
    rnorm(30), nrow = 5,
    dimnames = list(paste0("F", 1:5), paste0("S", 1:6))
  )
  meta <- data.frame(group = rep(c("A", "B"), each = 3),
                     row.names = paste0("S", 1:6),
                     stringsAsFactors = FALSE)
  feat <- data.frame(feature_id = paste0("F", 1:5),
                     row.names = paste0("F", 1:5),
                     stringsAsFactors = FALSE)
  input <- omics_input(expr, meta, feat, omics_type = "proteomics")
  omics_project("persist_demo", experiments = list(proteo = input))
}

test_that("save_project requires an omics_project", {
  tmp <- tempfile(fileext = ".omp")
  expect_error(save_project(list(), tmp), "omics_project")
})

test_that("save_project / load_project round-trips a project", {
  skip_if_not_installed("qs2")
  p <- make_persistence_project()
  tmp <- tempfile(fileext = ".omp")
  on.exit(unlink(tmp), add = TRUE)

  out <- save_project(p, tmp)
  expect_identical(out, tmp)
  expect_true(file.exists(tmp))

  q <- load_project(tmp)
  expect_true(is_omics_project(q))
  expect_equal(q$name, p$name)
  expect_equal(experiment_tags(q), experiment_tags(p))
  expect_equal(q$experiments$proteo$expr_mat,
               p$experiments$proteo$expr_mat)
})

test_that("save_project refuses to overwrite by default", {
  skip_if_not_installed("qs2")
  p <- make_persistence_project()
  tmp <- tempfile(fileext = ".omp")
  on.exit(unlink(tmp), add = TRUE)
  save_project(p, tmp)
  expect_error(save_project(p, tmp), "already exists")
  expect_silent(save_project(p, tmp, overwrite = TRUE))
})

test_that("load_project rejects non-omp files", {
  skip_if_not_installed("qs2")
  tmp <- tempfile(fileext = ".omp")
  on.exit(unlink(tmp), add = TRUE)
  qs2 <- asNamespace("qs2")
  qs2$qs_save(list(foo = 1), file = tmp)
  expect_error(load_project(tmp), "omicsCore project file")
})

test_that("load_project rejects future schema versions", {
  skip_if_not_installed("qs2")
  p <- make_persistence_project()
  tmp <- tempfile(fileext = ".omp")
  on.exit(unlink(tmp), add = TRUE)
  qs2 <- asNamespace("qs2")
  qs2$qs_save(list(schema_version = "999.0.0", payload = p), file = tmp)
  expect_error(load_project(tmp), "Unsupported project schema version")
})

test_that("load_project errors when file is missing", {
  expect_error(load_project(tempfile()), "does not exist")
})

test_that("save_project errors informatively without qs2", {
  if (requireNamespace("qs2", quietly = TRUE)) {
    skip("qs2 is installed; cannot test gate")
  }
  expect_error(save_project(make_persistence_project(), tempfile()),
               "qs2")
})
