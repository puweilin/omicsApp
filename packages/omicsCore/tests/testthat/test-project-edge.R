# Edge-case tests for omics_project, add_experiment, remove_experiment,
# experiment_tags, load_project, save_project.

make_proj_input <- function() {
  mat <- matrix(rnorm(40, mean = 10, sd = 2), nrow = 10, ncol = 4)
  rownames(mat) <- paste0("gene_", 1:10)
  colnames(mat) <- paste0("s", 1:4)
  meta <- data.frame(group = c("A", "A", "B", "B"),
                     row.names = paste0("s", 1:4))
  feat <- data.frame(feature_id = paste0("gene_", 1:10),
                     stringsAsFactors = FALSE)
  omics_input(mat, meta, feat, omics_type = "proteomics",
              assay_type = "normalized_intensity")
}

# ---- omics_project: construction ---------------------------------------

test_that("omics_project creates empty project", {
  p <- omics_project("test_proj")
  expect_s3_class(p, "omics_project")
  expect_equal(p$name, "test_proj")
  expect_equal(length(p$experiments), 0)
  expect_equal(experiment_tags(p), character(0))
})

test_that("omics_project errors on non-string name", {
  expect_error(omics_project(123), "single non-NA string")
  expect_error(omics_project(NA_character_), "single non-NA string")
  expect_error(omics_project(c("a", "b")), "single non-NA string")
})

test_that("omics_project errors when experiments are not named", {
  inp <- make_proj_input()
  expect_error(omics_project("p", experiments = list(inp)), "must be named")
})

test_that("omics_project errors when experiments have duplicate names", {
  inp <- make_proj_input()
  expect_error(
    omics_project("p", experiments = list(prot = inp, prot = inp)),
    "unique"
  )
})

test_that("omics_project errors when experiments contain non-omics_input", {
  inp <- make_proj_input()
  expect_error(
    omics_project("p", experiments = list(valid = inp, bad = list())),
    "omics_input"
  )
})

test_that("omics_project accepts valid experiments", {
  inp <- make_proj_input()
  p <- omics_project("multi", experiments = list(prot = inp))
  expect_equal(experiment_tags(p), "prot")
  expect_true(is_omics_project(p))
})

test_that("omics_project stores sample_link when valid", {
  inp <- make_proj_input()
  sl <- data.frame(tag = "prot", sample_id = "s1", donor_id = "D1",
                   stringsAsFactors = FALSE)
  p <- omics_project("p", experiments = list(prot = inp), sample_link = sl)
  expect_equal(nrow(p$sample_link), 1)
})

test_that("omics_project errors when sample_link references unknown tags", {
  sl <- data.frame(tag = "nonexistent", sample_id = "x", donor_id = "D1",
                   stringsAsFactors = FALSE)
  inp <- make_proj_input()
  expect_error(
    omics_project("p", experiments = list(prot = inp), sample_link = sl),
    "sample_link"
  )
})

test_that("omics_project errors when feature_link is not a data.frame", {
  expect_error(
    omics_project("p", feature_link = "not a df"),
    "data.frame"
  )
})

test_that("omics_project stores metadata", {
  p <- omics_project("p", metadata = list(author = "test", version = 1L))
  expect_equal(p$metadata$author, "test")
  expect_equal(p$metadata$version, 1L)
})

test_that("omics_project has correct structure fields", {
  p <- omics_project("test")
  expect_true("created_at" %in% names(p))
  expect_equal(p$schema_version, "1")
})

# ---- is_omics_project --------------------------------------------------

test_that("is_omics_project returns TRUE for omics_project objects", {
  p <- omics_project("test")
  expect_true(is_omics_project(p))
})

test_that("is_omics_project returns FALSE for other objects", {
  expect_false(is_omics_project(list()))
  expect_false(is_omics_project(NULL))
  expect_false(is_omics_project(data.frame()))
  expect_false(is_omics_project("string"))
})

# ---- add_experiment ----------------------------------------------------

test_that("add_experiment adds a new experiment", {
  p <- omics_project("test")
  inp <- make_proj_input()
  p2 <- add_experiment(p, "proteomics", inp)
  expect_equal(experiment_tags(p2), "proteomics")
  expect_equal(p2$experiments$proteomics, inp)
})

test_that("add_experiment errors on duplicate tag", {
  p <- omics_project("test")
  inp <- make_proj_input()
  p <- add_experiment(p, "prot", inp)
  expect_error(add_experiment(p, "prot", inp), "already exists")
})

test_that("add_experiment errors on invalid project", {
  inp <- make_proj_input()
  expect_error(add_experiment(list(), "prot", inp), "omics_project")
})

test_that("add_experiment errors on empty tag name", {
  p <- omics_project("test")
  inp <- make_proj_input()
  expect_error(add_experiment(p, "", inp), "non-empty string")
  expect_error(add_experiment(p, NA_character_, inp), "non-empty string")
})

test_that("add_experiment errors on non-omics_input", {
  p <- omics_project("test")
  expect_error(add_experiment(p, "prot", list()), "omics_input")
})

test_that("add_experiment validates the input before adding", {
  p <- omics_project("test")
  inp <- make_proj_input()
  inp$expr_mat <- NULL
  expect_error(add_experiment(p, "prot", inp), "expr_mat")
})

# ---- remove_experiment -------------------------------------------------

test_that("remove_experiment removes an existing experiment", {
  p <- omics_project("test")
  inp <- make_proj_input()
  p <- add_experiment(p, "prot", inp)
  p <- add_experiment(p, "rna", inp)
  expect_equal(length(experiment_tags(p)), 2)
  p2 <- remove_experiment(p, "prot")
  expect_equal(experiment_tags(p2), "rna")
  expect_null(p2$experiments$prot)
})

test_that("remove_experiment errors when tag not found", {
  p <- omics_project("test")
  expect_error(remove_experiment(p, "nonexistent"), "not found")
})

test_that("remove_experiment errors on invalid project", {
  expect_error(remove_experiment(list(), "x"), "omics_project")
})

# ---- experiment_tags ---------------------------------------------------

test_that("experiment_tags returns character(0) for empty project", {
  p <- omics_project("test")
  expect_equal(experiment_tags(p), character(0))
})

test_that("experiment_tags errors on invalid project", {
  expect_error(experiment_tags(list()), "omics_project")
})

# ---- print method ------------------------------------------------------

test_that("print.omics_project does not error", {
  p <- omics_project("test")
  expect_output(print(p), "omics_project")
  inp <- make_proj_input()
  p <- add_experiment(p, "prot", inp)
  expect_output(print(p), "prot")
  expect_output(print(p), "proteomics")
  expect_invisible(print(p))
})
