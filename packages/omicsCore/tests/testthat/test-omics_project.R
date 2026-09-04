make_proteo <- function() {
  set.seed(1)
  expr <- matrix(
    rnorm(15), nrow = 3,
    dimnames = list(c("P1", "P2", "P3"), c("S1", "S2", "S3", "S4", "S5"))
  )
  meta <- data.frame(
    group = c("A", "A", "B", "B", "B"),
    row.names = colnames(expr),
    stringsAsFactors = FALSE
  )
  feat <- data.frame(
    feature_id = rownames(expr),
    row.names = rownames(expr),
    stringsAsFactors = FALSE
  )
  omics_input(expr, meta, feat, omics_type = "proteomics",
              assay_type = "normalized_intensity")
}

make_rna <- function() {
  set.seed(2)
  expr <- matrix(
    rpois(20, lambda = 100), nrow = 4,
    dimnames = list(c("ENSG1", "ENSG2", "ENSG3", "ENSG4"),
                    c("R1", "R2", "R3", "R4", "R5"))
  )
  meta <- data.frame(
    group = c("A", "A", "B", "B", "B"),
    row.names = colnames(expr),
    stringsAsFactors = FALSE
  )
  feat <- data.frame(
    feature_id = rownames(expr),
    row.names = rownames(expr),
    stringsAsFactors = FALSE
  )
  omics_input(expr, meta, feat, omics_type = "rnaseq",
              assay_type = "raw_count")
}

test_that("omics_project() constructs an empty project", {
  p <- omics_project("demo")
  expect_true(is_omics_project(p))
  expect_equal(p$name, "demo")
  expect_length(p$experiments, 0)
  expect_null(p$sample_link)
})

test_that("omics_project() accepts pre-populated experiments", {
  p <- omics_project(
    "demo",
    experiments = list(proteomics = make_proteo(), rnaseq = make_rna())
  )
  expect_equal(experiment_tags(p), c("proteomics", "rnaseq"))
})

test_that("omics_project() rejects unnamed or non-input entries", {
  expect_error(
    omics_project("d", experiments = list(make_proteo())),
    "must be named"
  )
  expect_error(
    omics_project("d", experiments = list(proteomics = list(a = 1))),
    "omics_input"
  )
})

test_that("add_experiment / remove_experiment round-trip", {
  p <- omics_project("demo")
  p <- add_experiment(p, "proteomics", make_proteo())
  p <- add_experiment(p, "rnaseq", make_rna())
  expect_equal(experiment_tags(p), c("proteomics", "rnaseq"))

  expect_error(add_experiment(p, "proteomics", make_proteo()), "already exists")

  p <- remove_experiment(p, "rnaseq")
  expect_equal(experiment_tags(p), "proteomics")
  expect_error(remove_experiment(p, "rnaseq"), "not found")
})

test_that("sample_link is validated against experiment tags", {
  proteo <- make_proteo()
  link <- data.frame(
    tag = c("proteomics", "proteomics", "rnaseq"),
    sample_id = c("S1", "S2", "R1"),
    donor_id = c("D1", "D2", "D1"),
    stringsAsFactors = FALSE
  )
  expect_error(
    omics_project("d",
                  experiments = list(proteomics = proteo),
                  sample_link = link),
    "rnaseq"
  )
})

test_that("sample_link missing required cols errors", {
  proteo <- make_proteo()
  link <- data.frame(
    tag = "proteomics", sample_id = "S1",
    stringsAsFactors = FALSE
  )
  expect_error(
    omics_project("d",
                  experiments = list(proteomics = proteo),
                  sample_link = link),
    "donor_id"
  )
})
