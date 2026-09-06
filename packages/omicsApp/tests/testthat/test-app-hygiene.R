# What the demo leaves behind in a shared process.
#
# The example data seeded the global random stream with set.seed(),
# and every session in a Shiny process shares that stream: opening the
# demo in one tab reset it for an imputation running in another. The
# seed is now local to the function that asked for it, and the stream
# is put back on the way out.

seed_state <- function() {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  }
}

demo_calls <- list(
  proteomics = function() example_input("proteomics"),
  rnaseq = function() example_input("rnaseq"),
  project = function() example_project(),
  diff = function() example_diff_bundle(),
  enrich = function() example_enrich_bundle(),
  qc = function() example_qc_bundle(),
  integration = function() example_integration_bundle(),
  template_proteomics = function() import_template_sheets("proteomics"),
  template_rnaseq = function() import_template_sheets("rnaseq")
)

test_that("the demo leaves the caller's random stream where it was", {
  withr::defer(set.seed(1))
  for (nm in names(demo_calls)) {
    set.seed(20260906L)
    before <- seed_state()
    demo_calls[[nm]]()
    expect_identical(seed_state(), before, label = nm)
  }
})

test_that("the demo plants no seed in a session that had none", {
  withr::defer(set.seed(1))
  for (nm in c("proteomics", "rnaseq", "template_rnaseq")) {
    suppressWarnings(rm(".Random.seed", envir = globalenv()))
    demo_calls[[nm]]()
    expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE), label = nm)
  }
})

test_that("the demo is the same data whatever the stream was", {
  withr::defer(set.seed(1))
  set.seed(1); a <- example_input("rnaseq")
  set.seed(2); b <- example_input("rnaseq")
  expect_identical(a$expr_mat, b$expr_mat)
  set.seed(3); s1 <- import_template_sheets("proteomics")
  set.seed(4); s2 <- import_template_sheets("proteomics")
  expect_identical(s1, s2)
})

test_that("a frame that seeds twice restores to its first state", {
  withr::defer(set.seed(1))
  twice <- function() {
    local_seed(1L)
    stats::runif(1)
    local_seed(2L)
    stats::runif(1)
  }
  environment(twice) <- asNamespace("omicsApp")
  set.seed(20260906L)
  before <- seed_state()
  twice()
  expect_identical(seed_state(), before)
})

test_that("package code seeds the stream only through local_seed()", {
  r_dir <- file.path("..", "..", "R")
  skip_if_not(dir.exists(r_dir), "package source is not beside the tests")
  files <- list.files(r_dir, pattern = "\\.R$", full.names = TRUE)
  seeders <- Filter(function(f) {
    code <- grep("^\\s*#", readLines(f, warn = FALSE), value = TRUE, invert = TRUE)
    any(grepl("set\\.seed\\(", code))
  }, files)
  expect_identical(basename(seeders), "example_data.R")
})
