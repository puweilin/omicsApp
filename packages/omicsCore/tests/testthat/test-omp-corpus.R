# Every .omp file ever written must keep opening.
#
# A project file is the only durable record of a user's work, and the
# code that reads it changes every week. test-persistence-round-trip.R
# proves a file written by *this* version reads back; it cannot prove
# that a file written by last month's version still does, because it
# writes its fixture fresh each run. This corpus does not: every file
# under fixtures/omp/ was written by some past version and is committed
# as bytes. When the schema or the bundle layout changes, add a file
# rather than replacing one -- the point is that the old ones stay.
#
# To add the current version's file:
#
#   devtools::load_all("packages/omicsCore")
#   testthat::test_file("packages/omicsCore/tests/testthat/test-omp-corpus.R")
#   write_omp_corpus_fixture()      # defined below; sourced by the test
#
# and commit the new file under tests/testthat/fixtures/omp/.

corpus_dir <- testthat::test_path("fixtures", "omp")

corpus_file_name <- function() {
  sprintf("schema-%s__omicsCore-%s.omp",
          OMP_SCHEMA_VERSION,
          as.character(utils::packageVersion("omicsCore")))
}

write_omp_corpus_fixture <- function(dir = corpus_dir) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  inp <- realistic_input(n_per_group = 3L)
  diff <- run_diff(inp, method = "ttest", analysis_type = "group",
                   group_col = "group", control_group = "G1",
                   case_group = "G2")
  project <- omics_project(
    name = "corpus",
    experiments = list(proteomics = inp)
  )
  project$bundles <- list(diff = diff)
  path <- file.path(dir, corpus_file_name())
  save_project(project, path, overwrite = TRUE)
  path
}

corpus_files <- function() {
  list.files(corpus_dir, pattern = "\\.omp$", full.names = TRUE)
}

test_that("the corpus holds a file written under the current schema version", {
  skip_if_not_installed("qs2")
  expected <- sprintf("schema-%s__", OMP_SCHEMA_VERSION)
  present <- basename(corpus_files())
  expect_true(
    any(startsWith(present, expected)),
    info = sprintf(
      "No corpus file for schema %s. Run write_omp_corpus_fixture() and commit %s.",
      OMP_SCHEMA_VERSION, corpus_file_name()
    )
  )
})

for (f in corpus_files()) {
  test_that(paste("still opens:", basename(f)), {
    skip_if_not_installed("qs2")
    project <- load_project(f)
    expect_true(is_omics_project(project))
    expect_gt(length(project$experiments), 0L)
    for (layer in project$experiments) {
      expect_true(is_omics_input(layer))
      expect_silent(suppressWarnings(validate_omics_input(layer)))
      expect_true(is.matrix(layer$expr_mat))
      expect_identical(rownames(layer$meta_df), colnames(layer$expr_mat))
    }
  })

  test_that(paste("its analyses are still usable:", basename(f)), {
    skip_if_not_installed("qs2")
    project <- load_project(f)
    for (b in project$bundles) {
      expect_true(is_analysis_bundle(b))
    }
    diff <- project$bundles$diff
    skip_if(is.null(diff), "no diff bundle in this file")
    expect_s3_class(plot_volcano(diff), "ggplot")
    kept <- filter_diff_results(diff$results$diff_result_df, p_cutoff = 1)
    expect_identical(nrow(kept), nrow(diff$results$diff_result_df))
    script <- export_script(project)
    expect_no_error(parse(text = script))
  })
}
