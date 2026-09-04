# Everything the deployment relies on has to survive `.omp`. A field
# that quietly does not round-trip is invisible until a user restores a
# session and finds a feature has forgotten what it knew -- the import
# fingerprint stops recognising a re-upload, or the exported script
# loses the path to the file it was generated from.

rt_input <- function(assay_type = "normalized_intensity", ...) {
  mat <- matrix(as.numeric(1:24), nrow = 6,
                dimnames = list(paste0("g", 1:6), paste0("s", 1:4)))
  meta <- data.frame(group = c("A", "A", "B", "B"),
                     row.names = paste0("s", 1:4))
  feat <- data.frame(feature_id = paste0("g", 1:6),
                     feature_symbol = paste0("SYM", 1:6))
  omics_input(mat, meta, feat, omics_type = "proteomics",
              assay_type = assay_type, ...)
}

save_and_load <- function(project) {
  path <- tempfile(fileext = ".omp")
  on.exit(unlink(path), add = TRUE)
  save_project(project, path)
  load_project(path)
}

test_that("the import fingerprint survives a save/load cycle", {
  skip_if_not_installed("qs2")
  inp <- rt_input(assay_type = "normalized_intensity",
                  source_fingerprint = "abc123:proteomics:intensity")
  back <- save_and_load(omics_project("p", experiments = list(proteomics = inp)))
  # Without this, a restored session cannot tell a re-upload of the same
  # file from new data, and would clear the user's analyses for nothing.
  expect_identical(back$experiments$proteomics$source_fingerprint,
                   "abc123:proteomics:intensity")
})

test_that("the archived file path survives a save/load cycle", {
  skip_if_not_installed("qs2")
  inp <- rt_input(source_path = "/srv/omicsapp/users/x/raw/cheek__ab12.xlsx")
  back <- save_and_load(omics_project("p", experiments = list(proteomics = inp)))
  expect_identical(back$experiments$proteomics$source_path,
                   "/srv/omicsapp/users/x/raw/cheek__ab12.xlsx")
})

test_that("a project saved before these fields existed still loads", {
  skip_if_not_installed("qs2")
  # Simulates an older .omp: the fields are simply absent. They are
  # optional and stay out of validate_omics_input()'s required list
  # precisely so this keeps working.
  inp <- rt_input()
  inp$source_fingerprint <- NULL
  inp$source_path <- NULL
  proj <- omics_project("legacy", experiments = list(proteomics = inp))
  back <- save_and_load(proj)
  expect_true(is_omics_project(back))
  expect_null(back$experiments$proteomics$source_fingerprint)
  expect_no_error(validate_omics_input(back$experiments$proteomics))
})

test_that("bundles survive well enough to export a script from a restored project", {
  skip_if_not_installed("qs2")
  skip_if_not_installed("limma")
  inp <- rt_input(assay_type = "normalized_intensity", source_path = "raw/cheek.xlsx")
  dif <- run_diff(inp, method = "limma", analysis_type = "group",
                  group_col = "group", control_group = "A", case_group = "B")
  proj <- omics_project("restored", experiments = list(proteomics = inp))
  proj$bundles <- list(diff = dif)

  back <- save_and_load(proj)
  lines <- export_script(back)

  # The point of restoring a session is being able to carry on, which
  # includes explaining what was already done.
  expect_true(any(grepl('"raw/cheek.xlsx"', lines, fixed = TRUE)))
  expect_true(any(grepl('method        = "limma"', lines, fixed = TRUE)))
  expect_false(any(grepl("# NOTE:", lines, fixed = TRUE)))
  expect_no_error(parse(text = paste(lines, collapse = "\n")))
})

test_that("a restored bundle still plots", {
  skip_if_not_installed("qs2")
  skip_if_not_installed("limma")
  inp <- rt_input(assay_type = "normalized_intensity")
  dif <- run_diff(inp, method = "limma", analysis_type = "group",
                  group_col = "group", control_group = "A", case_group = "B")
  proj <- omics_project("p", experiments = list(proteomics = inp))
  proj$bundles <- list(diff = dif)
  back <- save_and_load(proj)
  expect_s3_class(plot_volcano(back$bundles$diff), "ggplot")
})

test_that("save_project refuses to clobber unless told to", {
  skip_if_not_installed("qs2")
  path <- tempfile(fileext = ".omp"); on.exit(unlink(path), add = TRUE)
  proj <- omics_project("p", experiments = list(proteomics = rt_input()))
  save_project(proj, path)
  expect_error(save_project(proj, path), "already exists")
  expect_no_error(save_project(proj, path, overwrite = TRUE))
})

test_that("an interrupted save leaves the previous file intact", {
  skip_if_not_installed("qs2")
  path <- tempfile(fileext = ".omp"); on.exit(unlink(path), add = TRUE)
  good <- omics_project("good", experiments = list(proteomics = rt_input()))
  save_project(good, path)
  # save_project() writes to `path.tmp` and renames, so a failure part
  # way through cannot truncate what is already on disk.
  expect_error(save_project("not a project", path, overwrite = TRUE))
  expect_identical(load_project(path)$name, "good")
})
