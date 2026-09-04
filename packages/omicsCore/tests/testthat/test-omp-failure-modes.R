# What load_project() says about a file that is not a project.
#
# The store hands these messages straight to the user ("Open failed:
# ..."), so each one must name the problem rather than the line of qs2
# that hit it. The future-schema case has a test already; these are the
# other ways a file arrives broken.

tiny_omp <- function() {
  inp <- realistic_input(n_per_group = 3L)
  p <- omics_project("p", experiments = list(proteomics = inp))
  f <- tempfile(fileext = ".omp")
  save_project(p, f)
  f
}

test_that("a truncated file fails with a message, not a crash", {
  skip_if_not_installed("qs2")
  f <- tiny_omp()
  bytes <- readBin(f, "raw", file.size(f))
  cut <- tempfile(fileext = ".omp")
  writeBin(bytes[seq_len(length(bytes) %/% 2L)], cut)
  err <- tryCatch(suppressWarnings(load_project(cut)),
                  error = function(e) conditionMessage(e))
  expect_type(err, "character")
  expect_gt(nchar(err), 10L)
})

test_that("a file that is not qs2 at all is refused", {
  skip_if_not_installed("qs2")
  txt <- tempfile(fileext = ".omp")
  writeLines("this is not a project", txt)
  expect_error(load_project(txt), "format|qs|read", ignore.case = TRUE)
  csv <- tempfile(fileext = ".omp")
  utils::write.csv(data.frame(a = 1:3), csv)
  expect_error(load_project(csv))
})

test_that("a qs2 file whose payload is not a project is refused by name", {
  skip_if_not_installed("qs2")
  f <- tempfile(fileext = ".omp")
  qs2::qs_save(list(schema_version = OMP_SCHEMA_VERSION, payload = list(a = 1)), f)
  expect_error(load_project(f), "omics_project")
  g <- tempfile(fileext = ".omp")
  qs2::qs_save(list(not_the = "shape"), g)
  expect_error(load_project(g))
})

test_that("an empty file and a missing file are told apart", {
  skip_if_not_installed("qs2")
  empty <- tempfile(fileext = ".omp")
  file.create(empty)
  expect_error(load_project(empty))
  expect_error(load_project(tempfile(fileext = ".omp")), "does not exist")
})

test_that("saving never leaves a partial file where a good one was", {
  skip_if_not_installed("qs2")
  f <- tiny_omp()
  before <- readBin(f, "raw", file.size(f))
  # A project that cannot be serialised: an environment with a cycle
  # is fine for qs2, so provoke the failure at the write instead by
  # making the target directory vanish mid-save.
  broken <- omics_project("broken", experiments = list(
    proteomics = realistic_input(n_per_group = 3L)))
  testthat::local_mocked_bindings(
    qs_save = function(object, file, ...) stop("disk full"),
    .package = "qs2"
  )
  expect_error(save_project(broken, f, overwrite = TRUE), "disk full")
  expect_identical(readBin(f, "raw", file.size(f)), before)
  expect_length(list.files(dirname(f), pattern = "\\.tmp$"), 0L)
})
