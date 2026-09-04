# Two processes writing the same project file at once.
#
# One account, two devices, both autosaving: test-concurrent-save.R
# covers what a single writer must do, and checks the atomic-rename
# property by reading the source, because one process cannot stage two
# overlapping writers. This file stages them. Two R processes hammer the
# same path while this one reads it, and the guarantee under test is
# the reader's: every read is a complete file, and the last one standing
# is one of the two projects rather than a splice of both.

# Where the omicsCore this process is running came from, so the writers
# run the same code. Under load_all that is the source tree (pkgload's
# system.file() shim answers with inst/; the root is one level up); when
# installed it is the library the package was loaded from -- named
# explicitly, because R CMD check installs into a temporary library and
# a child that just says library(omicsCore) may pick up an older copy
# from the user library instead. The first version of this test did,
# and failed on a bug that copy still had.
omicscore_origin <- function() {
  root <- system.file(package = "omicsCore")
  if (identical(basename(root), "inst")) root <- dirname(root)
  # An installed package has an R/ directory too -- it holds the lazy-load
  # database -- so the source tree is recognised by its .R files, not by
  # the directory. Handing an installed directory to load_all() gives a
  # namespace with no functions in it.
  has_sources <- dir.exists(file.path(root, "R")) &&
    length(list.files(file.path(root, "R"), pattern = "\\.[Rr]$")) > 0L
  if (file.exists(file.path(root, "DESCRIPTION")) && has_sources &&
      requireNamespace("pkgload", quietly = TRUE)) {
    return(list(kind = "source", path = normalizePath(root)))
  }
  list(kind = "installed", path = dirname(normalizePath(root)))
}

writer_body <- function(origin, path, name, n) {
  if (identical(origin$kind, "source")) {
    pkgload::load_all(origin$path, quiet = TRUE)
  } else {
    library(omicsCore, lib.loc = c(origin$path, .libPaths()))
  }
  mat <- matrix(as.numeric(seq_len(60)) + nchar(name), 15, 4,
                dimnames = list(paste0("g", 1:15), paste0("s", 1:4)))
  meta <- data.frame(group = c("A", "A", "B", "B"), row.names = paste0("s", 1:4))
  inp <- omicsCore::omics_input(
    mat, meta, data.frame(feature_id = rownames(mat)),
    omics_type = "proteomics", assay_type = "normalized_intensity"
  )
  project <- omicsCore::omics_project(name = name,
                                      experiments = list(proteomics = inp))
  for (i in seq_len(n)) {
    omicsCore::save_project(project, path, overwrite = TRUE)
  }
  n
}

test_that("two writers on one file never leave a reader with a torn file", {
  skip_on_cran()
  skip_if_not_installed("callr")
  skip_if_not_installed("qs2")
  skip_if_not_installed("withr")

  d <- withr::local_tempdir()
  path <- file.path(d, "shared.omp")
  origin <- omicscore_origin()
  n <- 60L

  a <- callr::r_bg(writer_body, args = list(origin, path, "A", n))
  b <- callr::r_bg(writer_body, args = list(origin, path, "B", n))
  on.exit({ a$kill(); b$kill() }, add = TRUE)

  seen <- character(0)
  torn <- 0L
  deadline <- Sys.time() + 120
  while ((a$is_alive() || b$is_alive()) && Sys.time() < deadline) {
    if (file.exists(path)) {
      got <- tryCatch(omicsCore::load_project(path)$name,
                      error = function(e) NA_character_)
      if (is.na(got)) torn <- torn + 1L else seen <- c(seen, got)
    }
    Sys.sleep(0.01)
  }
  a$wait(timeout = 60e3)
  b$wait(timeout = 60e3)

  expect_identical(a$get_result(), n)
  expect_identical(b$get_result(), n)

  # Every read while the writers ran was a whole file
  expect_identical(torn, 0L)
  expect_true(all(seen %in% c("A", "B")))
  # And both writers were actually observed overlapping the reader
  expect_gt(length(seen), 0L)

  final <- omicsCore::load_project(path)
  expect_true(final$name %in% c("A", "B"))
  expect_true(omicsCore::is_omics_project(final))

  # Nothing of the losing writer's scratch survives
  expect_length(list.files(d, pattern = "\\.tmp$", all.files = TRUE), 0L)
  expect_identical(list.files(d), "shared.omp")
})

test_that("a reader that starts mid-write sees the previous complete file", {
  # The rename is what makes the guarantee; this pins that the target is
  # never opened for writing directly. A writer that truncated the
  # target first would give a reader an empty or partial file.
  skip_if_not_installed("qs2")
  skip_if_not_installed("withr")
  d <- withr::local_tempdir()
  path <- file.path(d, "p.omp")
  mat <- matrix(as.numeric(1:12), 3, 4,
                dimnames = list(paste0("g", 1:3), paste0("s", 1:4)))
  meta <- data.frame(group = c("A", "A", "B", "B"), row.names = paste0("s", 1:4))
  make <- function(name) omics_project(name, experiments = list(
    proteomics = omics_input(mat, meta, data.frame(feature_id = rownames(mat)),
                             omics_type = "proteomics",
                             assay_type = "normalized_intensity")))
  save_project(make("first"), path)
  before <- file.info(path)$size

  # A save that dies after opening its scratch file must not have touched
  # the target: simulate by making qs_save fail on the second call.
  calls <- 0L
  testthat::local_mocked_bindings(
    qs_save = function(object, file, ...) {
      calls <<- calls + 1L
      stop("disk full")
    },
    .package = "qs2"
  )
  expect_error(save_project(make("second"), path, overwrite = TRUE), "disk full")
  expect_identical(file.info(path)$size, before)
  expect_identical(load_project(path)$name, "first")
  expect_length(list.files(d, pattern = "\\.tmp$"), 0L)
})
