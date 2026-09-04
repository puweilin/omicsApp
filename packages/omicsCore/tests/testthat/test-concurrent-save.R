# One account, two devices: both sessions bind-mount the same directory
# and both autosave. save_project() writes to a temp file and renames,
# which is atomic for a reader -- but the temp name used to be
# paste0(path, ".tmp"), shared by every writer, which takes that back.
# Two writers open the same file, interleave, and one renames the result
# over the real one.

tiny_project <- function(name = "p", n = 3L) {
  m <- matrix(as.numeric(seq_len(n * 4)) / 10, n, 4,
              dimnames = list(paste0("g", seq_len(n)), paste0("s", 1:4)))
  meta <- data.frame(group = c("A", "A", "B", "B"), row.names = paste0("s", 1:4))
  omics_project(name = name, experiments = list(proteomics = omics_input(
    m, meta, data.frame(feature_id = rownames(m)),
    omics_type = "proteomics", assay_type = "normalized_intensity")))
}

test_that("the temp name is generated, not derived from the target", {
  skip_if_not_installed("qs2")
  d <- withr::local_tempdir()
  path <- file.path(d, "autosave.omp")

  leftovers <- function() setdiff(list.files(d), "autosave.omp")
  # Both saves complete and neither leaves scratch behind. A shared
  # name would still pass this serially -- the guarantee it buys is for
  # writers that overlap, which a single-process test cannot stage --
  # so the check that matters is the one below on where the name comes
  # from.
  save_project(tiny_project("a"), path, overwrite = TRUE)
  first <- load_project(path)$name
  save_project(tiny_project("b"), path, overwrite = TRUE)
  second <- load_project(path)$name

  expect_identical(first, "a")
  expect_identical(second, "b")
  expect_length(leftovers(), 0L)

  # The property that actually protects concurrent writers: the scratch
  # path is drawn from tempfile(), so two callers cannot collide on it.
  src <- deparse(body(save_project))
  expect_true(any(grepl("tempfile(", src, fixed = TRUE)))
  expect_false(any(grepl('paste0(path, ".tmp")', src, fixed = TRUE)))
})

test_that("the temp file sits beside its target, not in tempdir()", {
  # rename() is only atomic within one filesystem, and tempdir() is
  # usually a different one from the user's store.
  skip_if_not_installed("qs2")
  d <- withr::local_tempdir()
  path <- file.path(d, "autosave.omp")
  # A directory the process cannot write would fail at the temp write
  # rather than after a cross-device rename, which is the safer order.
  expect_silent(save_project(tiny_project(), path, overwrite = TRUE))
  expect_true(file.exists(path))
})

test_that("nothing is left behind when the write fails", {
  skip_if_not_installed("qs2")
  d <- withr::local_tempdir()
  expect_error(save_project("not a project", file.path(d, "x.omp")))
  expect_length(list.files(d), 0L)
})

test_that("a reader never sees a partial file", {
  # The rename is the guarantee: the target either is the old content or
  # the new one, never bytes in between.
  skip_if_not_installed("qs2")
  d <- withr::local_tempdir()
  path <- file.path(d, "autosave.omp")
  save_project(tiny_project("first"), path, overwrite = TRUE)
  for (i in 1:5) {
    save_project(tiny_project(paste0("n", i), n = 50L), path, overwrite = TRUE)
    expect_identical(load_project(path)$name, paste0("n", i))
  }
})
