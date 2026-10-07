# Which omicsCore the suite runs against (see setup.R): the source tree
# beside omicsApp for a development run, the installed packages for a run
# against an installed omicsApp. setup.R writes its answer to
# OMICSAPP_TEST_CORE; the omicsCore actually loaded has to be that one,
# or the suite is not testing what it says it is.

test_that("the omicsCore loaded is the one setup.R settled on", {
  skip_if_not_installed("pkgload")
  mode <- Sys.getenv("OMICSAPP_TEST_CORE")
  expect_true(mode %in% c("installed", "source"))
  expect_identical(pkgload::is_dev_package("omicsCore"), identical(mode, "source"))
})

test_that("an installed omicsApp is tested against installed packages only", {
  skip_if_not_installed("pkgload")
  # A development run loads omicsApp from source, and this does not apply.
  skip_if(pkgload::is_dev_package("omicsApp"), "omicsApp is loaded from source")
  expect_identical(Sys.getenv("OMICSAPP_TEST_CORE"), "installed")
  expect_false(pkgload::is_dev_package("omicsCore"))
})

test_that("the browser tests run the app the suite tests", {
  where <- smoke_app_dir()
  from_source <- requireNamespace("pkgload", quietly = TRUE) &&
    pkgload::is_dev_package("omicsApp")
  # From source exactly when omicsApp is; app.R then leaves omicsCore
  # installed when setup.R chose the installed one.
  expect_identical(nzchar(where$dev_root), from_source)
  if (from_source) {
    expect_true(file.exists(file.path(where$dev_root, "omicsApp", "DESCRIPTION")))
  }
})
