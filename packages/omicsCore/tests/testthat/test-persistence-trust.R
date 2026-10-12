# Which .omp files load_project() will read, and how an old one is
# brought up to date.
#
# A qs2 file can hold any R object, so a reader that checks a file only
# after deserialising it has already let the file decide what to build.
# These pin down the three lines drawn against that: a signed file is
# verified before it is read, an unsigned one from elsewhere is checked
# for anything that is not data, and a file in an older format is
# upgraded by the registered migrations rather than misread.

trust_project <- function(name = "trust") {
  mat <- matrix(as.numeric(1:24), nrow = 6,
                dimnames = list(paste0("g", 1:6), paste0("s", 1:4)))
  meta <- data.frame(group = c("A", "A", "B", "B"), row.names = paste0("s", 1:4))
  feat <- data.frame(feature_id = paste0("g", 1:6))
  omics_project(name, experiments = list(
    proteomics = omics_input(mat, meta, feat, omics_type = "proteomics",
                             assay_type = "normalized_intensity")))
}

KEY <- "a-test-key-of-reasonable-length-0123456789"

# ---- signing -----------------------------------------------------------------

test_that("a signed project round-trips under its key and reads unchanged without one", {
  skip_if_not_installed("qs2")
  skip_if_not_installed("openssl")
  p <- trust_project()
  f <- withr::local_tempfile(fileext = ".omp")
  save_project(p, f, signing_key = KEY)
  back <- load_project(f, signing_key = KEY)
  expect_identical(back$experiments, p$experiments)
  expect_identical(back$name, p$name)
  # A reader that asks for no signature -- an older release after a
  # rollback, the restore drill -- still opens it: qs2 ignores the
  # trailer.
  expect_identical(load_project(f)$experiments, p$experiments)
  expect_identical(qs2::qs_read(f)$payload$experiments, p$experiments)
})

test_that("a tampered file is refused before it is read", {
  skip_if_not_installed("qs2")
  skip_if_not_installed("openssl")
  f <- withr::local_tempfile(fileext = ".omp")
  save_project(trust_project(), f, signing_key = KEY)
  bytes <- readBin(f, "raw", file.size(f))
  bytes[40] <- xor(bytes[40], as.raw(0x01))
  writeBin(bytes, f)
  # If the body were deserialised first, this would be qs2's error.
  testthat::local_mocked_bindings(
    qs_deserialize = function(...) stop("deserialised an unverified file"),
    .package = "qs2")
  err <- expect_error(load_project(f, signing_key = KEY), class = "omp_signature_error")
  expect_match(conditionMessage(err), "changed since this app saved it")
})

test_that("a file signed under another key is refused", {
  skip_if_not_installed("qs2")
  skip_if_not_installed("openssl")
  f <- withr::local_tempfile(fileext = ".omp")
  save_project(trust_project(), f, signing_key = "some-other-deployment's-key")
  expect_error(load_project(f, signing_key = KEY), class = "omp_signature_error")
})

test_that("an unsigned foreign file is refused when a signature is required", {
  skip_if_not_installed("qs2")
  skip_if_not_installed("openssl")
  f <- withr::local_tempfile(fileext = ".omp")
  save_project(trust_project(), f)
  err <- expect_error(load_project(f, signing_key = KEY), class = "omp_unsigned_error")
  expect_match(conditionMessage(err), "not saved by this app")
  # A file too short to carry a trailer is unsigned, not a crash.
  g <- withr::local_tempfile(fileext = ".omp")
  writeBin(as.raw(1:10), g)
  expect_error(load_project(g, signing_key = KEY), class = "omp_unsigned_error")
})

test_that("a stripped or truncated signature is refused", {
  skip_if_not_installed("qs2")
  skip_if_not_installed("openssl")
  f <- withr::local_tempfile(fileext = ".omp")
  save_project(trust_project(), f, signing_key = KEY)
  bytes <- readBin(f, "raw", file.size(f))
  writeBin(bytes[seq_len(length(bytes) - 1L)], f)
  expect_error(load_project(f, signing_key = KEY), class = "omp_untrusted_error")
})

test_that("the signing key is checked by name", {
  skip_if_not_installed("qs2")
  f <- withr::local_tempfile(fileext = ".omp")
  for (bad in list(NA, "", 1, character(0), list("k"))) {
    expect_error(save_project(trust_project(), f, overwrite = TRUE, signing_key = bad),
                 "`signing_key` must be NULL", fixed = TRUE)
  }
  save_project(trust_project(), f, overwrite = TRUE)
  expect_error(load_project(f, untrusted = NA), "`untrusted` must be TRUE or FALSE", fixed = TRUE)
})

# ---- what an untrusted file may hold --------------------------------------------

write_raw_omp <- function(payload) {
  f <- tempfile(fileext = ".omp")
  qs2::qs_save(list(schema_version = OMP_SCHEMA_VERSION, payload = payload), f)
  f
}

test_that("an untrusted file carrying a function is refused, and names where", {
  skip_if_not_installed("qs2")
  p <- trust_project()
  p$metadata$hook <- function() system("echo pwned")
  f <- write_raw_omp(p)
  on.exit(unlink(f), add = TRUE)
  err <- expect_error(load_project(f, untrusted = TRUE), class = "omp_unsafe_error")
  expect_match(conditionMessage(err), "a function")
  expect_match(conditionMessage(err), "metadata$hook", fixed = TRUE)
  # The trusted path is unchanged: what the caller asked for is what runs.
  expect_true(is_omics_project(load_project(f)))
})

test_that("environments, external pointers and code are refused wherever they sit", {
  skip_if_not_installed("qs2")
  cases <- list(
    environment = function(p) { p$experiments$proteomics$cache <- new.env(); p },
    attribute   = function(p) { attr(p$experiments$proteomics$expr_mat, "x") <- new.env(); p },
    pointer     = function(p) { p$metadata$ptr <- methods::new("externalptr"); p },
    formula     = function(p) { p$metadata$f <- y ~ x; p },
    call        = function(p) { p$metadata$cl <- quote(system("id")); p },
    s3_class    = function(p) { p$metadata$x <- structure(list(), class = "evil"); p },
    bundle      = function(p) {
      p$bundles <- list(diff = structure(list(results = list(fit = mean)),
                                         class = "analysis_bundle"))
      p
    }
  )
  for (nm in names(cases)) {
    f <- write_raw_omp(cases[[nm]](trust_project()))
    expect_error(load_project(f, untrusted = TRUE), class = "omp_unsafe_error", info = nm)
    unlink(f)
  }
})

test_that("a project of real analyses passes the structure check", {
  skip_if_not_installed("qs2")
  skip_if_not_installed("limma")
  inp <- realistic_input("proteomics")
  p <- omics_project("real", experiments = list(proteomics = inp))
  p$bundles <- list(
    qc = run_qc(inp),
    diff = run_diff(inp, method = "limma", analysis_type = "group", group_col = "group",
                    control_group = "G1", case_group = "G2"),
    ttest = run_diff(inp, method = "ttest", analysis_type = "group", group_col = "group",
                     control_group = "G1", case_group = "G2"))
  p$visited_steps <- c("import", "qc", "diff")
  f <- withr::local_tempfile(fileext = ".omp")
  save_project(p, f)
  back <- load_project(f, untrusted = TRUE)
  expect_identical(back$bundles$diff$results$diff_result_df,
                   p$bundles$diff$results$diff_result_df)
})

test_that("absurdly deep nesting is refused rather than recursed into", {
  skip_if_not_installed("qs2")
  p <- trust_project()
  deep <- list()
  for (i in 1:100) deep <- list(deep)
  p$metadata$deep <- deep
  f <- write_raw_omp(p)
  on.exit(unlink(f), add = TRUE)
  expect_error(load_project(f, untrusted = TRUE), "nested too deeply")
})

# ---- adopting and re-signing files ------------------------------------------------

test_that("sign_project_file adopts a clean unsigned file without rewriting its archive", {
  skip_if_not_installed("qs2")
  skip_if_not_installed("openssl")
  f <- withr::local_tempfile(fileext = ".omp")
  save_project(trust_project(), f)
  before <- readBin(f, "raw", file.size(f))
  expect_identical(sign_project_file(f, KEY)$status, "refused")   # not without consent
  expect_identical(sign_project_file(f, KEY, adopt_unsigned = TRUE)$status, "signed")
  after <- readBin(f, "raw", file.size(f))
  expect_identical(after[seq_along(before)], before)
  expect_true(is_omics_project(load_project(f, signing_key = KEY)))
  expect_identical(sign_project_file(f, KEY, adopt_unsigned = TRUE)$status, "valid")
})

test_that("sign_project_file will not adopt a file carrying a closure", {
  skip_if_not_installed("qs2")
  skip_if_not_installed("openssl")
  p <- trust_project()
  p$metadata$hook <- function() NULL
  f <- write_raw_omp(p)
  on.exit(unlink(f), add = TRUE)
  before <- readBin(f, "raw", file.size(f))
  res <- sign_project_file(f, KEY, adopt_unsigned = TRUE)
  expect_identical(res$status, "refused")
  expect_match(res$message, "a function")
  expect_identical(readBin(f, "raw", file.size(f)), before)
  expect_error(load_project(f, signing_key = KEY), class = "omp_unsigned_error")
})

test_that("sign_project_file moves a file to a new key, and only from a key it was given", {
  skip_if_not_installed("qs2")
  skip_if_not_installed("openssl")
  f <- withr::local_tempfile(fileext = ".omp")
  save_project(trust_project(), f, signing_key = "old-key")
  expect_identical(sign_project_file(f, KEY)$status, "refused")
  expect_identical(sign_project_file(f, KEY, verify_key = "wrong")$status, "refused")
  expect_identical(sign_project_file(f, KEY, verify_key = "old-key")$status, "resigned")
  expect_true(is_omics_project(load_project(f, signing_key = KEY)))
  expect_error(load_project(f, signing_key = "old-key"), class = "omp_signature_error")
})

# ---- format versions and migrations -------------------------------------------------

test_that("saved projects carry the format version, and one without it is format 1", {
  skip_if_not_installed("qs2")
  f <- withr::local_tempfile(fileext = ".omp")
  save_project(trust_project(), f)
  env <- qs2::qs_read(f)
  expect_identical(env$format_version, OMP_FORMAT_VERSION)
  # Still readable by a release that only knows schema_version.
  expect_identical(env$schema_version, "1.0.0")
  # A file from before format_version existed loads unchanged.
  old <- write_raw_omp(trust_project("old"))
  on.exit(unlink(old), add = TRUE)
  expect_null(qs2::qs_read(old)$format_version)
  expect_identical(load_project(old)$name, "old")
})

test_that("an old-format project is upgraded by its registered migration", {
  skip_if_not_installed("qs2")
  # Pretend format 2 renamed project$name to project$title's successor:
  # a format-1 file is the current layout, and the migration adds a
  # field format 2 requires.
  testthat::local_mocked_bindings(current_omp_format_version = function() 2L)
  calls <- 0L
  prev <- register_omp_migration(1L, function(project) {
    calls <<- calls + 1L
    project$metadata$migrated_from <- 1L
    project
  })
  withr::defer(if (is.null(prev)) rm("1", envir = .omp_migrations)
               else assign("1", prev, envir = .omp_migrations))

  old <- write_raw_omp(trust_project("legacy"))      # no format_version: format 1
  on.exit(unlink(old), add = TRUE)
  back <- load_project(old)
  expect_identical(calls, 1L)
  expect_identical(back$metadata$migrated_from, 1L)
  expect_identical(back$name, "legacy")

  # A file written now is format 2 and is not migrated again.
  f <- withr::local_tempfile(fileext = ".omp")
  save_project(back, f)
  expect_identical(qs2::qs_read(f)$format_version, 2L)
  load_project(f)
  expect_identical(calls, 1L)
})

test_that("a missing migration and a newer format are both refused by name", {
  skip_if_not_installed("qs2")
  testthat::local_mocked_bindings(current_omp_format_version = function() 3L)
  old <- write_raw_omp(trust_project())
  on.exit(unlink(old), add = TRUE)
  expect_error(load_project(old), "No upgrade is available for project file format 1")

  testthat::local_mocked_bindings(current_omp_format_version = function() 1L)
  f <- tempfile(fileext = ".omp")
  on.exit(unlink(f), add = TRUE)
  qs2::qs_save(list(schema_version = "1.0.0", format_version = 5L,
                    payload = trust_project()), f)
  expect_error(load_project(f), "saved by a newer version of the app")
})
