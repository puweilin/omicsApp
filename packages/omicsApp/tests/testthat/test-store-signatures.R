# The store opens only files it wrote.
#
# Projects are qs2 files, and qs2 will build whatever R object a file
# describes. The store signs everything it saves and verifies before it
# deserialises; a store that predates signing is migrated once, signing
# the files that read as plain data; after that, an unsigned file is one
# something else put there, and it is refused. See the signatures block
# in R/project_store.R.

local_store <- function(env = parent.frame(), key = NA) {
  dir <- withr::local_tempdir(.local_envir = env)
  withr::local_envvar(OMICSAPP_DATA_DIR = dir, OMICSAPP_QUOTA_GB = NA,
                      OMICSAPP_SIGNING_KEY = key, .local_envir = env)
  dir
}

sig_project <- function(name = "sig") {
  mat <- matrix(as.numeric(1:12), nrow = 3,
                dimnames = list(paste0("g", 1:3), paste0("s", 1:4)))
  meta <- data.frame(group = c("A", "A", "B", "B"), row.names = paste0("s", 1:4))
  feat <- data.frame(feature_id = paste0("g", 1:3))
  inp <- omicsCore::omics_input(mat, meta, feat, omics_type = "proteomics",
                                assay_type = "normalized_intensity",
                                source_fingerprint = "abcdef123456:proteomics:x")
  omicsCore::omics_project(name = name, experiments = list(proteomics = inp))
}

ends_signed <- function(path) {
  bytes <- readBin(path, "raw", file.size(path))
  identical(rawToChar(utils::tail(bytes, 8L)), "OMPSIG01")
}

test_that("a saved project is signed and opens again", {
  skip_if_not_installed("openssl")
  dir <- local_store()
  res <- store_save_project(sig_project(), "mine", dir = dir)
  expect_true(res$ok)
  expect_true(ends_signed(res$path))
  back <- store_load_project("mine", dir = dir)
  expect_true(back$ok)
  expect_identical(back$project$name, "sig")
  # The generated key is kept in the store, readable by its owner only.
  key_file <- file.path(dir, SIGNING_KEY_FILE)
  expect_true(file.exists(key_file))
  if (.Platform$OS.type == "unix") {
    expect_identical(format(file.info(key_file)$mode), "600")
  }
  expect_true(file.exists(file.path(dir, SIGNING_MARKER_FILE)))
})

test_that("a tampered project is refused with a plain message", {
  skip_if_not_installed("openssl")
  dir <- local_store()
  path <- store_save_project(sig_project(), "mine", dir = dir)$path
  bytes <- readBin(path, "raw", file.size(path))
  bytes[60] <- xor(bytes[60], as.raw(0xff))
  writeBin(bytes, path)
  res <- store_load_project("mine", dir = dir)
  expect_false(res$ok)
  expect_null(res$project)
  expect_match(res$message, "changed since this app saved it")
  expect_false(grepl("qs2|HMAC|deserial", res$message, ignore.case = TRUE))
})

test_that("an unsigned file dropped into a signed store is refused", {
  skip_if_not_installed("openssl")
  dir <- local_store()
  store_save_project(sig_project(), "mine", dir = dir)     # store now migrated
  omicsCore::save_project(sig_project("foreign"), file.path(dir, "foreign.omp"))
  res <- store_load_project("foreign", dir = dir)
  expect_false(res$ok)
  expect_match(res$message, "not saved by this app")
  # Nor does it sneak in as a session to restore.
  omicsCore::save_project(sig_project("snap"), file.path(dir, "_autosave-zzz.omp"))
  expect_null(store_read_autosave(dir, path = file.path(dir, "_autosave-zzz.omp")))
})

test_that("a project carrying a closure is never signed or opened", {
  skip_if_not_installed("openssl")
  dir <- local_store()
  p <- sig_project("trojan")
  p$metadata$on_open <- function() stop("this ran")
  # Written as a pre-signing store would have held it.
  omicsCore::save_project(p, file.path(dir, "trojan.omp"))
  res <- store_load_project("trojan", dir = dir)
  expect_false(res$ok)
  expect_false(ends_signed(file.path(dir, "trojan.omp")))
  expect_match(res$message, "not saved by this app")
})

test_that("an existing store's unsigned projects are signed once and keep opening", {
  skip_if_not_installed("openssl")
  dir <- local_store()
  # A store written by a release before signing: plain save_project().
  omicsCore::save_project(sig_project("old"), file.path(dir, "old.omp"))
  omicsCore::save_project(sig_project("snap"), autosave_path(dir, "abc"))
  before <- readBin(file.path(dir, "old.omp"), "raw", file.size(file.path(dir, "old.omp")))
  expect_false(file.exists(file.path(dir, SIGNING_MARKER_FILE)))

  res <- store_load_project("old", dir = dir)
  expect_true(res$ok)
  expect_identical(res$project$name, "old")
  expect_true(ends_signed(file.path(dir, "old.omp")))
  expect_true(ends_signed(autosave_path(dir, "abc")))
  expect_identical(store_read_autosave(dir)$name, "snap")
  # The project itself was not rewritten, only signed.
  after <- readBin(file.path(dir, "old.omp"), "raw", file.size(file.path(dir, "old.omp")))
  expect_identical(after[seq_along(before)], before)
  expect_true(file.exists(file.path(dir, SIGNING_MARKER_FILE)))

  # Deleting the marker re-runs the migration, which is how a store gets
  # back files written by an older release after a rollback.
  omicsCore::save_project(sig_project("rolled-back"), file.path(dir, "rb.omp"))
  expect_false(store_load_project("rb", dir = dir)$ok)
  unlink(file.path(dir, SIGNING_MARKER_FILE))
  expect_true(store_load_project("rb", dir = dir)$ok)
})

test_that("setting OMICSAPP_SIGNING_KEY later moves the store to it", {
  skip_if_not_installed("openssl")
  dir <- local_store()
  store_save_project(sig_project(), "mine", dir = dir)
  expect_true(file.exists(file.path(dir, SIGNING_KEY_FILE)))

  withr::local_envvar(OMICSAPP_SIGNING_KEY = "deployment-key-0123456789abcdef0123456789")
  res <- store_load_project("mine", dir = dir)
  expect_true(res$ok)
  # The generated key is retired once nothing needs it.
  expect_false(file.exists(file.path(dir, SIGNING_KEY_FILE)))
  expect_true(is.list(omicsCore::load_project(
    file.path(dir, "mine.omp"),
    signing_key = "deployment-key-0123456789abcdef0123456789")))
})

test_that("autosave snapshots are signed and verified too", {
  skip_if_not_installed("openssl")
  dir <- local_store(key = "k-0123456789abcdef0123456789abcdef")
  expect_true(store_autosave(sig_project("snap"), dir = dir, id = "s1"))
  expect_true(ends_signed(autosave_path(dir, "s1")))
  expect_identical(store_read_autosave(dir)$name, "snap")
  path <- autosave_path(dir, "s1")
  bytes <- readBin(path, "raw", file.size(path))
  bytes[50] <- xor(bytes[50], as.raw(1))
  writeBin(bytes, path)
  expect_null(store_read_autosave(dir))
})

test_that("an upload is kept while a project that might use it will not open", {
  skip_if_not_installed("openssl")
  dir <- local_store()
  store_save_project(sig_project("a"), "a", dir = dir)
  raw <- raw_dir(dir)
  upload <- file.path(raw, "data__abcdef123456.xlsx")
  writeLines("x", upload)
  # Another project, unreadable: it may refer to the same upload.
  writeLines("not a project", file.path(dir, "b.omp"))
  res <- store_delete_project("a", dir = dir)
  expect_true(res$ok)
  expect_true(file.exists(upload))
})
