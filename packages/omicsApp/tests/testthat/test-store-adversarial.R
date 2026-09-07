# The project store takes a filename from a text box and a quota from a
# deployment config, and writes to a directory shared by one user's whole
# working life. Both inputs are worth being unpleasant to.

adv_store <- function(code, quota_gb = NULL) {
  dir <- file.path(tempdir(), paste0("adv-", as.integer(runif(1, 1, 1e9))))
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  old_dir <- Sys.getenv("OMICSAPP_DATA_DIR", unset = NA)
  old_quota <- Sys.getenv("OMICSAPP_QUOTA_GB", unset = NA)
  Sys.setenv(OMICSAPP_DATA_DIR = dir)
  if (is.null(quota_gb)) Sys.unsetenv("OMICSAPP_QUOTA_GB")
  else Sys.setenv(OMICSAPP_QUOTA_GB = as.character(quota_gb))
  on.exit({
    if (is.na(old_dir)) Sys.unsetenv("OMICSAPP_DATA_DIR")
    else Sys.setenv(OMICSAPP_DATA_DIR = old_dir)
    if (is.na(old_quota)) Sys.unsetenv("OMICSAPP_QUOTA_GB")
    else Sys.setenv(OMICSAPP_QUOTA_GB = old_quota)
    unlink(dir, recursive = TRUE)
  }, add = TRUE)
  force(code)
}

adv_project <- function(name = "adv") {
  mat <- matrix(as.numeric(1:12), nrow = 3,
                dimnames = list(paste0("g", 1:3), paste0("s", 1:4)))
  meta <- data.frame(group = c("A", "A", "B", "B"),
                     row.names = paste0("s", 1:4))
  feat <- data.frame(feature_id = paste0("g", 1:3))
  inp <- omicsCore::omics_input(mat, meta, feat, omics_type = "proteomics", assay_type = "normalized_intensity")
  omicsCore::omics_project(name = name,
                           experiments = list(proteomics = inp))
}

# ---- names that try to leave the directory ----------------------------

test_that("no name escapes the store", {
  adv_store({
    root <- normalizePath(omicsapp_data_dir(), mustWork = TRUE)
    hostile <- c(
      "../../../../etc/passwd", "..\\..\\windows\\system32",
      "/etc/shadow", "....//....//etc", "a/../../b",
      "~/.ssh/authorized_keys", "$HOME/x", "%SYSTEMROOT%",
      ".hidden", "..", ".", "...", "./x", "con", "NUL"
    )
    for (name in hostile) {
      slug <- project_slug(name)
      if (is.na(slug)) next
      resolved <- normalizePath(dirname(project_path(slug)), mustWork = FALSE)
      expect_identical(resolved, root, info = name)
      expect_false(grepl("[/\\\\]", slug), info = name)
    }
  })
})

test_that("a name of only punctuation resolves to nothing usable", {
  for (name in c("...", "___", "///", "\\\\", "..", ". . .", "___...")) {
    expect_true(is.na(project_slug(name)), info = name)
  }
})

test_that("the stored filename never doubles a dot", {
  # `stem.` + ".omp" gives `stem..omp` -- the shape stripped elsewhere
  # for safety, and the length cap can land exactly on a dot.
  for (name in c(paste0(strrep("a", 79), ".tail"),
                 paste0(strrep("b", 80), "."),
                 "trailing.")) {
    slug <- project_slug(name)
    if (is.na(slug)) next
    expect_false(grepl("\\.\\.", basename(project_path(slug))), info = name)
  }
})

test_that("a name is capped in bytes, which is what filesystems count", {
  for (name in list(strrep("a", 500), strrep("\u86cb", 200),
                    paste0(strrep("\u86cb", 100), strrep("z", 100)))) {
    slug <- project_slug(name)
    expect_lte(nchar(slug, type = "bytes"), 80L)
    expect_false(is.na(iconv(slug, "UTF-8", "UTF-8")))
  }
})

# ---- quota -------------------------------------------------------------

test_that("a malformed quota means unlimited, not locked out", {
  # A typo in a deployment config must not leave a user unable to save
  # their own work.
  for (bad in c("abc", "", "-1", "0", "NaN", "1e400x", " ")) {
    adv_store(expect_identical(omicsapp_quota_gb(), Inf), quota_gb = bad)
  }
})

test_that("a fractional quota is honoured", {
  adv_store(expect_identical(omicsapp_quota_gb(), 0.5), quota_gb = "0.5")
})

test_that("usage counts the raw archive as well as the projects", {
  skip_if_not_installed("qs2")
  adv_store({
    store_save_project(adv_project(), slug = "p")
    before <- data_dir_usage_bytes()
    src <- tempfile(fileext = ".csv"); on.exit(unlink(src), add = TRUE)
    writeLines(strrep("x", 4096), src)
    store_raw_upload(src, "big.csv", "aabbccddeeff:proteomics:x")
    # Archived uploads live in a sub-directory; a usage figure that
    # ignored them would let the store grow past a quota it claims to
    # enforce.
    expect_gt(data_dir_usage_bytes(), before)
  })
})

test_that("a quota that is already exceeded still allows tidying up", {
  skip_if_not_installed("qs2")
  adv_store({
    expect_true(store_save_project(adv_project(), slug = "keep")$ok)
    expect_true(quota_exceeded())
    # Over quota, a user must still be able to overwrite and delete;
    # otherwise the only way out of a full directory is an admin.
    expect_true(store_save_project(adv_project(), slug = "keep",
                                   overwrite = TRUE)$ok)
    expect_true(store_delete_project("keep")$ok)
  }, quota_gb = 1e-9)
})

# ---- store integrity ---------------------------------------------------

test_that("a corrupt project file fails to open without taking the app down", {
  adv_store({
    writeLines("this is not a qs2 archive", project_path("broken"))
    res <- store_load_project("broken")
    expect_false(res$ok)
    expect_null(res$project)
    # ...and it still lists, so the user can see it and delete it.
    expect_true("broken" %in% list_saved_projects()$slug)
    expect_true(store_delete_project("broken")$ok)
  })
})

test_that("a foreign file in the store is ignored rather than listed", {
  adv_store({
    writeLines("notes", file.path(omicsapp_data_dir(), "README.txt"))
    dir.create(file.path(omicsapp_data_dir(), "subdir"))
    expect_equal(nrow(list_saved_projects()), 0L)
  })
})

test_that("saving into a store whose directory has vanished recreates it", {
  skip_if_not_installed("qs2")
  adv_store({
    unlink(omicsapp_data_dir(), recursive = TRUE)
    expect_true(store_save_project(adv_project(), slug = "revived")$ok)
  })
})

test_that("autosave survives the store directory vanishing mid-session", {
  skip_if_not_installed("qs2")
  adv_store({
    unlink(omicsapp_data_dir(), recursive = TRUE)
    expect_true(store_autosave(adv_project()))
    expect_true(omicsCore::is_omics_project(store_read_autosave()))
  })
})

# ---- raw archive -------------------------------------------------------

test_that("two files with the same name but different content both survive", {
  adv_store({
    a <- tempfile(); b <- tempfile()
    on.exit(unlink(c(a, b)), add = TRUE)
    writeLines("content A", a); writeLines("content B", b)
    ra <- store_raw_upload(a, "data.csv", "1111aaaa2222:proteomics:x")
    rb <- store_raw_upload(b, "data.csv", "3333bbbb4444:proteomics:x")
    expect_true(ra$ok); expect_true(rb$ok)
    # Named from the fingerprint, so the second does not silently
    # replace the first and leave a project pointing at the wrong file.
    expect_false(identical(ra$path, rb$path))
    expect_equal(length(list.files(raw_dir())), 2L)
  })
})

test_that("an upload with a hostile filename is archived safely", {
  adv_store({
    src <- tempfile(); on.exit(unlink(src), add = TRUE)
    writeLines("x", src)
    res <- store_raw_upload(src, "../../evil.xlsx", "abcdef123456:p:x")
    expect_true(res$ok)
    expect_identical(normalizePath(dirname(res$path)),
                     normalizePath(raw_dir()))
  })
})

# --- the slug is a client-supplied value --------------------------------

# The existing traversal case above uses `../../etc/passwd`, which addresses
# `passwd.omp` -- a file that does not exist, so "nothing happened" proved
# nothing about containment. These build a real .omp outside the store and
# assert it is still there afterwards.
test_that("a real .omp outside the store can be neither opened nor deleted", {
  skip_if_not_installed("withr")
  store   <- withr::local_tempdir()
  outside <- withr::local_tempdir()
  victim  <- file.path(outside, "victim.omp")
  omicsCore::save_project(omicsCore::omics_project("secret"), victim)

  rel <- file.path("..", basename(outside), "victim")
  open_res <- store_load_project(rel, dir = store)
  del_res  <- store_delete_project(rel, dir = store)

  expect_false(isTRUE(open_res$ok))
  expect_null(open_res$project)
  expect_false(isTRUE(del_res$ok))
  expect_true(file.exists(victim))
  # A path-shaped slug gets one fixed answer whether or not the target is
  # there, so a probe cannot map the filesystem by reading the message.
  expect_equal(open_res$message, "Unknown project.")
  expect_equal(del_res$message, "Unknown project.")
  absent <- store_load_project(file.path("..", "nothing-here", "x"), dir = store)
  expect_equal(absent$message, "Unknown project.")
})

test_that("an absolute path outside the store is refused", {
  skip_if_not_installed("withr")
  store   <- withr::local_tempdir()
  outside <- withr::local_tempdir()
  victim  <- file.path(outside, "victim.omp")
  omicsCore::save_project(omicsCore::omics_project("secret"), victim)

  abs_slug <- sub("\\.omp$", "", victim)
  expect_false(isTRUE(store_load_project(abs_slug, dir = store)$ok))
  expect_false(isTRUE(store_delete_project(abs_slug, dir = store)$ok))
  expect_true(file.exists(victim))
})

test_that("a symlink in the store pointing out of it is refused", {
  skip_if_not_installed("withr")
  skip_on_os("windows")
  store   <- withr::local_tempdir()
  outside <- withr::local_tempdir()
  victim  <- file.path(outside, "victim.omp")
  omicsCore::save_project(omicsCore::omics_project("secret"), victim)
  # Listed by the store, but resolves outside it.
  link_ok <- file.symlink(victim, file.path(store, "innocent.omp"))
  skip_if_not(link_ok, "could not create a symlink")

  expect_false(isTRUE(store_load_project("innocent", dir = store)$ok))
  expect_false(isTRUE(store_delete_project("innocent", dir = store)$ok))
  expect_true(file.exists(victim))
})

test_that("a normal project, including a CJK name, still opens and deletes", {
  skip_if_not_installed("withr")
  store <- withr::local_tempdir()
  for (nm in c("my project", "蛋白组项目")) {
    saved <- store_save_project(omicsCore::omics_project(nm),
                               slug = project_slug(nm), dir = store)
    expect_true(isTRUE(saved$ok), info = nm)
    slug <- project_slug(nm)
    opened <- store_load_project(slug, dir = store)
    expect_true(isTRUE(opened$ok), info = nm)
    expect_equal(opened$project$name, nm)
    expect_true(isTRUE(store_delete_project(slug, dir = store)$ok), info = nm)
  }
})

test_that("a slug that vanishes between listing and use is reported as such", {
  skip_if_not_installed("withr")
  store <- withr::local_tempdir()
  store_save_project(omicsCore::omics_project("gone"), slug = "gone", dir = store)
  path <- file.path(store, "gone.omp")
  # resolve_project_path() lists, then normalizes; simulate the file being
  # removed by another session in between.
  expect_true(isTRUE(resolve_project_path("gone", store)$ok))
  unlink(path)
  res <- resolve_project_path("gone", store)
  expect_false(isTRUE(res$ok))
  expect_match(res$message, "no longer exists")
})
