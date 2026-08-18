# Unit tests for the per-user project store (R/project_store.R).
#
# Every test points OMICSAPP_DATA_DIR at a fresh temp directory so the
# suite never touches a real deployment's `/data` mount, and restores
# the previous environment on exit.

with_temp_store <- function(code, quota_gb = NULL) {
  dir <- file.path(tempdir(), paste0("store-", as.integer(runif(1, 1, 1e9))))
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  old_dir <- Sys.getenv("OMICSAPP_DATA_DIR", unset = NA)
  old_quota <- Sys.getenv("OMICSAPP_QUOTA_GB", unset = NA)
  Sys.setenv(OMICSAPP_DATA_DIR = dir)
  if (is.null(quota_gb)) {
    Sys.unsetenv("OMICSAPP_QUOTA_GB")
  } else {
    Sys.setenv(OMICSAPP_QUOTA_GB = as.character(quota_gb))
  }
  on.exit({
    if (is.na(old_dir)) Sys.unsetenv("OMICSAPP_DATA_DIR")
    else Sys.setenv(OMICSAPP_DATA_DIR = old_dir)
    if (is.na(old_quota)) Sys.unsetenv("OMICSAPP_QUOTA_GB")
    else Sys.setenv(OMICSAPP_QUOTA_GB = old_quota)
    unlink(dir, recursive = TRUE)
  }, add = TRUE)
  force(code)
}

tiny_project <- function(name = "unit-test") {
  mat <- matrix(as.numeric(1:12), nrow = 3,
                dimnames = list(paste0("g", 1:3), paste0("s", 1:4)))
  meta <- data.frame(group = c("A", "A", "B", "B"),
                     row.names = paste0("s", 1:4))
  feat <- data.frame(feature_id = paste0("g", 1:3))
  inp <- omicsCore::omics_input(mat, meta, feat, omics_type = "proteomics")
  omicsCore::omics_project(name = name, experiments = list(proteomics = inp))
}

# ---- slug sanitisation ------------------------------------------------

test_that("project_slug keeps readable names and collapses separators", {
  expect_equal(project_slug("cheek G2 vs G1"), "cheek_G2_vs_G1")
  expect_equal(project_slug("  padded  "), "padded")
  expect_equal(project_slug("keep-dashes_and_underscores"),
               "keep-dashes_and_underscores")
})

test_that("project_slug neutralises path traversal", {
  for (evil in c("../../etc/passwd", "/etc/passwd", "..\\..\\windows",
                 "a/b/c")) {
    slug <- project_slug(evil)
    expect_false(is.na(slug))
    expect_false(grepl("[/\\]", slug))
    expect_false(grepl("\\.\\.", slug))
  }
})

test_that("project_slug refuses names that resolve to nothing usable", {
  expect_true(is.na(project_slug("")))
  expect_true(is.na(project_slug("   ")))
  expect_true(is.na(project_slug("...")))
  expect_true(is.na(project_slug(NA_character_)))
  expect_true(is.na(project_slug(NULL)))
})

test_that("no user-supplied name can reach the autosave slot", {
  # The invariant that matters is the resulting *path*, not the slug:
  # leading underscores are stripped, so "_autosave" lands on
  # "autosave.omp" — a different, harmless file.
  with_temp_store({
    adversarial <- c("_autosave", "__autosave__", "_autosave.omp",
                     "../_autosave", "/_autosave", ".autosave", "autosave")
    for (name in adversarial) {
      slug <- project_slug(name)
      if (is.na(slug)) next
      expect_false(normalizePath(project_path(slug), mustWork = FALSE) ==
                     normalizePath(autosave_path(), mustWork = FALSE),
                   info = name)
    }
  })
})

test_that("project_slug caps very long names", {
  expect_lte(nchar(project_slug(strrep("a", 500))), 80L)
})

test_that("project_slug strips control characters", {
  expect_equal(project_slug(paste0("ctl", "\x01\x02", "name")), "ctlname")
})

test_that("project_slug preserves non-Latin names in any locale", {
  # Regression guard: POSIX classes such as [[:cntrl:]] match byte-wise
  # under LC_CTYPE=C and shred multi-byte characters. Container base
  # images do not reliably set a UTF-8 locale, so this must hold
  # regardless of where the app runs. The literals are \u escapes so the
  # test source itself is locale-independent (and stays ASCII).
  slug <- project_slug("\u86cb\u767d\u7ec4 G2 vs G1")
  expect_equal(slug, "\u86cb\u767d\u7ec4_G2_vs_G1")
  expect_false(is.na(iconv(slug, "UTF-8", "UTF-8")))
})

test_that("truncating a long non-Latin name still yields valid UTF-8", {
  slug <- project_slug(strrep("\u86cb", 200))
  expect_false(is.na(iconv(slug, "UTF-8", "UTF-8")))
  expect_lte(nchar(slug, type = "bytes"), 80L)
})

# ---- quota ------------------------------------------------------------

test_that("quota is unlimited when unset or malformed", {
  with_temp_store(expect_identical(omicsapp_quota_gb(), Inf))
  with_temp_store(expect_identical(omicsapp_quota_gb(), Inf), quota_gb = "abc")
  with_temp_store(expect_identical(omicsapp_quota_gb(), Inf), quota_gb = "-5")
  with_temp_store(expect_identical(omicsapp_quota_gb(), Inf), quota_gb = "0")
})

test_that("quota reads a valid value", {
  with_temp_store(expect_identical(omicsapp_quota_gb(), 20), quota_gb = 20)
})

test_that("usage starts at zero and grows with saved files", {
  with_temp_store({
    expect_identical(data_dir_usage_bytes(), 0)
    writeBin(raw(4096), file.path(omicsapp_data_dir(), "blob.bin"))
    expect_gte(data_dir_usage_bytes(), 4096)
  })
})

test_that("quota_exceeded fires only once the directory is over budget", {
  with_temp_store({
    expect_false(quota_exceeded())
  }, quota_gb = 100)
  # A quota small enough that any file exceeds it.
  with_temp_store({
    writeBin(raw(2048), file.path(omicsapp_data_dir(), "blob.bin"))
    expect_true(quota_exceeded())
  }, quota_gb = 1e-9)
})

# ---- save / load / delete round trip ---------------------------------

test_that("a project round-trips through the store", {
  skip_if_not_installed("qs2")
  with_temp_store({
    res <- store_save_project(tiny_project(), slug = "roundtrip")
    expect_true(res$ok)
    expect_true(file.exists(res$path))

    back <- store_load_project("roundtrip")
    expect_true(back$ok)
    expect_true(omicsCore::is_omics_project(back$project))
    expect_equal(back$project$name, "unit-test")
  })
})

test_that("saving refuses to clobber without overwrite", {
  skip_if_not_installed("qs2")
  with_temp_store({
    expect_true(store_save_project(tiny_project(), slug = "dup")$ok)
    again <- store_save_project(tiny_project(), slug = "dup")
    expect_false(again$ok)
    expect_match(again$message, "already exists")
    expect_true(store_save_project(tiny_project(), slug = "dup",
                                   overwrite = TRUE)$ok)
  })
})

test_that("saving rejects an unusable name and a non-project", {
  with_temp_store({
    bad_name <- store_save_project(tiny_project(), slug = project_slug(""))
    expect_false(bad_name$ok)
    expect_match(bad_name$message, "project name")

    not_project <- store_save_project(list(a = 1), slug = "nope")
    expect_false(not_project$ok)
    expect_match(not_project$message, "Nothing to save")
  })
})

test_that("a new save is blocked at quota but an overwrite is not", {
  skip_if_not_installed("qs2")
  with_temp_store({
    expect_true(store_save_project(tiny_project(), slug = "first")$ok)
    # The directory now holds a file, so the tiny quota is exceeded.
    blocked <- store_save_project(tiny_project(), slug = "second")
    expect_false(blocked$ok)
    expect_match(blocked$message, "quota")
    # Replacing an existing project cannot grow the store without bound.
    expect_true(store_save_project(tiny_project(), slug = "first",
                                   overwrite = TRUE)$ok)
  }, quota_gb = 1e-9)
})

test_that("delete removes the file and reports a missing one", {
  skip_if_not_installed("qs2")
  with_temp_store({
    store_save_project(tiny_project(), slug = "goner")
    res <- store_delete_project("goner")
    expect_true(res$ok)
    expect_false(file.exists(project_path("goner")))
    expect_false(store_delete_project("goner")$ok)
  })
})

test_that("loading a missing project fails cleanly", {
  with_temp_store({
    res <- store_load_project("never-existed")
    expect_false(res$ok)
    expect_null(res$project)
    expect_match(res$message, "no longer exists")
  })
})

# ---- listing ----------------------------------------------------------

test_that("list_saved_projects is empty for a fresh store", {
  with_temp_store({
    expect_equal(nrow(list_saved_projects()), 0L)
  })
})

test_that("list_saved_projects hides the autosave snapshot", {
  skip_if_not_installed("qs2")
  with_temp_store({
    store_save_project(tiny_project(), slug = "visible")
    expect_true(store_autosave(tiny_project()))
    listed <- list_saved_projects()
    expect_equal(listed$slug, "visible")
    expect_false(paste0(AUTOSAVE_SLUG, ".omp") %in% basename(listed$path))
  })
})

# ---- autosave ---------------------------------------------------------

test_that("autosave round-trips and reports its timestamp", {
  skip_if_not_installed("qs2")
  with_temp_store({
    expect_null(store_read_autosave())
    expect_null(autosave_mtime())

    expect_true(store_autosave(tiny_project("recovered")))
    restored <- store_read_autosave()
    expect_true(omicsCore::is_omics_project(restored))
    expect_equal(restored$name, "recovered")
    expect_s3_class(autosave_mtime(), "POSIXct")
  })
})

test_that("autosave ignores non-projects and never throws", {
  with_temp_store({
    expect_false(store_autosave(NULL))
    expect_false(store_autosave(list(a = 1)))
  })
})

test_that("autosave still writes when the store is over quota", {
  skip_if_not_installed("qs2")
  with_temp_store({
    writeBin(raw(2048), file.path(omicsapp_data_dir(), "ballast.bin"))
    expect_true(quota_exceeded())
    # The snapshot replaces itself, so it must not be gated by quota —
    # otherwise recovery dies exactly when the user needs it.
    expect_true(store_autosave(tiny_project()))
  }, quota_gb = 1e-9)
})

# ---- raw upload archive -----------------------------------------------

write_upload <- function(content = "a,b\n1,2\n") {
  path <- tempfile(fileext = ".csv")
  writeLines(content, path)
  path
}

test_that("an upload is archived under the raw directory", {
  with_temp_store({
    src <- write_upload(); on.exit(unlink(src), add = TRUE)
    res <- store_raw_upload(src, "cheek.csv", "abc123def456789:proteomics:x")
    expect_true(res$ok)
    expect_true(file.exists(res$path))
    expect_equal(dirname(res$path), raw_dir())
    expect_match(basename(res$path), "^cheek__abc123def456\\.csv$")
  })
})

test_that("archiving the same upload twice does not duplicate it", {
  with_temp_store({
    src <- write_upload(); on.exit(unlink(src), add = TRUE)
    fp <- "abc123def456789:proteomics:x"
    first <- store_raw_upload(src, "cheek.csv", fp)
    second <- store_raw_upload(src, "cheek.csv", fp)
    expect_true(second$ok)
    expect_equal(second$path, first$path)
    expect_length(list.files(raw_dir()), 1L)
  })
})

test_that("archiving is skipped at quota but reports why", {
  with_temp_store({
    src <- write_upload(); on.exit(unlink(src), add = TRUE)
    writeBin(raw(2048), file.path(omicsapp_data_dir(), "ballast.bin"))
    res <- store_raw_upload(src, "cheek.csv", "deadbeefcafe:proteomics:x")
    expect_false(res$ok)
    expect_match(res$message, "quota")
  }, quota_gb = 1e-9)
})

test_that("archiving fails soft on a missing file or fingerprint", {
  with_temp_store({
    src <- write_upload(); on.exit(unlink(src), add = TRUE)
    expect_false(store_raw_upload(tempfile(), "x.csv", "fp")$ok)
    expect_false(store_raw_upload(src, "x.csv", NULL)$ok)
    # An import must never be blocked by the archive step.
    expect_no_error(store_raw_upload(NULL, NULL, NULL))
  })
})

test_that("archived uploads do not show up as saved projects", {
  skip_if_not_installed("qs2")
  with_temp_store({
    src <- write_upload(); on.exit(unlink(src), add = TRUE)
    store_raw_upload(src, "cheek.csv", "abc123def456:proteomics:x")
    store_save_project(tiny_project(), slug = "real-project")
    expect_equal(list_saved_projects()$slug, "real-project")
    # ...but they do count against the quota, since they really do grow
    # the store.
    expect_gt(data_dir_usage_bytes(), 0)
  })
})

# ---- autosave wiring --------------------------------------------------
#
# `wire_autosave()` is the only path deciding whether a recycled
# container costs the user their work, so it is driven here through a
# real reactive graph rather than by calling `store_autosave()` directly.

test_that("a project change reaches the store", {
  skip_if_not_installed("qs2")
  with_temp_store({
    proj <- tiny_project("recovered-by-wiring")
    shiny::testServer(
      function(input, output, session) {
        current_project <- shiny::reactiveVal(NULL)
        wire_autosave(current_project)
        shiny::observeEvent(input$import, current_project(proj))
      },
      {
        expect_null(store_read_autosave())
        session$setInputs(import = 1)
        restored <- store_read_autosave()
        expect_true(omicsCore::is_omics_project(restored))
        expect_equal(restored$name, "recovered-by-wiring")
      }
    )
  })
})

test_that("the write lands in the same flush, with no clock to advance", {
  skip_if_not_installed("qs2")
  with_temp_store({
    proj <- tiny_project("no-window")
    shiny::testServer(
      function(input, output, session) {
        current_project <- shiny::reactiveVal(NULL)
        wire_autosave(current_project)
        shiny::observeEvent(input$import, current_project(proj))
      },
      {
        session$setInputs(import = 1)
        # No session$elapse(). A debounce here would leave the change in
        # memory for its whole window, and ShinyProxy stops containers
        # with a signal that runs no R handler -- anything still only in
        # memory at that moment is gone.
        expect_false(is.null(store_read_autosave()))
      }
    )
  })
})

test_that("a cascade leaves the latest state on disk", {
  written <- list()
  recording_writer <- function(project) {
    written[[length(written) + 1L]] <<- project
    TRUE
  }
  first <- tiny_project("just-imported")
  second <- tiny_project("with-bundles")
  shiny::testServer(
    function(input, output, session) {
      current_project <- shiny::reactiveVal(NULL)
      wire_autosave(current_project, writer = recording_writer)
      # Mimics app_server(): importing adds the experiment, then the
      # bundle-attach observer replaces the project again.
      shiny::observeEvent(input$import, {
        current_project(first)
        current_project(second)
      })
    },
    {
      session$setInputs(import = 1)
      # Intermediate states may be written -- they are immediately
      # superseded. What matters is that the last word is the newest
      # state, never a stale one.
      expect_gte(length(written), 1L)
      expect_equal(written[[length(written)]]$name, "with-bundles")
    }
  )
})

test_that("a NULL project is never autosaved", {
  calls <- 0L
  counting_writer <- function(project) {
    calls <<- calls + 1L
    TRUE
  }
  shiny::testServer(
    function(input, output, session) {
      current_project <- shiny::reactiveVal(NULL)
      wire_autosave(current_project, writer = counting_writer)
    },
    {
      expect_equal(calls, 0L)
    }
  )
})
