# The individual pieces are covered elsewhere. These are the seams
# between them, which is where the bugs found so far actually lived: the
# stale-bundle bug was not in either observer, it was in what the two of
# them did to each other.

fi_store <- function(code) {
  dir <- file.path(tempdir(), paste0("fi-", as.integer(runif(1, 1, 1e9))))
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  old <- Sys.getenv("OMICSAPP_DATA_DIR", unset = NA)
  Sys.setenv(OMICSAPP_DATA_DIR = dir)
  on.exit({
    if (is.na(old)) Sys.unsetenv("OMICSAPP_DATA_DIR")
    else Sys.setenv(OMICSAPP_DATA_DIR = old)
    unlink(dir, recursive = TRUE)
  }, add = TRUE)
  force(code)
}

fi_input <- function(fingerprint = NULL, source_path = NULL, n = 6L) {
  mat <- matrix(as.numeric(seq_len(n * 4L)), nrow = n,
                dimnames = list(paste0("g", seq_len(n)), paste0("s", 1:4)))
  meta <- data.frame(group = c("A", "A", "B", "B"),
                     row.names = paste0("s", 1:4))
  feat <- data.frame(feature_id = paste0("g", seq_len(n)),
                     feature_symbol = paste0("SYM", seq_len(n)))
  omicsCore::omics_input(mat, meta, feat, omics_type = "proteomics",
                         assay_type = "normalized_intensity",
                         source_fingerprint = fingerprint,
                         source_path = source_path)
}

fi_project <- function(input, bundles = list()) {
  proj <- omicsCore::omics_project("interaction",
                                   experiments = list(proteomics = input))
  proj$bundles <- bundles
  proj
}

# ---- autosave and the layer-replacement guard -------------------------

test_that("the autosave reflects a replacement rather than the state before it", {
  skip_if_not_installed("qs2")
  fi_store({
    before <- fi_project(fi_input("fp-A"), bundles = list(diff = "OLD"))
    after  <- fi_project(fi_input("fp-B"))
    shiny::testServer(
      function(input, output, session) {
        current_project <- shiny::reactiveVal(NULL)
        wire_autosave(current_project)
        shiny::observeEvent(input$first, current_project(before))
        shiny::observeEvent(input$replace, current_project(after))
      },
      {
        session$setInputs(first = 1)
        expect_length(store_read_autosave()$bundles, 1L)
        session$setInputs(replace = 1)
        # A recovery copy still holding results computed on data the
        # project no longer has is the stale-bundle bug, persisted.
        expect_length(store_read_autosave()$bundles, 0L)
        expect_identical(
          store_read_autosave()$experiments$proteomics$source_fingerprint,
          "fp-B")
      }
    )
  })
})

test_that("a restored autosave still knows which file it came from", {
  skip_if_not_installed("qs2")
  fi_store({
    proj <- fi_project(fi_input("fp-A", source_path = "raw/cheek__ab12.xlsx"))
    expect_true(store_autosave(proj))
    back <- store_read_autosave()
    # Both are needed after a restore: the fingerprint to recognise a
    # re-upload of the same file, the path to export a runnable script.
    expect_identical(back$experiments$proteomics$source_fingerprint, "fp-A")
    expect_identical(back$experiments$proteomics$source_path,
                     "raw/cheek__ab12.xlsx")
  })
})

test_that("an explicit save and the autosave do not tread on each other", {
  skip_if_not_installed("qs2")
  fi_store({
    named <- fi_project(fi_input("fp-named"))
    rolling <- fi_project(fi_input("fp-rolling"))
    expect_true(store_save_project(named, slug = "my-work")$ok)
    expect_true(store_autosave(rolling))
    expect_identical(
      store_load_project("my-work")$project$experiments$proteomics$source_fingerprint,
      "fp-named")
    expect_identical(
      store_read_autosave()$experiments$proteomics$source_fingerprint,
      "fp-rolling")
    expect_equal(list_saved_projects()$slug, "my-work")
  })
})

# ---- the import gate, end to end --------------------------------------

test_that("re-picking the same file leaves the analyses in place", {
  proj <- fi_project(fi_input("fp-A"), bundles = list(diff = "OLD"))
  shiny::testServer(
    import_view_server,
    args = list(current_project = shiny::reactiveVal(proj)),
    {
      parsed(list(input = fi_input("fp-A"), report = NULL))
      session$setInputs(confirm = 1)
      expect_null(confirmed_input())
    }
  )
  # Nothing was committed, so app_server never bumps the generation and
  # the views keep what they computed.
  expect_length(proj$bundles, 1L)
})

test_that("the confirm gate reads the layer matching the incoming omics type", {
  # A project holding proteomics receives an RNA-seq file: that extends
  # the project rather than replacing anything, so no dialog and no
  # clearing.
  proj <- fi_project(fi_input("fp-A"), bundles = list(diff = "OLD"))
  rnaseq <- fi_input("fp-rna")
  rnaseq$omics_type <- "rnaseq"
  shiny::testServer(
    import_view_server,
    args = list(current_project = shiny::reactiveVal(proj)),
    {
      parsed(list(input = rnaseq, report = NULL))
      session$setInputs(confirm = 1)
      expect_false(is.null(confirmed_input()))
      expect_identical(confirmed_input()$omics_type, "rnaseq")
    }
  )
})

# ---- the export, end to end -------------------------------------------

test_that("a project assembled by the app exports a script that names its file", {
  skip_if_not_installed("qs2")
  fi_store({
    src <- tempfile(fileext = ".csv"); on.exit(unlink(src), add = TRUE)
    writeLines("a,b\n1,2\n", src)
    archived <- store_raw_upload(src, "cheek.csv", "abcdef123456:proteomics:x")
    expect_true(archived$ok)

    inp <- fi_input("abcdef123456:proteomics:x", source_path = archived$path)
    proj <- fi_project(inp)
    lines <- omicsCore::export_script(proj)

    # The archive, the fingerprint and the script are three features that
    # only pay off together: the script is runnable because the file it
    # names is one the app kept.
    expect_true(any(grepl(basename(archived$path), lines, fixed = TRUE)))
    expect_false(any(grepl("was not archived", lines, fixed = TRUE)))
    expect_no_error(parse(text = paste(lines, collapse = "\n")))
  })
})

test_that("a demo project exports a script that says it has no data file", {
  proj <- fi_project(fi_input())
  lines <- omicsCore::export_script(proj)
  # Honest about the gap rather than emitting a path that does not exist.
  expect_true(any(grepl("was not archived", lines, fixed = TRUE)))
  expect_true(any(grepl("<path-to-", lines, fixed = TRUE)))
})
