# Deleting a project used to remove a few MB of project file and leave
# the 170 MB workbook behind it on disk -- so the one thing a user does
# to free space freed none of it. And a project was all-or-nothing, so a
# mistaken RNA-seq layer cost the proteomics work sitting next to it.

del_input <- function(fingerprint = NULL, ids = c("S1", "S2")) {
  m <- matrix(seq_len(4 * length(ids)) * 1000, nrow = 4,
              dimnames = list(paste0("F", 1:4), ids))
  inp <- omicsCore::omics_input(
    m,
    data.frame(sample_id = ids, condition = "G1", row.names = ids,
               stringsAsFactors = FALSE),
    data.frame(feature_id = rownames(m), row.names = rownames(m),
               stringsAsFactors = FALSE),
    omics_type = "proteomics", assay_type = "raw_intensity")
  if (!is.null(fingerprint)) inp$source_fingerprint <- fingerprint
  inp
}

local_store <- function(env = parent.frame()) {
  dir <- withr::local_tempdir(.local_envir = env)
  withr::local_envvar(c(OMICSAPP_DATA_DIR = dir), .local_envir = env)
  dir
}

fake_upload <- function(bytes = 500000) {
  p <- withr::local_tempfile(fileext = ".xlsx", .local_envir = parent.frame())
  writeBin(as.raw(rep(0L, bytes)), p)
  p
}

test_that("deleting a project takes its archived upload with it", {
  local_store()
  up <- fake_upload()

  store_save_project(
    omicsCore::omics_project("A", list(x = del_input("aaaaaaaaaaaa:1"))), "a")
  store_raw_upload(up, "big.xlsx", "aaaaaaaaaaaa:1")
  expect_length(list.files(raw_dir()), 1L)

  res <- store_delete_project("a")
  expect_true(res$ok)
  expect_length(list.files(raw_dir()), 0L)
  # The freed size is the point of the message: it is what tells someone
  # whether deleting actually gave them room back.
  expect_match(res$message, "freed")
})

test_that("an upload another project still uses is kept", {
  local_store()
  up <- fake_upload()

  for (s in c("b", "c")) {
    store_save_project(
      omicsCore::omics_project(s, list(x = del_input("bbbbbbbbbbbb:1"))), s)
  }
  store_raw_upload(up, "shared.xlsx", "bbbbbbbbbbbb:1")

  res <- store_delete_project("b")
  expect_true(res$ok)
  expect_length(list.files(raw_dir()), 1L)
  expect_no_match(res$message, "freed")

  # Now nothing refers to it.
  store_delete_project("c")
  expect_length(list.files(raw_dir()), 0L)
})

test_that("a project with no archived upload deletes quietly", {
  local_store()
  store_save_project(omicsCore::omics_project("D", list(x = del_input())), "d")
  res <- store_delete_project("d")
  expect_true(res$ok)
  expect_identical(res$message, "Deleted 'd'.")
})

test_that("a layer can be removed without touching the others", {
  proj <- shiny::reactiveVal(omicsCore::omics_project(
    "p", list(prot = del_input(), rna = del_input())))

  shiny::testServer(project_view_server,
                    args = list(current_project = proj), {
    session$flushReact()
    session$setInputs(drop_layer_2 = 1)
    session$flushReact()
    # Confirmed, not immediate: View and Remove are one word apart, and
    # Remove throws away every result computed on that layer.
    expect_identical(pending_drop(), "rna")
    expect_length(current_project()$experiments, 2L)

    session$setInputs(confirm_drop_layer = 1)
    session$flushReact()
    expect_identical(names(current_project()$experiments), "prot")
  })
})

test_that("removing a layer drops the sample links that named it", {
  # A link to a layer that is gone would pair samples to nothing.
  p <- omicsCore::omics_project("p",
                                list(prot = del_input(), rna = del_input()))
  p$sample_link <- data.frame(
    tag = c("prot", "prot", "rna", "rna"),
    sample_id = c("S1", "S2", "S1", "S2"),
    donor_id = c("D1", "D2", "D1", "D2"),
    stringsAsFactors = FALSE)
  proj <- shiny::reactiveVal(p)

  shiny::testServer(project_view_server,
                    args = list(current_project = proj), {
    session$flushReact()
    session$setInputs(drop_layer_2 = 1)
    session$flushReact()
    session$setInputs(confirm_drop_layer = 1)
    session$flushReact()

    link <- current_project()$sample_link
    expect_identical(unique(link$tag), "prot")
    expect_equal(nrow(link), 2L)
  })
})

test_that("cancelling leaves the project alone", {
  proj <- shiny::reactiveVal(omicsCore::omics_project(
    "p", list(prot = del_input(), rna = del_input())))

  shiny::testServer(project_view_server,
                    args = list(current_project = proj), {
    session$flushReact()
    session$setInputs(drop_layer_1 = 1)
    session$flushReact()
    expect_identical(pending_drop(), "prot")
    # No confirm click: the modal was dismissed.
    expect_length(current_project()$experiments, 2L)
  })
})
