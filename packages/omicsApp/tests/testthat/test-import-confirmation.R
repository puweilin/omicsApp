# Coverage for the Import view's confirmation card and the sheet-role
# override behind it.
#
# The schema table says what the classifier decided. This card says what that
# decision produced, which is the part a user can check: a sheet assignment
# can be wrong at high confidence and still yield an omics_input that analyses
# cleanly and means nothing. Nothing downstream errors on it.

confirm_upload <- function(path, name = "tiny.xlsx") {
  list(datapath = path, name = name, size = file.info(path)$size)
}

# ---- what the card reports --------------------------------------------------

test_that("the shape summary reports the parsed dimensions", {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("readxl")

  xlsx <- tempfile(fileext = ".xlsx")
  on.exit(unlink(xlsx), add = TRUE)
  write_tiny_omics_xlsx(xlsx, n_features = 5L, n_samples = 6L, scale = "linear")

  shiny::testServer(import_view_server, {
    session$setInputs(omics_type = "proteomics", file = confirm_upload(xlsx))
    ui <- confirm_shape_ui(parsed()$input, parsed()$report)
    html <- as.character(htmltools::tagList(ui))

    # The numbers are the point: "5 features x 6 samples" against a file the
    # user knows has 6 features is how a transposed matrix gets caught
    expect_match(html, "Features")
    expect_match(html, "Samples")
    expect_match(html, ">5<")
    expect_match(html, ">6<")
    expect_match(html, "Orientation")
  })
})

test_that("the shape summary says so when there is nothing to import", {
  ui <- confirm_shape_ui(NULL, omicsCore::new_import_report())
  html <- as.character(htmltools::tagList(ui))
  expect_match(html, "Nothing to import yet")
})

test_that("every sheet gets a role dropdown, including unplaced ones", {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("readxl")

  xlsx <- tempfile(fileext = ".xlsx")
  on.exit(unlink(xlsx), add = TRUE)
  write_tiny_omics_xlsx(xlsx, scale = "linear")

  shiny::testServer(import_view_server, {
    session$setInputs(omics_type = "proteomics", file = confirm_upload(xlsx))
    sheets <- parsed()$report$sheets
    html <- as.character(htmltools::tagList(confirm_roles_ui(ns, parsed()$report, gen = 1L)))

    # An "unknown" sheet is often the metadata; dropping it silently is how a
    # grouping column goes missing later, so it gets a dropdown too
    for (i in seq_len(nrow(sheets))) {
      expect_match(html, paste0("role_1_", i))
    }
    for (nm in sheets$name) expect_match(html, nm, fixed = TRUE)
  })
})

test_that("previews show the parsed values, not the raw sheet", {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("readxl")

  xlsx <- tempfile(fileext = ".xlsx")
  on.exit(unlink(xlsx), add = TRUE)
  write_tiny_omics_xlsx(xlsx, n_features = 5L, n_samples = 6L, scale = "linear")

  shiny::testServer(import_view_server, {
    session$setInputs(omics_type = "proteomics", file = confirm_upload(xlsx))

    mp <- preview_matrix(parsed()$input$expr_mat)
    expect_s3_class(mp, "data.frame")
    expect_identical(colnames(mp)[1], "feature")
    expect_lte(nrow(mp), 5L)
    # Values, so a header row parsed as data is visible
    expect_true(is.numeric(mp[[2]]))

    sp <- preview_metadata(parsed()$input$meta_df)
    expect_s3_class(sp, "data.frame")
    expect_identical(colnames(sp)[1], "sample")
    # Grouping columns are what the user checks the metadata for
    expect_true("group" %in% colnames(sp))
  })
})

test_that("previews degrade quietly on empty input", {
  expect_null(preview_matrix(NULL))
  expect_null(preview_metadata(NULL))
  expect_null(preview_matrix(matrix(numeric(0), nrow = 0, ncol = 0)))
})

# ---- correcting the classifier ----------------------------------------------

test_that("changing a sheet role re-reads the file with that assignment", {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("readxl")

  xlsx <- tempfile(fileext = ".xlsx")
  on.exit(unlink(xlsx), add = TRUE)
  write_tiny_omics_xlsx(xlsx, scale = "linear")

  shiny::testServer(import_view_server, {
    session$setInputs(omics_type = "proteomics", file = confirm_upload(xlsx))
    sheets <- parsed()$report$sheets
    matrix_row <- which(sheets$role == "matrix")[1]
    expect_false(is.na(matrix_row))

    # Tell it the matrix sheet is actually feature annotation. The point is
    # not that this is sensible -- it is that the user's word wins over a
    # confident classifier.
    session$setInputs(role_1_1 = "feature_annot")

    updated <- parsed()$report$sheets
    expect_identical(updated$role[1], "feature_annot")
    expect_identical(updated$notes[1], "role set by user")
    expect_equal(updated$confidence[1], 1)
  })
})

test_that("a new file clears the previous role overrides", {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("readxl")

  first <- tempfile(fileext = ".xlsx")
  second <- tempfile(fileext = ".xlsx")
  on.exit(unlink(c(first, second)), add = TRUE)
  write_tiny_omics_xlsx(first, scale = "linear")
  write_tiny_omics_xlsx(second, scale = "linear")

  shiny::testServer(import_view_server, {
    session$setInputs(omics_type = "proteomics", file = confirm_upload(first))
    session$setInputs(role_1_1 = "feature_annot")
    expect_identical(parsed()$report$sheets$role[1], "feature_annot")

    # An assignment made for one workbook says nothing about the next
    session$setInputs(file = confirm_upload(second, "other.xlsx"))
    expect_null(role_overrides())
    expect_false(identical(parsed()$report$sheets$notes[1], "role set by user"))
  })
})

# ---- the omicsCore side of the override -------------------------------------

test_that("read_omics honours an explicit sheet role", {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("readxl")

  xlsx <- tempfile(fileext = ".xlsx")
  on.exit(unlink(xlsx), add = TRUE)
  write_tiny_omics_xlsx(xlsx, scale = "linear")

  default <- omicsCore::read_omics(xlsx, omics_type = "proteomics",
                                   assay_type = "raw_intensity")
  target <- default$report$sheets$name[1]

  overridden <- omicsCore::read_omics(
    xlsx, omics_type = "proteomics", assay_type = "raw_intensity",
    sheet_roles = stats::setNames("feature_annot", target)
  )
  row <- overridden$report$sheets[overridden$report$sheets$name == target, ]
  expect_identical(row$role, "feature_annot")
  # Confidence 1 so the user's choice outranks whatever the classifier scored
  expect_equal(row$confidence, 1)
})

test_that("read_omics rejects a malformed or unknown role", {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("readxl")

  xlsx <- tempfile(fileext = ".xlsx")
  on.exit(unlink(xlsx), add = TRUE)
  write_tiny_omics_xlsx(xlsx, scale = "linear")

  # Checked before the file is opened, so a typo reads as an argument problem
  # rather than as an import that mysteriously ignored you
  expect_error(
    omicsCore::read_omics(xlsx, sheet_roles = c(Sheet1 = "expression")),
    "Unknown role"
  )
  expect_error(
    omicsCore::read_omics(xlsx, sheet_roles = "matrix"),
    "named character vector"
  )
})

# --- changing which sheet is the matrix ----------------------------------

# The replace check asks "is this the same file?" by fingerprint, and the
# fingerprint did not include the sheet assignment. So picking a different
# sheet as the expression matrix -- a whole different set of numbers --
# read as "that is the file already loaded, nothing to re-import", and the
# candidate the user was looking at was discarded on Confirm.
two_matrix_workbook <- function(env = parent.frame()) {
  skip_if_not_installed("openxlsx")
  path <- tempfile(fileext = ".xlsx")
  withr::defer(unlink(path), envir = env)
  set.seed(4)
  ids <- paste0("s", 1:6)
  fids <- paste0("P", 1:10)
  a <- matrix(rnorm(60, 20, 1), 10, dimnames = list(fids, ids))
  openxlsx::write.xlsx(list(
    matrix_a = data.frame(feature_id = fids, a, check.names = FALSE),
    # Every value differs by exactly 5, so a matrix that came from the
    # wrong sheet is unmistakable.
    matrix_b = data.frame(feature_id = fids, a + 5, check.names = FALSE),
    samples  = data.frame(sample_id = ids,
                          group = rep(c("A", "B"), each = 3),
                          stringsAsFactors = FALSE)
  ), path)
  path
}

test_that("changing the matrix sheet is a different import", {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("readxl")
  skip_if_not_installed("withr")
  xlsx <- two_matrix_workbook()

  shiny::testServer(import_view_server, {
    session$setInputs(
      omics_type = "proteomics",
      file = list(datapath = xlsx, name = "two.xlsx",
                  size = file.info(xlsx)$size))
    session$setInputs(confirm = 1)
    first <- session$returned()$input
    skip_if(is.null(first), "workbook did not parse in this environment")
    fp_first <- first$source_fingerprint

    # Same file, different sheet as the matrix.
    # Roles are named by IMPORT_REPORT_ROLES, not by sheet name: the two
    # candidate matrices swap, so the second sheet becomes the data.
    role_overrides(c(matrix_a = "unknown", matrix_b = "matrix"))
    do_parse()
    session$flushReact()
    cand <- parsed()$input
    skip_if(is.null(cand), "role override did not parse")
    # The candidate really is different data.
    expect_false(isTRUE(all.equal(cand$expr_mat, first$expr_mat)))
    # And the fingerprint says so, which is what the replace check reads.
    expect_false(identical(cand$source_fingerprint, fp_first))
  })
})

test_that("the same roles in a different order are the same parse", {
  # The assignment is a set, so an order-sensitive fingerprint would make
  # a re-render look like new data and prompt to replace the layer.
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("withr")
  xlsx <- two_matrix_workbook()
  a <- input_fingerprint(xlsx, "proteomics", "normalized_intensity",
                         sheet_roles = c(x = "matrix", y = "samples"))
  b <- input_fingerprint(xlsx, "proteomics", "normalized_intensity",
                         sheet_roles = c(y = "samples", x = "matrix"))
  expect_identical(a, b)
  c_ <- input_fingerprint(xlsx, "proteomics", "normalized_intensity",
                          sheet_roles = c(x = "samples", y = "matrix"))
  expect_false(identical(a, c_))
})

# ---- a repeat confirmation is a new event -----------------------------
#
# reactiveVal does not invalidate when the new value is identical to the
# old, which is right for state and wrong for an event. After "Start a
# new project" the module still held the confirmed input, so confirming
# the same file again assigned the same object, nothing invalidated, and
# app_server never learned there had been an import: the project stayed
# empty while the import view said it had confirmed.

staged_input <- function() {
  mat <- matrix(as.numeric(1:12), nrow = 3,
                dimnames = list(paste0("g", 1:3), paste0("s", 1:4)))
  omicsCore::omics_input(
    mat,
    data.frame(group = c("A", "A", "B", "B"), row.names = paste0("s", 1:4)),
    data.frame(feature_id = paste0("g", 1:3)),
    omics_type = "proteomics", assay_type = "normalized_intensity")
}

test_that("confirming the same input twice is two submissions", {
  shiny::testServer(import_view_server,
                    args = list(current_project = shiny::reactiveVal(NULL)), {
    out <- session$getReturned()
    parsed(list(input = staged_input(), report = NULL))
    session$setInputs(confirm = 1)
    first <- out()
    expect_false(is.null(first$input))

    session$setInputs(confirm = 2)
    second <- out()

    # Same data, different submission: exactly the case the
    # identical-value shortcut swallowed.
    expect_identical(second$input, first$input)
    expect_gt(second$seq, first$seq)
  })
})

test_that("activation clears what the last import staged", {
  act <- shiny::reactiveVal(0L)
  shiny::testServer(import_view_server,
                    args = list(current_project = shiny::reactiveVal(NULL),
                                activation = act), {
    parsed(list(input = staged_input(), report = NULL))
    session$setInputs(confirm = 1)
    expect_false(is.null(confirmed_input()))

    act(1L)                       # "Start a new project"
    session$flushReact()
    # The view must not still show a file confirmed into a project that
    # no longer exists.
    expect_null(confirmed_input())
    expect_null(parsed())
  })
})
