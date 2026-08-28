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
