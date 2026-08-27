# testServer coverage for the scale declaration and normalization step the
# Import view gained in Phase B.
#
# The point of these tests is that the layer reaching the project is on a scale
# the differential backends can use. limma applies no transform of its own, so
# importing linear intensities unchanged means limma runs on raw instrument
# output and still returns a full, plausible, wrong result table. Nothing
# downstream errors, which is why the check has to live here.

skip_unless_xlsx <- function() {
  testthat::skip_if_not_installed("openxlsx")
  testthat::skip_if_not_installed("readxl")
}

tiny_upload <- function(path, name = "tiny.xlsx") {
  list(datapath = path, name = name, size = file.info(path)$size)
}

# ---- scale inference --------------------------------------------------------

test_that("a linear proteomics upload is labelled raw_intensity", {
  skip_unless_xlsx()
  xlsx <- tempfile(fileext = ".xlsx")
  on.exit(unlink(xlsx), add = TRUE)
  write_tiny_omics_xlsx(xlsx, scale = "linear")

  shiny::testServer(import_view_server, {
    session$setInputs(omics_type = "proteomics", file = tiny_upload(xlsx))
    expect_true(parse_ok())
    expect_identical(parsed()$input$assay_type, "raw_intensity")
    expect_true(normalizable())
  })
})

test_that("an already-log2 proteomics upload is labelled normalized_intensity", {
  skip_unless_xlsx()
  xlsx <- tempfile(fileext = ".xlsx")
  on.exit(unlink(xlsx), add = TRUE)
  write_tiny_omics_xlsx(xlsx, scale = "log2")

  shiny::testServer(import_view_server, {
    session$setInputs(omics_type = "proteomics", file = tiny_upload(xlsx))
    expect_identical(parsed()$input$assay_type, "normalized_intensity")
    # Nothing to normalize, so the step is not offered
    expect_false(normalizable())
  })
})

test_that("RNA-seq is labelled raw_count and never offered normalization", {
  skip_unless_xlsx()
  xlsx <- tempfile(fileext = ".xlsx")
  on.exit(unlink(xlsx), add = TRUE)
  write_tiny_omics_xlsx(xlsx, scale = "linear")

  shiny::testServer(import_view_server, {
    session$setInputs(omics_type = "rnaseq", file = tiny_upload(xlsx))
    expect_identical(parsed()$input$assay_type, "raw_count")
    # DESeq2 and edgeR model counts directly; the generic backends
    # log-transform raw_count themselves
    expect_false(normalizable())
  })
})

# ---- user override ----------------------------------------------------------

test_that("choosing a different scale rewrites the label and drops confirmation", {
  skip_unless_xlsx()
  xlsx <- tempfile(fileext = ".xlsx")
  on.exit(unlink(xlsx), add = TRUE)
  write_tiny_omics_xlsx(xlsx, scale = "linear")

  shiny::testServer(import_view_server, {
    session$setInputs(omics_type = "proteomics", file = tiny_upload(xlsx))
    session$setInputs(confirm = 1)
    expect_s3_class(session$returned(), "omics_input")

    # The inference is a starting point; the user gets the last word
    session$setInputs(assay_type = "imputed_intensity")
    expect_identical(parsed()$input$assay_type, "imputed_intensity")
    # Relabelling changes what the analyses will do, so the previous
    # confirmation no longer stands
    expect_null(session$returned())
    expect_false(normalizable())
  })
})

# ---- normalization on commit ------------------------------------------------

test_that("confirming a linear upload normalizes it by default", {
  skip_unless_xlsx()
  xlsx <- tempfile(fileext = ".xlsx")
  on.exit(unlink(xlsx), add = TRUE)
  write_tiny_omics_xlsx(xlsx, scale = "linear")

  shiny::testServer(import_view_server, {
    session$setInputs(omics_type = "proteomics", file = tiny_upload(xlsx))
    # normalize defaults to TRUE; method defaults to vsn, but log2 keeps the
    # test independent of whether vsn is installed
    session$setInputs(normalize = TRUE, normalize_method = "log2")
    session$setInputs(confirm = 1)

    inp <- session$returned()
    expect_s3_class(inp, "omics_input")
    expect_identical(inp$assay_type, "normalized_intensity")
    # The layer that lands in the project is now on a scale limma can use.
    # check_assay_scale() is the semantic form of "values match the label".
    expect_silent(omicsCore::check_assay_scale(inp))
    expect_lt(max(inp$expr_mat, na.rm = TRUE), 100)
  })
})

test_that("normalization keeps the pre-normalization matrix", {
  skip_unless_xlsx()
  xlsx <- tempfile(fileext = ".xlsx")
  on.exit(unlink(xlsx), add = TRUE)
  write_tiny_omics_xlsx(xlsx, scale = "linear")

  shiny::testServer(import_view_server, {
    session$setInputs(omics_type = "proteomics", file = tiny_upload(xlsx))
    session$setInputs(normalize = TRUE, normalize_method = "log2")
    session$setInputs(confirm = 1)

    inp <- session$returned()
    expect_false(is.null(inp$raw_mat))
    expect_identical(dim(inp$raw_mat), dim(inp$expr_mat))
    # raw_mat holds the linear values, expr_mat the transformed ones
    expect_gt(max(inp$raw_mat, na.rm = TRUE), max(inp$expr_mat, na.rm = TRUE))
  })
})

test_that("unticking normalize imports the linear values unchanged", {
  skip_unless_xlsx()
  xlsx <- tempfile(fileext = ".xlsx")
  on.exit(unlink(xlsx), add = TRUE)
  write_tiny_omics_xlsx(xlsx, scale = "linear")

  shiny::testServer(import_view_server, {
    session$setInputs(omics_type = "proteomics", file = tiny_upload(xlsx))
    session$setInputs(normalize = FALSE)
    session$setInputs(confirm = 1)

    inp <- session$returned()
    expect_identical(inp$assay_type, "raw_intensity")
    expect_gt(max(inp$expr_mat, na.rm = TRUE), 100)
  })
})

test_that("an already-log2 upload is committed untouched", {
  skip_unless_xlsx()
  xlsx <- tempfile(fileext = ".xlsx")
  on.exit(unlink(xlsx), add = TRUE)
  write_tiny_omics_xlsx(xlsx, scale = "log2")

  shiny::testServer(import_view_server, {
    session$setInputs(omics_type = "proteomics", file = tiny_upload(xlsx))
    session$setInputs(confirm = 1)

    inp <- session$returned()
    expect_identical(inp$assay_type, "normalized_intensity")
    # No second transform: normalizing twice would compress the dynamic range
    expect_null(inp$raw_mat)
  })
})

test_that("RNA-seq counts are committed without transformation", {
  skip_unless_xlsx()
  xlsx <- tempfile(fileext = ".xlsx")
  on.exit(unlink(xlsx), add = TRUE)
  write_tiny_omics_xlsx(xlsx, scale = "linear")

  shiny::testServer(import_view_server, {
    session$setInputs(omics_type = "rnaseq", file = tiny_upload(xlsx))
    session$setInputs(confirm = 1)

    inp <- session$returned()
    expect_identical(inp$assay_type, "raw_count")
    expect_null(inp$raw_mat)
  })
})

test_that("vsn is used when asked for and available", {
  skip_unless_xlsx()
  testthat::skip_if_not_installed("vsn")

  xlsx <- tempfile(fileext = ".xlsx")
  on.exit(unlink(xlsx), add = TRUE)
  # vsn needs more rows than its minDataPointsPerStratum default of 42
  write_tiny_omics_xlsx(xlsx, n_features = 60L, n_samples = 6L, scale = "linear")

  shiny::testServer(import_view_server, {
    session$setInputs(omics_type = "proteomics", file = tiny_upload(xlsx))
    session$setInputs(normalize = TRUE, normalize_method = "vsn")
    session$setInputs(confirm = 1)

    inp <- session$returned()
    expect_identical(inp$assay_type, "normalized_intensity")
    expect_identical(dim(inp$expr_mat), c(60L, 6L))
    expect_silent(omicsCore::check_assay_scale(inp))
  })
})

# ---- fingerprint ------------------------------------------------------------

test_that("the normalization choice is part of the upload's identity", {
  skip_unless_xlsx()
  xlsx <- tempfile(fileext = ".xlsx")
  on.exit(unlink(xlsx), add = TRUE)
  write_tiny_omics_xlsx(xlsx, scale = "linear")

  # Re-importing the same file with normalization off yields different
  # numbers, so it has to read as new data rather than "nothing changed"
  fp_on  <- input_fingerprint(xlsx, "proteomics", "raw_intensity", "vsn")
  fp_off <- input_fingerprint(xlsx, "proteomics", "raw_intensity", "none")
  expect_false(identical(fp_on, fp_off))

  # Same settings, same identity
  expect_identical(
    fp_on, input_fingerprint(xlsx, "proteomics", "raw_intensity", "vsn"))
})
