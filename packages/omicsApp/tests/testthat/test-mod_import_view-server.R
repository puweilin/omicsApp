# testServer harness for the Import view (slice 3A).
#
# Drives the module with a synthetic xlsx workbook produced by
# `write_tiny_omics_xlsx()` (in helper-tiny-xlsx.R) and asserts that:
#   * uploading a parseable file populates the `parsed` reactive
#     with a valid `omics_input`;
#   * clicking Confirm exposes that input as the module's return
#     value (the contract slice 3B will consume);
#   * a missing file path produces a warning instead of crashing.
#
# Note on shinytest2: the slice 2F smoke harness covers boot + nav
# but not parameter-driven flows. testServer() runs purely in R, so
# it's the right tool for slice 3A.

test_that("import view exposes an omics_input after Confirm", {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("readxl")

  xlsx <- tempfile(fileext = ".xlsx")
  on.exit(unlink(xlsx), add = TRUE)
  write_tiny_omics_xlsx(xlsx)

  shiny::testServer(import_view_server, {
    expect_null(session$returned())

    session$setInputs(
      omics_type = "proteomics",
      file       = list(datapath = xlsx,
                        name     = "tiny.xlsx",
                        size     = file.info(xlsx)$size)
    )
    # Parsing should have succeeded.
    expect_true(parse_ok())
    # But not yet confirmed.
    expect_null(session$returned())

    session$setInputs(confirm = 1)
    inp <- session$returned()
    expect_s3_class(inp, "omics_input")
    expect_equal(nrow(inp$expr_mat), 5L)
    expect_equal(ncol(inp$expr_mat), 6L)
    expect_equal(inp$omics_type, "proteomics")
  })
})

test_that("import view surfaces read_omics warnings without crashing", {
  shiny::testServer(import_view_server, {
    session$setInputs(
      omics_type = "proteomics",
      file       = list(datapath = "/no/such/file.xlsx",
                        name     = "broken.xlsx",
                        size     = 0)
    )
    # parsed() is populated (with the synthetic error-report); parse_ok
    # is FALSE so Confirm stays inert.
    expect_false(parse_ok())
    rep <- parsed()$report
    expect_s3_class(rep, "ImportReport")
    expect_true(any(grepl("does not exist|failed", rep$warnings)))

    # Pressing Confirm anyway must not flip the returned reactive.
    session$setInputs(confirm = 1)
    expect_null(session$returned())
  })
})

test_that("changing omics_type after upload re-parses and resets confirm", {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("readxl")

  xlsx <- tempfile(fileext = ".xlsx")
  on.exit(unlink(xlsx), add = TRUE)
  write_tiny_omics_xlsx(xlsx)

  shiny::testServer(import_view_server, {
    session$setInputs(
      omics_type = "proteomics",
      file       = list(datapath = xlsx,
                        name     = "tiny.xlsx",
                        size     = file.info(xlsx)$size)
    )
    session$setInputs(confirm = 1)
    expect_s3_class(session$returned(), "omics_input")
    expect_equal(session$returned()$omics_type, "proteomics")

    # Flip the radio. Should reset the confirmed input and rebuild the
    # candidate with the new omics_type.
    session$setInputs(omics_type = "rnaseq")
    expect_null(session$returned())
    expect_true(parse_ok())
    expect_equal(parsed()$input$omics_type, "rnaseq")

    session$setInputs(confirm = 2)
    expect_equal(session$returned()$omics_type, "rnaseq")
  })
})
