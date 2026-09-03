# The QC view's two controls used to be inert on the demo project: the
# view handed back a bundle computed at fixed settings, while the
# slider and the radio stayed enabled. A control that is enabled and
# does nothing reads as a broken app, which is how it was reported.

test_that("the outlier radio reaches run_qc on the demo project", {
  shiny::testServer(qc_view_server, args = list(), {
    session$setInputs(missing_threshold = 0.5, outlier_method = "iqr")
    expect_identical(last_bundle()$params$outlier_method, "iqr")

    session$setInputs(outlier_method = "pca")
    expect_identical(last_bundle()$params$outlier_method, "pca")

    session$setInputs(outlier_method = "connectivity")
    expect_identical(last_bundle()$params$outlier_method, "connectivity")
  })
})

test_that("the missing-rate slider changes what the demo keeps", {
  # Not just that the parameter is passed through, but that it has an
  # effect: the demo carries ~5% NA, so a 2% cutoff has to drop
  # features a 50% cutoff keeps.
  shiny::testServer(qc_view_server, args = list(), {
    session$setInputs(missing_threshold = 0.5, outlier_method = "iqr")
    kept_loose <- nrow(last_bundle()$results$cleaned_input$expr_mat)

    session$setInputs(missing_threshold = 0.02)
    kept_tight <- nrow(last_bundle()$results$cleaned_input$expr_mat)

    expect_lt(kept_tight, kept_loose)
  })
})

test_that("the demo input is stable across calls", {
  # It is cached and seeded, so two reads must give the same NA
  # pattern -- otherwise the panel would shift under the user every
  # time a control moved.
  a <- example_qc_input()
  b <- example_qc_input()
  expect_identical(which(is.na(a$expr_mat)), which(is.na(b$expr_mat)))
  expect_true(any(is.na(a$expr_mat)))
})

test_that("example_qc_bundle still answers at its documented settings", {
  # Other callers (the Report view) take the bundle rather than the
  # input, and expect the fixed 0.5 / iqr settings.
  b <- example_qc_bundle()
  expect_identical(b$params$outlier_method, "iqr")
  expect_equal(b$params$missing_threshold, 0.5)
})
