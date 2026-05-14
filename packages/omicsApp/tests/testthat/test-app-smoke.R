# shinytest2 smoke harness for slice 2F. Boots the installed
# `omicsApp` app, cycles through the 7 sidebar nav items, asserts
# the page-title text for each view, and captures a portable
# input-register snapshot once on the Report view as a Phase-2 lock.
#
# Gates: shinytest2 + chromote are Suggests-only, and CRAN can't
# install Chromium. Skip when either is missing.
#
# What this DOES catch:
#   - silent UI breakage on view-switch (missing title, dropped
#     widget, accidental view rename)
#   - the input register changing shape (renamed `ns()`, dropped
#     `actionButton`, etc.)
#
# What this DOES NOT catch (deferred to Phase 4):
#   - plot pixel changes (output snapshots are ~16 MB of base64
#     PNG and false-positive across font stacks)
#   - input-driven flow regression (sliders -> volcano point
#     count, etc.)

test_that("omicsApp boots and all 7 views render", {
  skip_on_cran()
  skip_if_not_installed("shinytest2")
  skip_if_not_installed("chromote")

  app_dir <- system.file("app", package = "omicsApp")
  if (!nzchar(app_dir)) {
    skip("omicsApp is not installed (system.file('app') is empty).")
  }

  app <- shinytest2::AppDriver$new(
    app_dir,
    name         = "omicsApp-smoke",
    load_timeout = 20000
  )
  on.exit(app$stop(), add = TRUE)

  expected <- list(
    nav_project     = "Project overview",
    nav_import      = "Import data",
    nav_qc          = "Quality control",
    nav_diff        = "Differential",
    nav_enrich      = "Pathway enrichment",
    nav_integration = "Multi-omics integration",
    nav_report      = "Report"
  )

  for (id in names(expected)) {
    app$click(selector = paste0("#", id))
    app$wait_for_idle(timeout = 1000)
    # The hidden tabsetPanel mounts all 7 view bodies at once, so we
    # must scope to `.tab-pane.active` to read only the visible title.
    title <- app$get_text(".tab-pane.active .page-title")
    expect_match(title, expected[[id]], fixed = TRUE,
                 info = paste("view:", id))
  }

  # Lock the input register on the Report view. Inputs-only because
  # output snapshots include base64-encoded plot bytes and aren't
  # portable across font stacks.
  app$expect_values(output = FALSE, input = TRUE, export = FALSE)
})
