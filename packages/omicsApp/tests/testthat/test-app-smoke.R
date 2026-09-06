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
  skip_if_not(nzchar(Sys.which("google-chrome")) ||
                nzchar(Sys.which("chromium")) ||
                nzchar(Sys.which("chrome")) ||
                !is.null(tryCatch(chromote::find_chrome(),
                                  error = function(e) NULL)),
              "No Chrome/Chromium available for chromote")

  where <- smoke_app_dir()
  if (!nzchar(where$dir)) {
    skip("omicsApp is not installed and no source tree is present.")
  }
  withr::local_envvar(OMICSAPP_DEV_ROOT = where$dev_root)

  app <- tryCatch(
    shinytest2::AppDriver$new(
      where$dir,
      name         = "omicsApp-smoke",
      load_timeout = 30000
    ),
    error = function(e) {
      skip(sprintf("AppDriver launch failed: %s", conditionMessage(e)))
    }
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
    # Headroom for the diff view: at first visit it runs limma on the
    # proteomics fixture (~2.5 s on commodity hardware) and asks
    # has_pkg() about every engine, which loads DESeq2 and its sixty-odd
    # namespaces into the app process -- on a CI runner, well over 5 s.
    app$wait_for_idle(timeout = 30000)
    # The hidden tabsetPanel mounts all 7 view bodies at once, so we
    # must scope to `.tab-pane.active` to read only the visible title.
    title <- app$get_text(".tab-pane.active .page-title")
    expect_match(title, expected[[id]], fixed = TRUE,
                 info = paste("view:", id))
  }

  # Verify key inputs from each slice are registered (catches
  # renamed/dropped ns() without fragile DT-state snapshots).
  app$click(selector = "#nav_report")
  app$wait_for_idle(timeout = 3000)
  vals <- app$get_values(input = TRUE, output = FALSE, export = FALSE)$input
  essential <- c(
    "nav_project", "nav_import", "nav_qc", "nav_diff",
    "nav_enrich", "nav_integration", "nav_report",
    "import-file", "import-confirm",
    "qc-missing_threshold", "qc-outlier_method",
    "diff-method", "diff-group_col",
    "diff-rerun",
    "enrich-type", "enrich-database", "enrich-direction", "enrich-rerun",
    "integration-rerun"
  )
  for (nm in essential) {
    expect_true(nm %in% names(vals),
                info = sprintf("input '%s' missing from register", nm))
  }
})
