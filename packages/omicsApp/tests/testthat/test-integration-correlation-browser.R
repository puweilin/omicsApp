# The correlation scatter in a real browser: a saved two-layer project
# with paired samples, opened from My projects, switched to sample-level
# correlation; a point hovered shows its card, and the figure downloads
# as a dated PNG.
#
# The saved test projects only ever had concordance results, so this
# figure's hover card and download had never been seen in a browser.
# testServer() checks the card's text (test-plot-hover.R) and the
# download's file (test-plot-download.R); only a browser checks that the
# pointer, over a drawn point, brings the card up there.
#
# Same gates as the journey test: shinytest2, chromote, a Chrome, and
# the source tree.

# Two layers measured on the same ten donors; a third of the genes share
# a donor-level signal between protein and RNA, so some correlate.
browser_corr_project <- function() {
  set.seed(11)
  n_genes <- 150L
  n <- 10L
  donor <- sprintf("D%02d", seq_len(n))
  grp <- rep(c("Control", "Treated"), each = n / 2L)
  sym <- sprintf("GENE%03d", seq_len(n_genes))
  shared <- matrix(stats::rnorm(n_genes * n), n_genes)
  load <- ifelse(seq_len(n_genes) <= 50L, 1.5, 0)
  layer <- function(prefix, omics) {
    m <- 10 + load * shared + matrix(stats::rnorm(n_genes * n, 0, 0.5), n_genes)
    ids <- paste0(prefix, seq_len(n_genes))
    samples <- paste0(prefix, "_", donor)
    dimnames(m) <- list(ids, samples)
    omicsCore::omics_input(
      m, data.frame(group = grp, donor = donor, row.names = samples),
      data.frame(feature_id = ids, feature_symbol = sym, row.names = ids),
      omics_type = omics,
      assay_type = if (omics == "rnaseq") "logcpm" else "normalized_intensity")
  }
  omicsCore::omics_project("Paired donors",
                           list(proteomics = layer("p", "proteomics"),
                                rnaseq = layer("r", "rnaseq")))
}

test_that("the correlation scatter's points show a card on hover, and it downloads", {
  skip_on_cran()
  skip_if_not_installed("shinytest2")
  skip_if_not_installed("chromote")
  skip_if_not_installed("withr")
  skip_if_not(!is.null(tryCatch(chromote::find_chrome(), error = function(e) NULL)),
              "No Chrome/Chromium available for chromote")
  where <- smoke_app_dir()
  skip_if(!nzchar(where$dir), "no app directory to launch")

  store <- file.path(tempfile("corr-browser-"), "store")
  dir.create(store, recursive = TRUE)
  on.exit(unlink(dirname(store), recursive = TRUE), add = TRUE)
  withr::local_envvar(OMICSAPP_DEV_ROOT = where$dev_root, OMICSAPP_DATA_DIR = store)
  proj <- browser_corr_project()
  slug <- project_slug(proj$name)
  expect_true(store_save_project(proj, slug, dir = store)$ok)

  # The point to hover: the highest on the plot, which stands apart. The
  # app runs the same correlation on the same project, so its figure is
  # this one.
  g <- omicsCore::plot_integration(omicsCore::run_integration(proj, "correlation"),
                                   view = "scatter")
  built <- ggplot2::ggplot_build(g)
  pts <- built$data[[which(vapply(g$layers, function(l) inherits(l$geom, "GeomPoint"),
                                  logical(1)))[1]]]
  top <- which.max(pts$y)
  target <- list(x = pts$x[[top]], y = pts$y[[top]],
                 name = as.character(g$data$feature_symbol[[top]]))

  app <- tryCatch(
    shinytest2::AppDriver$new(where$dir, name = "omicsApp-corr-browser",
                              width = 1366, height = 900, load_timeout = 30000,
                              timeout = 30000, seed = 1),
    error = function(e) skip(sprintf("AppDriver launch failed: %s", conditionMessage(e))))
  on.exit(app$stop(), add = TRUE)

  # Opening a project and the view's first visit take seconds on a busy
  # runner: wait for the server to settle rather than for one output.
  # Picking a project changes no output, so nothing to wait for there.
  app$set_inputs(`project-saved_pick` = slug, wait_ = FALSE)
  app$wait_for_idle(timeout = 30000)
  app$click("project-open_project", wait_ = FALSE)
  app$wait_for_idle(timeout = 30000)
  app$run_js("document.getElementById('nav_integration').click()")
  app$wait_for_idle(timeout = 30000)
  app$set_inputs(`integration-method` = "correlation", wait_ = FALSE)
  app$wait_for_value(output = "integration-cor_scatter", timeout = 30000)
  app$wait_for_idle(timeout = 30000)

  # Data to screen through the coordmap Shiny sent with the image (in
  # the image's CSS pixels), then a real pointer move over that spot.
  cm <- app$get_value(output = "integration-cor_scatter")$coordmap
  d <- cm$panels[[1]]$domain
  r <- cm$panels[[1]]$range
  b <- app$get_chromote_session()
  js <- function(x) b$Runtime$evaluate(x, returnByValue = TRUE)$result$value
  js("document.getElementById('integration-cor_scatter').scrollIntoView({block: 'center'})")
  img <- js(paste0("(() => { const r = document.querySelector('#integration-cor_scatter img')",
                   ".getBoundingClientRect(); return {x: r.left, y: r.top, w: r.width}; })()"))
  scale <- img$w / cm$dims$width
  px <- img$x + scale * (r$left + (target$x - d$left) / (d$right - d$left) * (r$right - r$left))
  py <- img$y + scale * (r$top + (d$top - target$y) / (d$top - d$bottom) * (r$bottom - r$top))
  b$Input$dispatchMouseEvent(type = "mouseMoved", x = px - 4, y = py - 4)
  b$Input$dispatchMouseEvent(type = "mouseMoved", x = px, y = py)
  app$wait_for_js(
    "document.querySelector('#integration-cor_scatter_tip .plot-hover-tip') !== null",
    timeout = 15000)
  tip <- js("document.querySelector('#integration-cor_scatter_tip .plot-hover-tip').innerText")
  expect_match(tip, target$name, fixed = TRUE)
  expect_match(tip, "Spearman r", fixed = TRUE)
  expect_match(tip, "adjusted p", fixed = TRUE)

  # The figure's PNG, named for the project, the figure, the layers and
  # the day.
  f <- app$get_download("integration-cor_scatter_download-png")
  expect_identical(readBin(f, "raw", 8L),
                   as.raw(c(0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A)))
  expect_gt(file.size(f), 10000)
  expect_identical(basename(f), sprintf("Paired_donors_correlation_per_gene_proteomics_rnaseq_%s.png",
                                        format(Sys.Date(), "%Y%m%d")))
})
