# Every figure card's download menu (R/plot_download.R): PNG at 300 dpi,
# PDF and SVG, drawn from the figure's own ggplot at a print size.

# The figure the tests save: a title with the characters the app's own
# figures use ("≥" in a cut, "·" between facts) and a Chinese
# group name, which the plain pdf() device turned into dots.
download_test_plot <- function() {
  df <- data.frame(x = 1:4, y = c(2, 1, 4, 3),
                   g = c("Control", "Control", "\u5904\u7406", "\u5904\u7406"))
  ggplot2::ggplot(df, ggplot2::aes(x, y, colour = g)) +
    ggplot2::geom_point() +
    ggplot2::labs(title = "|log2FC| \u2265 1 \u00b7 adjusted p < 0.05")
}

# Every file name ends in the day it was saved, as the tables' do.
today <- function() format(Sys.Date(), "%Y%m%d")

png_signature <- as.raw(c(0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A))

# The pixels per metre a PNG's pHYs chunk records, or NULL without one.
png_ppm <- function(path) {
  raw <- readBin(path, "raw", 4096L)
  i <- grepRaw("pHYs", raw, fixed = TRUE)
  if (!length(i)) return(NULL)
  b <- as.integer(raw[(i + 4L):(i + 7L)])
  sum(b * 256^(3:0))
}

expect_valid_png <- function(path, width_in, height_in) {
  expect_true(file.exists(path))
  expect_gt(file.size(path), 1000)
  head <- readBin(path, "raw", 24L)
  expect_identical(head[1:8], png_signature)
  # Pixels = inches x 300: the width and height the IHDR chunk records
  # (read from the bytes, so the check needs no image package).
  big_endian <- function(b) sum(as.numeric(b) * 256^(3:0))
  expect_identical(c(big_endian(head[17:20]), big_endian(head[21:24])),
                   round(c(width_in, height_in) * 300))
  ppm <- png_ppm(path)
  # ragg records the resolution; cairo's PNG writes no pHYs at all.
  if (!is.null(ppm)) expect_equal(ppm / 39.3701, 300, tolerance = 0.01)
}

expect_valid_pdf <- function(path) {
  expect_gt(file.size(path), 1000)
  expect_identical(rawToChar(readBin(path, "raw", 5L)), "%PDF-")
  # cairo embeds its fonts, so the file reads the same on a machine
  # without them.
  expect_true(length(grepRaw("FontFile", readBin(path, "raw", file.size(path)))) > 0L)
}

# Parsed with xml2 where it is installed; it is no dependency of the app,
# so it is looked up rather than called with `::`.
expect_valid_svg <- function(path) {
  skip_if_not_installed("xml2")
  xml <- function(f) getExportedValue("xml2", f)
  expect_gt(file.size(path), 1000)
  doc <- xml("read_xml")(path)
  expect_identical(xml("xml_name")(doc), "svg")
  find <- xml("xml_find_all")
  if (has_pkg("svglite")) {
    # svglite keeps the labels as text, the title among them.
    text <- xml("xml_text")(find(doc, "//*[local-name()='text']"))
    expect_gt(length(text), 3L)
    return(invisible())
  }
  # cairo writes text as glyph outlines, each glyph defined once and
  # placed with <use>: text is drawn when glyphs are placed.
  glyphs <- find(doc, "//*[local-name()='g'][starts-with(@id, 'glyph')]")
  uses <- find(doc, "//*[local-name()='use']")
  expect_gt(length(glyphs), 10L)
  expect_gt(length(uses), 20L)
}

test_that("each format is a valid file of the asked size, through the menu's handlers", {
  shiny::testServer(plot_download_server, args = list(
    plot_reactive = shiny::reactive(download_test_plot()),
    filename_stem = function() c("Cheek \u00b7 \u4e2d\u6587", "pca", "proteomics"),
    width_in = 6, height_in = function() 4
  ), {
    session$flushReact()
    expect_valid_png(output$png, 6, 4)
    # The file is named for its parts, made safe, Chinese kept, and dated.
    expect_identical(basename(output$png),
                     sprintf("Cheek_\u00b7_\u4e2d\u6587_pca_proteomics_%s.png", today()))
    expect_identical(basename(output$svg),
                     sprintf("Cheek_\u00b7_\u4e2d\u6587_pca_proteomics_%s.svg", today()))
    expect_valid_pdf(output$pdf)
    expect_valid_svg(output$svg)
  })
})

test_that("the cairo PNG is used where ragg is not installed, at the same size", {
  local_mocked_bindings(has_pkg = function(pkg) !identical(pkg, "ragg"),
                        .package = "omicsApp")
  f <- withr::local_tempfile(fileext = ".png")
  save_plot_file(download_test_plot(), f, "png", width_in = 3, height_in = 2)
  expect_valid_png(f, 3, 2)
})

test_that("the SVG goes through svglite, its text kept editable, where svglite is installed", {
  # svglite is not installed on every machine the tests run on, so a
  # stand-in with its arguments records the call and draws the file with
  # cairo, as svglite would draw it.
  seen <- NULL
  fake_svglite <- function(filename = "Rplot%03d.svg", width = 10, height = 8,
                           bg = "white", pointsize = 12, standalone = TRUE,
                           fix_text_size = TRUE, file) {
    seen <<- as.list(match.call())[-1]
    grDevices::svg(if (missing(file)) filename else file, width = width,
                   height = height, bg = bg)
  }
  local_mocked_bindings(has_pkg = function(pkg) TRUE,
                        svglite_device = function() fake_svglite,
                        .package = "omicsApp")
  f <- withr::local_tempfile(fileext = ".svg")
  save_plot_file(download_test_plot(), f, "svg", width_in = 5, height_in = 3)
  expect_identical(seen$file, f)
  expect_equal(c(seen$width, seen$height), c(5, 3))
  expect_identical(seen$bg, "white")
  # Labels left at their own width, so a retyped one is not stretched.
  expect_false(seen$fix_text_size)
  expect_gt(file.size(f), 1000)

  # An svglite before 2.0, without fix_text_size, is not given it.
  old_svglite <- function(file = "Rplot%03d.svg", width = 10, height = 8,
                          bg = "white", pointsize = 12, standalone = TRUE) {
    seen <<- as.list(match.call())[-1]
    grDevices::svg(file, width = width, height = height, bg = bg)
  }
  local_mocked_bindings(svglite_device = function() old_svglite, .package = "omicsApp")
  save_plot_file(download_test_plot(), f, "svg", width_in = 5, height_in = 3)
  expect_identical(seen$file, f)
  expect_null(seen$fix_text_size)
})

test_that("the cairo SVG is used where svglite is not installed", {
  local_mocked_bindings(has_pkg = function(pkg) !identical(pkg, "svglite"),
                        svglite_device = function() stop("svglite is not installed"),
                        .package = "omicsApp")
  expect_identical(plot_svg_device(), grDevices::svg)
  f <- withr::local_tempfile(fileext = ".svg")
  save_plot_file(download_test_plot(), f, "svg", width_in = 5, height_in = 3)
  expect_valid_svg(f)
})

test_that("a patchwork figure (the QC quality panels) saves in every format", {
  p <- omicsCore::plot_qc(example_qc_bundle(), view = "missing")
  expect_s3_class(p, "patchwork")
  shiny::testServer(plot_download_server, args = list(
    plot_reactive = shiny::reactive(p), filename_stem = function() "qc",
    width_in = 8, height_in = 4.5
  ), {
    session$flushReact()
    expect_valid_png(output$png, 8, 4.5)
    expect_valid_pdf(output$pdf)
    expect_valid_svg(output$svg)
  })
})

test_that("the menu offers the three formats, labelled, while there is a figure", {
  shiny::testServer(plot_download_server, args = list(
    plot_reactive = shiny::reactive(download_test_plot()), filename_stem = function() "x"
  ), {
    session$flushReact()
    html <- paste(as.character(output$menu$html), collapse = "")
    for (fmt in c("png", "pdf", "svg")) {
      expect_match(html, sprintf('id="%s"', session$ns(fmt)), fixed = TRUE)
    }
    expect_match(html, "PNG (300 dpi)", fixed = TRUE)
    expect_match(html, 'aria-label="Download figure"', fixed = TRUE)
    expect_match(html, 'data-bs-toggle="dropdown"', fixed = TRUE)
    expect_match(html, "shiny-download-link[^\"]*dropdown-item")
  })
})

test_that("no menu before there is a figure, for an 'empty' figure, or when the caller says so", {
  waiting <- shiny::reactiveVal(NULL)
  shiny::testServer(plot_download_server, args = list(
    plot_reactive = shiny::reactive({ shiny::req(waiting()); waiting() }),
    filename_stem = function() "x"
  ), {
    session$flushReact()
    expect_identical(paste(as.character(output$menu$html), collapse = ""), "")
    # omicsCore's "No pathways to plot." figure is a sentence, not a result.
    waiting(omicsCore:::empty_plot("No pathways to plot."))
    session$flushReact()
    expect_identical(paste(as.character(output$menu$html), collapse = ""), "")
    waiting(download_test_plot())
    session$flushReact()
    expect_match(paste(as.character(output$menu$html), collapse = ""), "dropdown-item")
  })
  shiny::testServer(plot_download_server, args = list(
    plot_reactive = shiny::reactive(download_test_plot()), filename_stem = function() "x",
    available = shiny::reactive(FALSE)
  ), {
    session$flushReact()
    expect_identical(paste(as.character(output$menu$html), collapse = ""), "")
  })
})

test_that("an empty figure is told from a real one", {
  expect_true(is_empty_plot(omicsCore:::empty_plot("No integration rows to plot.")))
  expect_false(is_empty_plot(download_test_plot()))
  expect_false(is_empty_plot(omicsCore::plot_qc(example_qc_bundle(), view = "missing")))
  expect_false(plot_ready(function() stop("not yet")))
  expect_true(plot_ready(function() download_test_plot()))
})

test_that("file names are the project, the figure and the comparison, made safe, and dated", {
  day <- as.Date("2026-10-09")
  expect_identical(plot_download_name(c("My project", "volcano", "TreatA_vs_Control"), day),
                   "My_project_volcano_TreatA_vs_Control_20261009")
  # Chinese kept whole; slashes and the characters Windows refuses taken out.
  expect_identical(plot_download_name(c("\u4e2d\u6587 \u9879\u76ee", "heatmap", "A/B: C?"), day),
                   "\u4e2d\u6587_\u9879\u76ee_heatmap_A_B_C_20261009")
  # Empty and missing parts are dropped; nothing at all is "figure".
  expect_identical(plot_download_name(list(NULL, "pca", NA, ""), day), "pca_20261009")
  expect_identical(plot_download_name(NULL, day), "figure_20261009")
  # Today's date by default, in the tables' format; none when asked.
  expect_identical(plot_download_name("pca"), paste0("pca_", today()))
  expect_identical(plot_download_name("pca", date = NULL), "pca")
  # A long name is shortened before the date goes on: the date survives.
  long <- plot_download_name(c(strrep("a", 90), strrep("b", 90)), day)
  expect_lte(nchar(long, type = "bytes"), 150L)
  expect_match(long, "^a+_b+_20261009$")
  long_cn <- plot_download_name(strrep("\u4e2d", 80), day)
  expect_lte(nchar(long_cn, type = "bytes"), 150L)
  expect_true(validUTF8(long_cn))
  expect_match(long_cn, "_20261009$")
  expect_identical(plot_download_size(function() 6.5, 7), 6.5)
  expect_identical(plot_download_size(function() stop("no"), 7), 7)
  expect_identical(plot_download_size(NA, 4.5), 4.5)
})

test_that("the QC view offers its figures for a project and not for the demo", {
  demo <- shiny::reactiveVal(NULL)
  shiny::testServer(qc_view_server, args = list(current_project = demo), {
    session$setInputs(missing_threshold = 0.5, outlier_method = "iqr")
    expect_true(inherits(pca_plot(), "ggplot"))
    expect_identical(paste(as.character(output[["pca_download-menu"]]$html), collapse = ""), "")
  })
  proj <- shiny::reactiveVal(example_project())
  shiny::testServer(qc_view_server, args = list(current_project = proj), {
    session$setInputs(missing_threshold = 0.5, outlier_method = "iqr")
    session$flushReact()
    expect_match(paste(as.character(output[["pca_download-menu"]]$html), collapse = ""),
                 "dropdown-item")
    expect_match(paste(as.character(output[["missing_download-menu"]]$html), collapse = ""),
                 "dropdown-item")
    f <- output[["pca_download-png"]]
    expect_valid_png(f, 7, 5)
    expect_identical(basename(f), sprintf("%s_pca_%s_%s.png", project_slug(example_project()$name),
                                          active()$tag, today()))
    # The quality panel is named for the panel shown.
    session$setInputs(quality_view = "depth")
    expect_match(basename(output[["missing_download-pdf"]]),
                 sprintf("_intensity_proteomics_%s\\.pdf$", today()))
  })
})

# A small three-group project with clear hits, for the Differential view.
dl_diff_project <- function() {
  set.seed(7)
  groups <- rep(c("DrugA", "DMSO", "DrugB"), each = 4)
  samp <- paste0("S", seq_along(groups))
  ids <- paste0("P", 1:60)
  m <- matrix(stats::rnorm(60 * length(groups), 12, 0.3), 60, dimnames = list(ids, samp))
  m[1:8, groups == "DrugA"] <- m[1:8, groups == "DrugA"] + 2.5
  m[9:14, groups == "DrugA"] <- m[9:14, groups == "DrugA"] - 2.5
  inp <- omicsCore::omics_input(
    m, data.frame(treatment = groups, row.names = samp, stringsAsFactors = FALSE),
    data.frame(feature_id = ids, feature_symbol = paste0("G", 1:60), row.names = ids,
               stringsAsFactors = FALSE),
    omics_type = "proteomics", assay_type = "normalized_intensity")
  omicsCore::omics_project("Drug screen", list(proteomics = inp))
}

test_that("the Differential view's figures are saved at their sizes, with the comparison's name", {
  shiny::testServer(diff_view_server,
                    args = list(current_project = shiny::reactiveVal(dl_diff_project())), {
    menu <- function(id) paste(as.character(output[[paste0(id, "-menu")]]$html), collapse = "")
    session$flushReact()
    # Before a run there is nothing to save.
    expect_identical(menu("volcano_download"), "")
    expect_identical(menu("heatmap_download"), "")

    session$setInputs(group_col = "treatment", control = "DMSO", case = c("DrugA", "DrugB"),
                      method = "limma", fdr_cut = 0.05, fc_cut = 0.263, p_kind = "adj",
                      label_top = TRUE, rerun = 1)
    session$elapse(300)
    session$flushReact()
    expect_match(menu("volcano_download"), "dropdown-item")
    expect_match(menu("heatmap_download"), "dropdown-item")
    expect_match(menu("contrast_plot_download"), "Download the hits per comparison")
    expect_match(menu("overlap_plot_download"), "dropdown-item")
    # Nothing selected: no feature figure to save.
    expect_identical(menu("feature_plot_download"), "")

    # The volcano as plot_volcano() draws it, the top hits named.
    v <- results$volcano_plot()
    expect_true(any(vapply(v$layers, function(l) grepl("Repel", class(l$geom)[1]), logical(1))))
    f <- output[["volcano_download-pdf"]]
    expect_valid_pdf(f)
    cmp <- shown_bundle()$params$comparison
    expect_identical(basename(f), sprintf("Drug_screen_volcano_%s_%s.pdf", cmp, today()))

    # The heatmap as tall as the card draws it for its rows.
    n <- length(detail$heatmap_hits()$ids)
    f <- output[["heatmap_download-png"]]
    expect_valid_png(f, 8, heatmap_height(n) / PLOT_RES)

    session$setInputs(overlap_dir = "up")
    expect_match(basename(output[["overlap_plot_download-svg"]]),
                 sprintf("^Drug_screen_overlap_up_proteomics_%s\\.svg$", today()))

    # A selected feature's figure is named for the feature.
    session$setInputs(hits_rows_selected = 1L)
    session$flushReact()
    expect_match(menu("feature_plot_download"), "dropdown-item")
    sym <- results$hits_df()$feature_symbol[[1]]
    expect_identical(basename(output[["feature_plot_download-png"]]),
                     sprintf("Drug_screen_%s_%s_%s.png", sym, cmp, today()))
  })
})

test_that("the demo project's figures have no download", {
  shiny::testServer(diff_view_server, args = list(current_project = shiny::reactiveVal(NULL)), {
    session$flushReact()
    expect_true(isTRUE(active()$is_demo))
    expect_identical(paste(as.character(output[["volcano_download-menu"]]$html), collapse = ""), "")
  })
  shiny::testServer(integration_view_server,
                    args = list(current_project = shiny::reactiveVal(NULL)), {
    session$flushReact()
    expect_true(inherits(figures$scatter(), "ggplot"))
    expect_identical(paste(as.character(output[["scatter_download-menu"]]$html), collapse = ""), "")
  })
})
