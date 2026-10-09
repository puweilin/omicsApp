# The Differential view's Heatmap card is drawn once it comes into view,
# not when the view opens (R/mod_diff_detail.R): it sits below the
# volcano and the table, and clustering and drawing it added about 0.4 s
# to opening the view. Its download is offered all the same, and draws
# the figure when asked for.

lazy_project <- function() {
  set.seed(11)
  groups <- rep(c("Ctrl", "Treat"), each = 4)
  samp <- paste0("S", seq_along(groups))
  ids <- paste0("P", 1:40)
  m <- matrix(stats::rnorm(40 * length(groups), 12, 0.3), 40, dimnames = list(ids, samp))
  m[1:10, groups == "Treat"] <- m[1:10, groups == "Treat"] + 2
  inp <- omicsCore::omics_input(
    m, data.frame(group = groups, row.names = samp, stringsAsFactors = FALSE),
    data.frame(feature_id = ids, feature_symbol = paste0("G", 1:40), row.names = ids,
               stringsAsFactors = FALSE),
    omics_type = "proteomics", assay_type = "normalized_intensity")
  omicsCore::omics_project("Lazy", list(proteomics = inp))
}

test_that("the heatmap is not drawn until its card is in view; its download still works", {
  calls <- 0L
  real <- omicsCore::plot_heatmap
  testthat::local_mocked_bindings(
    plot_heatmap = function(...) {
      calls <<- calls + 1L
      real(...)
    },
    .package = "omicsCore")
  shiny::testServer(diff_view_server,
                    args = list(current_project = shiny::reactiveVal(lazy_project())), {
    session$setInputs(group_col = "group", control = "Ctrl", case = "Treat",
                      method = "limma", fdr_cut = 0.05, fc_cut = 0.263, p_kind = "adj",
                      rerun = 1)
    session$elapse(300)
    session$flushReact()
    expect_true(detail$has_heatmap())
    # The card's panel and note are there; the figure is not drawn.
    expect_match(output$heatmap_note$html, "scaled to its own mean")
    expect_error(output$heatmap)
    expect_identical(calls, 0L)
    # The download menu is offered without drawing the figure ...
    menu <- paste(as.character(output[["heatmap_download-menu"]]$html), collapse = "")
    expect_match(menu, "dropdown-item")
    expect_identical(calls, 0L)
    # ... and the file is drawn when asked for.
    f <- output[["heatmap_download-png"]]
    expect_true(file.exists(f))
    expect_identical(readBin(f, "raw", 8L)[2:4], charToRaw("PNG"))
    expect_identical(calls, 1L)

    # Scrolled to: the card draws it.
    session$setInputs(heatmap_visible = TRUE)
    expect_match(output$heatmap$src, "^data:image/png")
    # Once seen it follows the selection like any other figure.
    session$setInputs(hits_rows_selected = 1L)
    expect_match(output$heatmap$src, "^data:image/png")
  })
})

test_that("nothing to draw: no download, and no wait for the card either", {
  shiny::testServer(diff_view_server,
                    args = list(current_project = shiny::reactiveVal(lazy_project())), {
    session$setInputs(group_col = "group", control = "Ctrl", case = "Treat",
                      method = "limma", fdr_cut = 1e-300, fc_cut = 0.263, p_kind = "adj",
                      rerun = 1)
    session$elapse(300)
    session$flushReact()
    expect_false(detail$heatmap_available())
    expect_identical(paste(as.character(output[["heatmap_download-menu"]]$html), collapse = ""), "")
  })
})

test_that("the card tells the server when it comes into view, once", {
  html <- as.character(diff_view_ui("diff"))
  expect_match(html, 'id="diff-heatmap_card"', fixed = TRUE)
  js <- as.character(when_visible_script("diff-heatmap_card", "diff-heatmap_visible"))
  expect_match(js, "IntersectionObserver", fixed = TRUE)
  expect_match(js, "getElementById('diff-heatmap_card')", fixed = TRUE)
  expect_match(js, "Shiny.setInputValue('diff-heatmap_visible', true)", fixed = TRUE)
  # Sent again after a reconnect, whose new session has not heard it.
  expect_match(js, "shiny:connected", fixed = TRUE)
  # The script is part of the view, not escaped into text.
  expect_match(html, "new IntersectionObserver(function(entries)", fixed = TRUE)
})
