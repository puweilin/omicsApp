# The Project view's Experiments table used to render a bare <a> with
# no href and no handler: it looked like a link, was styled like one,
# and did nothing. These cover the wiring that replaced it.

demo_two_layer <- function() {
  mat <- function(n) {
    m <- matrix(as.numeric(seq_len(n * 4)), nrow = n,
                dimnames = list(paste0("f", seq_len(n)), paste0("s", 1:4)))
    m
  }
  meta <- data.frame(group = c("A", "A", "B", "B"),
                     row.names = paste0("s", 1:4))
  prot <- omicsCore::omics_input(mat(6), meta,
            data.frame(feature_id = paste0("f", 1:6)),
            omics_type = "proteomics", assay_type = "normalized_intensity")
  rna <- omicsCore::omics_input(round(mat(5)), meta,
            data.frame(feature_id = paste0("f", 1:5)),
            omics_type = "rnaseq", assay_type = "raw_count")
  omicsCore::omics_project(name = "two",
                           experiments = list(proteomics = prot, rnaseq = rna))
}

# ---- the card ---------------------------------------------------------

test_that("the View cell is an actual control when a namespace is given", {
  exps <- demo_two_layer()$experiments
  html <- as.character(project_experiments_card(exps, ns = function(x) paste0("p-", x)))
  # An id per row, keyed by position rather than by tag: a tag is
  # user-supplied text and would make an unsafe input id.
  expect_true(grepl('id="p-view_layer_1"', html, fixed = TRUE))
  expect_true(grepl('id="p-view_layer_2"', html, fixed = TRUE))
  expect_true(grepl("action-button", html, fixed = TRUE))
})

test_that("without a namespace the cell is inert rather than a fake link", {
  # Rendered where there is no session to click in, the cell must not
  # advertise an action it cannot perform.
  html <- as.character(project_experiments_card(demo_two_layer()$experiments))
  expect_false(grepl("action-button", html, fixed = TRUE))
  expect_false(grepl("view_layer_", html, fixed = TRUE))
})

test_that("the card carries no literal non-ASCII into the page", {
  # R CMD check fails on non-ASCII outside comments, and the arrow is
  # in a string.
  src <- readLines(system.file("R", package = "omicsApp") |> dirname() |>
                     file.path("R", "mod_project_view.R"), warn = FALSE)
  skip_if(length(src) == 0, "installed package has no sources")
  expect_false(any(grepl("[^\x01-\x7f]", src)))
})

# ---- the click --------------------------------------------------------

test_that("clicking a row reports that row's tag", {
  clicked <- NULL
  shiny::testServer(
    project_view_server,
    args = list(current_project = shiny::reactiveVal(demo_two_layer()),
                on_view_layer = function(tag) clicked <<- tag),
    {
      session$flushReact()
      session$setInputs(view_layer_2 = 1)
      expect_identical(clicked, "rnaseq")
      session$setInputs(view_layer_1 = 1)
      expect_identical(clicked, "proteomics")
    }
  )
})

test_that("a click on a row that no longer exists is ignored", {
  # The project can be replaced between the link being drawn and
  # clicked. Reporting a tag from the old project would send QC to a
  # layer that is not there.
  clicked <- "untouched"
  proj <- shiny::reactiveVal(demo_two_layer())
  shiny::testServer(
    project_view_server,
    args = list(current_project = proj,
                on_view_layer = function(tag) clicked <<- tag),
    {
      session$flushReact()
      one <- demo_two_layer()
      one$experiments <- one$experiments["proteomics"]
      proj(one)
      session$flushReact()
      session$setInputs(view_layer_2 = 1)
      expect_identical(clicked, "untouched")
    }
  )
})

# ---- what QC does with it ---------------------------------------------

test_that("QC shows the requested layer instead of its default pick", {
  # Left alone QC picks the first proteomics layer, so asking for
  # rnaseq is a request it would not have satisfied by accident.
  shiny::testServer(
    qc_view_server,
    args = list(current_project = shiny::reactiveVal(demo_two_layer()),
                requested_layer = shiny::reactiveVal("rnaseq")),
    {
      expect_identical(active()$tag, "rnaseq")
    }
  )
})

test_that("no request means the previous behaviour, unchanged", {
  shiny::testServer(
    qc_view_server,
    args = list(current_project = shiny::reactiveVal(demo_two_layer())),
    {
      expect_identical(active()$tag, "proteomics")
    }
  )
})

test_that("a request naming a layer that is gone falls back rather than emptying", {
  # Stale requests are the normal case after a re-import. Honouring one
  # would leave the view blank with nothing to say why.
  shiny::testServer(
    qc_view_server,
    args = list(current_project = shiny::reactiveVal(demo_two_layer()),
                requested_layer = shiny::reactiveVal("no_such_layer")),
    {
      expect_identical(active()$tag, "proteomics")
      expect_false(active()$is_demo)
    }
  )
})
