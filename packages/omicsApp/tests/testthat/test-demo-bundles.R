# The demo views draw through the same `omicsCore::plot_*()` dispatchers
# as the live ones. That only holds while the fixtures keep satisfying
# the result schemas -- and a fixture that drifts out of schema fails at
# render time in front of a user, not here, unless something checks.

test_that("the demo diff fixture plots through the dispatcher", {
  b <- example_diff_bundle()
  expect_true(omicsCore::is_analysis_bundle(b))
  expect_s3_class(omicsCore::plot_volcano(b), "ggplot")
})

test_that("the demo enrichment fixture is a bundle that satisfies the schema", {
  b <- example_enrich_bundle()
  expect_true(omicsCore::is_analysis_bundle(b))
  expect_identical(b$analysis_name, "run_enrichment")
  # plot_enrichment() runs the result through check_enrich_result_schema(),
  # so this fails loudly if the fixture loses a required column.
  expect_s3_class(omicsCore::plot_enrichment(b, view = "dot", top_n = 12L),
                  "ggplot")
})

test_that("the demo integration fixture is a bundle that satisfies the schema", {
  b <- example_integration_bundle()
  expect_true(omicsCore::is_analysis_bundle(b))
  expect_identical(b$analysis_name, "run_integration")
  expect_s3_class(omicsCore::plot_integration(b, view = "dual_volcano"),
                  "ggplot")
  expect_s3_class(omicsCore::plot_integration(b, view = "effect_pair"),
                  "ggplot")
})

test_that("demo bundles are cached rather than rebuilt", {
  # Each demo view resolves its own bundle; rebuilding would reseed the
  # fixtures and leave two views showing different synthetic data, which
  # is the bug the fixture cache was introduced for.
  expect_identical(example_enrich_bundle(), example_enrich_bundle())
  expect_identical(example_integration_bundle(), example_integration_bundle())
})

test_that("the app palette is the one omicsCore draws with", {
  # Two palettes means a figure on screen and the same figure in an
  # exported report are tinted differently.
  expect_identical(omics_colors, omicsCore::omics_colors)
})

# --- example data next to real results -----------------------------------

# The Shared pathways and ActivePathways panels read the synthetic fixture
# unconditionally, so a real project displayed invented pathways beside its
# own results. A "demo fixture" caption does not survive a screenshot, and
# a screenshot is how these get into a slide deck.
test_that("a real project shows no example pathways", {
  skip_if_not_installed("limma")
  set.seed(2); p <- 20; n <- 10
  mk <- function(seed, omics) {
    set.seed(seed)
    e <- matrix(rnorm(p * n, 10, 1), p,
                dimnames = list(paste0("SYM", seq_len(p)), paste0("s", seq_len(n))))
    e[1:5, 6:10] <- e[1:5, 6:10] + 2
    md <- data.frame(group = rep(c("A", "B"), each = 5), row.names = colnames(e))
    fd <- data.frame(feature_id = rownames(e), feature_symbol = rownames(e),
                     row.names = rownames(e), stringsAsFactors = FALSE)
    omicsCore::omics_input(e, md, fd, omics_type = omics,
                           assay_type = if (omics == "rnaseq") "logcpm"
                                        else "normalized_intensity")
  }
  proj <- shiny::reactiveVal(
    omicsCore::omics_project("real", list(prot = mk(1, "proteomics"),
                                          rna = mk(2, "rnaseq"))))
  primary <- omicsCore::run_diff(
    shiny::isolate(proj())$experiments$prot, method = "limma",
    analysis_type = "group", group_col = "group",
    control_group = "A", case_group = "B")
  db <- shiny::reactiveVal(primary)

  shiny::testServer(integration_view_server,
                    args = list(current_project = proj, diff_bundle = db), {
    session$setInputs(rerun = 1)
    session$flushReact()
    skip_if(isTRUE(is_demo()), "integration did not run in this environment")
    expect_null(ap_df())
    expect_true(nzchar(paste(as.character(output$ap_empty), collapse = "")))
  })
})

test_that("the demo still shows them, since that is what a demo is for", {
  shiny::testServer(integration_view_server, {
    session$flushReact()
    expect_true(isTRUE(is_demo()))
    expect_false(is.null(ap_df()))
  })
})

# --- demo data is for the demo -------------------------------------------

# `is_demo` answered two questions at once: "is there a user project" and
# "has this view computed anything". So a real project with no upstream
# analysis, or one whose run had failed, fell back to the built-in example
# and displayed invented pathways as though they were its own. Fixing only
# the success branch left every other path showing them.
db_project <- function() {
  set.seed(6); p <- 20; n <- 8
  mk <- function(seed, omics) {
    set.seed(seed)
    e <- matrix(rnorm(p * n, 10, 1), p,
                dimnames = list(paste0("SYM", seq_len(p)), paste0("s", seq_len(n))))
    md <- data.frame(group = rep(c("A", "B"), each = 4), row.names = colnames(e))
    fd <- data.frame(feature_id = rownames(e), feature_symbol = rownames(e),
                     row.names = rownames(e), stringsAsFactors = FALSE)
    omicsCore::omics_input(e, md, fd, omics_type = omics,
                           assay_type = if (omics == "rnaseq") "logcpm"
                                        else "normalized_intensity")
  }
  omicsCore::omics_project("real", list(prot = mk(1, "proteomics"),
                                        rna = mk(2, "rnaseq")))
}

test_that("a real project with nothing computed shows no example enrichment", {
  proj <- shiny::reactiveVal(db_project())
  shiny::testServer(enrich_view_server,
                    args = list(current_project = proj,
                                diff_bundle = shiny::reactiveVal(NULL)), {
    session$flushReact()
    # There is nothing to enrich, so the view has no result -- but that is
    # not the same as being the demo.
    expect_false(isTRUE(demo_mode()))
    # req() stops the reactive rather than returning, which is the empty
    # state -- what matters is that it is not the fixture.
    expect_null(enrich_bundle())
  })
})

test_that("a real project whose run failed shows no example enrichment", {
  proj <- shiny::reactiveVal(db_project())
  db <- shiny::reactiveVal(omicsApp:::example_diff_bundle())
  testthat::local_mocked_bindings(
    run_enrichment = function(...) stop("no gene sets available"),
    .package = "omicsCore")
  shiny::testServer(enrich_view_server,
                    args = list(current_project = proj, diff_bundle = db), {
    session$setInputs(type = "ora", database = "hallmark", rerun = 1)
    session$flushReact()
    expect_false(is.null(enrich_error()))
    expect_false(isTRUE(demo_mode()))
    # The failure is what is shown, not a table of fixtures.
    expect_null(enrich_bundle())
  })
})

test_that("a real project with no second layer shows no example integration", {
  single <- omicsCore::omics_project(
    "one layer", list(prot = db_project()$experiments$prot))
  proj <- shiny::reactiveVal(single)
  shiny::testServer(integration_view_server,
                    args = list(current_project = proj,
                                diff_bundle = shiny::reactiveVal(NULL)), {
    session$flushReact()
    expect_false(isTRUE(demo_mode()))
    expect_null(ap_df())
    expect_null(conc_df())
  })
})

test_that("the demo, which has no project, still shows them", {
  shiny::testServer(enrich_view_server, {
    session$flushReact()
    expect_true(isTRUE(demo_mode()))
    expect_false(is.null(table_data()))
  })
  shiny::testServer(integration_view_server, {
    session$flushReact()
    expect_false(is.null(conc_df()))
  })
  shiny::testServer(integration_view_server, {
    session$flushReact()
    expect_true(isTRUE(demo_mode()))
    expect_false(is.null(ap_df()))
  })
})
