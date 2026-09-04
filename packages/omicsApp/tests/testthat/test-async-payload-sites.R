# What each view actually hands to the background worker.
#
# future serialises the function it is given, environment and all. A
# closure written inline in a module server therefore carried that
# server's whole scope -- previous bundle, project, demo fixtures -- and
# on a real workbook reached 527 MB, which future refused to export,
# from inside an observer, which ended the session. detached_call()
# exists to stop that, and test-async-payload.R proves it works in
# isolation. This file checks the three places it is *used*: if someone
# adds a fourth run_async() call site and forgets it, the payload grows
# silently and nothing else here would notice.

payload_sizes <- function(code) {
  sizes <- list()
  testthat::local_mocked_bindings(
    run_async = function(func, on_success, on_error, message = "Running...",
                         .future = NULL) {
      # Without source references: under load_all a function carries
      # its whole file's parse data (300 KB here), which an installed
      # package does not. The worker pays for the data, not the srcref.
      sizes[[length(sizes) + 1L]] <<- length(serialize(utils::removeSource(func), NULL))
      result <- tryCatch(func(), error = function(e) e)
      if (inherits(result, "error")) on_error(conditionMessage(result))
      else on_success(result)
      invisible(NULL)
    },
    .package = "omicsApp",
    .env = parent.frame()
  )
  force(code)
  unlist(sizes)
}

bytes <- function(x) length(serialize(x, NULL))

tiny_project <- function() {
  xlsx <- tempfile(fileext = ".xlsx")
  write_tiny_omics_xlsx(xlsx, n_features = 30L, n_samples = 8L)
  inp <- omicsCore::read_omics(xlsx, omics_type = "proteomics",
                               assay_type = "normalized_intensity")$input
  unlink(xlsx)
  omicsCore::omics_project(name = "payload",
                           experiments = list(proteomics = inp))
}

# A payload that carries only what it needs is about the size of the
# data plus the arguments. One that carries the module scope is that
# plus every fixture the module can see -- the demo project alone is
# hundreds of kilobytes here and hundreds of megabytes on real data.
# The bound is data-relative so it does not drift with the fixture.
expect_carries_only <- function(sizes, data, slack = 200e3) {
  expect_gt(length(sizes), 0L)
  for (s in sizes) {
    expect_lt(s, bytes(data) + slack)
  }
}

test_that("the differential view sends the layer and its arguments, not its scope", {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("readxl")
  proj <- tiny_project()
  inp <- proj$experiments$proteomics
  sizes <- payload_sizes(
    shiny::testServer(
      diff_view_server,
      args = list(current_project = shiny::reactiveVal(proj)),
      {
        session$setInputs(method = "ttest", group_col = "group",
                          control = "G1", case = "G2", rerun = 1)
        expect_s3_class(diff_bundle(), "analysis_bundle")
      }
    )
  )
  expect_carries_only(sizes, inp)
})

test_that("the enrichment view sends the bundle and its thresholds, not its scope", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("msigdbr")
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("readxl")
  proj <- tiny_project()
  bundle <- omicsCore::run_diff(
    proj$experiments$proteomics, method = "ttest", analysis_type = "group",
    group_col = "group", control_group = "G1", case_group = "G2"
  )
  sizes <- payload_sizes(
    shiny::testServer(
      enrich_view_server,
      args = list(diff_bundle = shiny::reactiveVal(bundle)),
      {
        session$setInputs(type = "ora", database = "hallmark",
                          direction = "both", rerun = 1)
      }
    )
  )
  expect_carries_only(sizes, bundle)
})

test_that("the integration view sends the project, not its scope", {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("readxl")
  xlsx <- tempfile(fileext = ".xlsx")
  on.exit(unlink(xlsx), add = TRUE)
  write_tiny_omics_xlsx(xlsx, n_features = 30L, n_samples = 8L)
  prot <- omicsCore::read_omics(xlsx, omics_type = "proteomics",
                                assay_type = "normalized_intensity")$input
  rna <- omicsCore::read_omics(xlsx, omics_type = "rnaseq",
                               assay_type = "logcpm")$input
  proj <- omicsCore::omics_project(
    name = "two", experiments = list(proteomics = prot, rnaseq = rna)
  )
  primary <- omicsCore::run_diff(
    prot, method = "ttest", analysis_type = "group",
    group_col = "group", control_group = "G1", case_group = "G2"
  )
  sizes <- payload_sizes(
    shiny::testServer(
      integration_view_server,
      args = list(current_project = shiny::reactiveVal(proj),
                  diff_bundle = shiny::reactiveVal(primary)),
      {
        session$setInputs(rerun = 0)
        expect_true(isTRUE(can_run()$ok))
        session$setInputs(rerun = 1)
        expect_s3_class(integration_bundle(), "analysis_bundle")
      }
    )
  )
  # The worker gets the secondary layer and the primary bundle, which
  # together are about the project plus one result.
  expect_carries_only(sizes, list(proj, primary))
})

test_that("every run_async() call site in the package goes through detached_call()", {
  # The three tests above cover the call sites that exist today. This
  # keeps a fourth from arriving without one.
  src_dir <- system.file("R", package = "omicsApp") |> dirname() |>
    file.path("R")
  src <- list.files(src_dir, pattern = "\\.R$", full.names = TRUE)
  skip_if(length(src) == 0L, "package sources not available")
  for (f in src) {
    if (basename(f) == "async_helpers.R") next
    lines <- readLines(f, warn = FALSE)
    lines[grepl("^\\s*#", lines)] <- ""   # comments and roxygen mention it too
    calls <- grep("run_async(", lines, fixed = TRUE)
    for (i in calls) {
      window <- paste(lines[i:min(i + 12L, length(lines))], collapse = "\n")
      expect_true(grepl("detached_call(", window, fixed = TRUE),
                  info = sprintf("%s:%d", basename(f), i))
    }
  }
})
