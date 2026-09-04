# The first design goal in the README is that omicsCore installs light:
# the Bioconductor stack is Suggests, and a feature whose package is
# missing tells the user what to install instead of failing somewhere
# inside a stack trace.
#
# Nothing had ever run that path. Every developer machine and every CI
# runner has the full stack, and the tests that were meant to cover the
# gates all said "X is installed; cannot test error path" -- seven
# skips, on every run, everywhere. The gates now ask `is_installed()`,
# which a test can answer, so both branches run on the same machine.

absent <- function(...) {
  gone <- c(...)
  testthat::local_mocked_bindings(
    is_installed = function(pkg) {
      !(pkg %in% gone) && requireNamespace(pkg, quietly = TRUE)
    },
    .package = "omicsCore",
    .env = parent.frame()
  )
}

tiny_input <- function(omics_type = "proteomics") {
  set.seed(11)
  n <- 8L
  mat <- matrix(stats::rnorm(20 * n, mean = 18, sd = 1), 20, n,
                dimnames = list(sprintf("F%02d", 1:20), sprintf("S%d", 1:n)))
  assay <- "normalized_intensity"
  if (omics_type == "rnaseq") {
    mat <- matrix(stats::rpois(20 * n, 200), 20, n, dimnames = dimnames(mat))
    storage.mode(mat) <- "integer"
    assay <- "raw_count"
  }
  meta <- data.frame(group = rep(c("A", "B"), each = n / 2),
                     age = seq_len(n) + 30,
                     row.names = colnames(mat))
  feat <- data.frame(feature_id = rownames(mat),
                     feature_symbol = rownames(mat))
  omics_input(mat, meta, feat, omics_type = omics_type, assay_type = assay)
}

# ---- the differential engines -----------------------------------------

test_that("each differential engine names itself and its install group when absent", {
  cases <- list(
    list(pkg = "limma",  method = "limma",  omics = "proteomics", group = "proteomics"),
    list(pkg = "DESeq2", method = "deseq2", omics = "rnaseq",     group = "rnaseq"),
    list(pkg = "edgeR",  method = "edger",  omics = "rnaseq",     group = "rnaseq")
  )
  for (cs in cases) {
    absent(cs$pkg)
    err <- tryCatch(
      run_diff(tiny_input(cs$omics), method = cs$method,
               analysis_type = "group", group_col = "group",
               control_group = "A", case_group = "B"),
      error = function(e) conditionMessage(e)
    )
    expect_type(err, "character")
    expect_match(err, cs$pkg, fixed = TRUE, info = cs$pkg)
    expect_match(err, paste0("install_optional('", cs$group, "')"),
                 fixed = TRUE, info = cs$pkg)
  }
})

test_that("method = 'auto' steps down to an engine that is present", {
  # This is the promise the lightweight install makes: the analysis
  # still runs, on a simpler engine, rather than refusing.
  absent("limma")
  b <- run_diff(tiny_input(), method = "auto", analysis_type = "group",
                group_col = "group", control_group = "A", case_group = "B")
  expect_identical(b$params$method, "ttest")

  absent("DESeq2", "edgeR", "limma")
  b <- run_diff(tiny_input("rnaseq"), method = "auto",
                analysis_type = "group", group_col = "group",
                control_group = "A", case_group = "B")
  expect_false(b$params$method %in% c("deseq2", "edger", "limma"))
  expect_gt(nrow(b$results$diff_result_df), 0L)
})

# ---- enrichment, scoring, normalization, persistence ------------------

test_that("the enrichment layer lists every missing package at once", {
  absent("clusterProfiler", "msigdbr")
  b <- run_diff(tiny_input(), method = "ttest", analysis_type = "group",
                group_col = "group", control_group = "A", case_group = "B")
  err <- tryCatch(run_enrichment(b, type = "ora", database = "hallmark"),
                  error = function(e) conditionMessage(e))
  expect_match(err, "clusterProfiler", fixed = TRUE)
  expect_match(err, "msigdbr", fixed = TRUE)
  expect_match(err, "install_optional('enrichment')", fixed = TRUE)
})

test_that("run_gsva() refuses without GSVA and says which group installs it", {
  absent("GSVA")
  err <- tryCatch(run_gsva(tiny_input(), database = "hallmark"),
                  error = function(e) conditionMessage(e))
  expect_match(err, "GSVA", fixed = TRUE)
  expect_match(err, "install_optional('enrichment')", fixed = TRUE)
})

test_that("normalize_omics() refuses without vsn", {
  absent("vsn")
  inp <- tiny_input()
  inp$expr_mat <- 2^inp$expr_mat
  inp$assay_type <- "raw_intensity"
  err <- tryCatch(normalize_omics(inp, method = "vsn"),
                  error = function(e) conditionMessage(e))
  expect_match(err, "vsn", fixed = TRUE)
})

test_that("save_project() and load_project() refuse without qs2", {
  absent("qs2")
  p <- omics_project("p", experiments = list(proteomics = tiny_input()))
  path <- tempfile(fileext = ".omp")
  err <- tryCatch(save_project(p, path), error = function(e) conditionMessage(e))
  expect_match(err, "qs2", fixed = TRUE)
  expect_match(err, "persistence", fixed = TRUE)
  expect_false(file.exists(path))
})

test_that("export_report() refuses without rmarkdown", {
  absent("rmarkdown")
  p <- omics_project("p", experiments = list(proteomics = tiny_input()))
  err <- tryCatch(export_report(p, tempfile(fileext = ".html")),
                  error = function(e) conditionMessage(e))
  expect_match(err, "rmarkdown", fixed = TRUE)
})

test_that("each imputation method names its package when absent", {
  inp <- tiny_input()
  inp$expr_mat[1:3, 1:2] <- NA
  cases <- c(knn = "impute", missforest = "missForest", bpca = "pcaMethods")
  for (method in names(cases)) {
    absent(cases[[method]])
    err <- tryCatch(impute_matrix(inp$expr_mat, method = method),
                    error = function(e) conditionMessage(e))
    expect_type(err, "character")
    expect_match(err, cases[[method]], fixed = TRUE, info = method)
  }
})

test_that("ActivePathways integration refuses without the package", {
  absent("ActivePathways")
  err <- tryCatch(ensure_active_pathways(),
                  error = function(e) conditionMessage(e))
  expect_match(err, "ActivePathways", fixed = TRUE)
})

# ---- plots degrade rather than fail -----------------------------------

test_that("the volcano still draws without ggrepel", {
  absent("ggrepel")
  b <- run_diff(tiny_input(), method = "ttest", analysis_type = "group",
                group_col = "group", control_group = "A", case_group = "B")
  p <- plot_volcano(b, top_n = 5)
  expect_s3_class(p, "ggplot")
  classes <- unlist(lapply(p$layers, function(l) class(l$geom)))
  expect_false(any(grepl("Repel", classes)))
})

test_that("the heatmap falls back to ggplot without ComplexHeatmap", {
  absent("ComplexHeatmap", "circlize")
  p <- plot_heatmap(tiny_input(), n_top = 10)
  expect_s3_class(p, "ggplot")
  b <- run_diff(tiny_input(), method = "ttest", analysis_type = "group",
                group_col = "group", control_group = "A", case_group = "B")
  p2 <- plot_heatmap(b, input = tiny_input(), n_top = 10)
  expect_s3_class(p2, "ggplot")
})

test_that("the gates read is_installed(), so this file tests what it claims", {
  # If a gate goes back to calling requireNamespace() directly, the mock
  # above stops reaching it and every test here passes vacuously. The
  # only requireNamespace() in the package is the one inside
  # is_installed() itself.
  src <- list.files(system.file("R", package = "omicsCore") |>
                      dirname() |> file.path("R"),
                    pattern = "\\.R$", full.names = TRUE)
  skip_if(length(src) == 0L, "package sources not available")
  hits <- unlist(lapply(src, function(f) {
    lines <- readLines(f, warn = FALSE)
    code <- lines[!grepl("^\\s*#", lines)]
    if (any(grepl("requireNamespace(", code, fixed = TRUE))) basename(f)
  }))
  expect_identical(hits, "install-optional.R")
})
