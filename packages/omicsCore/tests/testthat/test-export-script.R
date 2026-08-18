# `export_script()` promises that the script it emits performs the same
# computation the report describes. A script that reads plausibly but
# computes something else is worse than no script, because a reader
# trusts it -- so the central test here does not inspect the text, it
# runs it in a fresh R process and compares the numbers that come back.

write_tiny_workbook <- function(path, n_features = 8L, n_samples = 6L,
                                seed = 2031L) {
  set.seed(seed)
  feat_ids <- sprintf("P%04d", seq_len(n_features))
  samp_ids <- sprintf("S%02d", seq_len(n_samples))
  expr <- matrix(rnorm(n_features * n_samples, mean = 18, sd = 1.2),
                 nrow = n_features,
                 dimnames = list(feat_ids, samp_ids))
  openxlsx::write.xlsx(
    list(
      expression = data.frame(feature_id = feat_ids, expr,
                              check.names = FALSE,
                              stringsAsFactors = FALSE),
      metadata = data.frame(sample_id = samp_ids,
                            group = rep(c("G1", "G2"), length.out = n_samples),
                            age = seq.int(30L, 30L + n_samples - 1L),
                            stringsAsFactors = FALSE),
      feature_annot = data.frame(feature_id = feat_ids,
                                 feature_symbol = paste0("SYM", seq_len(n_features)),
                                 stringsAsFactors = FALSE)
    ),
    file = path, overwrite = TRUE
  )
  invisible(path)
}

fixture_project <- function(source_path = NULL) {
  mat <- matrix(as.numeric(1:24), nrow = 6,
                dimnames = list(paste0("g", 1:6), paste0("s", 1:4)))
  meta <- data.frame(group = c("G1", "G1", "G2", "G2"),
                     row.names = paste0("s", 1:4))
  feat <- data.frame(feature_id = paste0("g", 1:6))
  inp <- omics_input(mat, meta, feat, omics_type = "proteomics",
                     assay_type = "intensity", source_path = source_path)
  proj <- omics_project("fixture", experiments = list(proteomics = inp))
  proj$bundles <- list()
  proj
}

# ---- rendering --------------------------------------------------------

test_that("export_script rejects anything that is not a project", {
  expect_error(export_script(list(a = 1)), "omics_project")
})

test_that("the header records what produced the script", {
  lines <- export_script(fixture_project("raw/x.xlsx"))
  expect_true(any(grepl("export_script", lines, fixed = TRUE)))
  expect_true(any(grepl("Project : fixture", lines, fixed = TRUE)))
  expect_true(any(grepl("omicsCore", lines, fixed = TRUE)))
  expect_true(any(grepl("^library\\(omicsCore\\)$", lines)))
})

test_that("the input line points at the archived file", {
  lines <- export_script(fixture_project("/srv/omicsapp/users/x/raw/cheek__ab12.xlsx"))
  txt <- paste(lines, collapse = "\n")
  expect_match(txt, 'read_omics(', fixed = TRUE)
  # Rewritten relative to the script, so it resolves wherever the user
  # unpacks the download.
  expect_match(txt, '"raw/cheek__ab12.xlsx"', fixed = TRUE)
  expect_match(txt, 'omics_type = "proteomics"', fixed = TRUE)
})

test_that("a project with no archived file says so instead of guessing", {
  lines <- export_script(fixture_project(NULL))
  txt <- paste(lines, collapse = "\n")
  expect_match(txt, "# NOTE:", fixed = TRUE)
  expect_match(txt, "was not archived", fixed = TRUE)
  # The note must precede the code a reader would otherwise run first.
  expect_lt(grep("# NOTE:", lines)[[1L]],
            grep("^library\\(omicsCore\\)$", lines)[[1L]])
})

test_that("derived values are not emitted as arguments", {
  proj <- fixture_project("raw/x.xlsx")
  proj$bundles <- list(diff = new_analysis_bundle(
    analysis_name = "run_diff",
    input_info = list(omics_type = "proteomics"),
    # `comparison` is an output of run_diff, not an argument to it.
    params = list(method = "limma", analysis_type = "group",
                  group_col = "group", control_group = "G1",
                  case_group = "G2", comparison = "G2_vs_G1")
  ))
  txt <- paste(export_script(proj), collapse = "\n")
  expect_match(txt, 'method        = "limma"', fixed = TRUE)
  expect_match(txt, 'control_group = "G1"', fixed = TRUE)
  expect_false(grepl("comparison =", txt, fixed = TRUE))
})

test_that("arguments follow the function signature order", {
  proj <- fixture_project("raw/x.xlsx")
  proj$bundles <- list(diff = new_analysis_bundle(
    analysis_name = "run_diff",
    input_info = list(omics_type = "proteomics"),
    # Deliberately shuffled relative to run_diff()'s signature.
    params = list(case_group = "G2", method = "limma",
                  group_col = "group", analysis_type = "group")
  ))
  lines <- export_script(proj)
  pos <- function(arg) grep(paste0("^\\s+", arg, "\\s"), lines)[[1L]]
  expect_lt(pos("method"), pos("analysis_type"))
  expect_lt(pos("analysis_type"), pos("group_col"))
  expect_lt(pos("group_col"), pos("case_group"))
})

test_that("bundles are emitted in dependency order", {
  proj <- fixture_project("raw/x.xlsx")
  proj$bundles <- list(
    enrich = new_analysis_bundle("run_enrichment",
               input_info = list(omics_type = "proteomics"),
               params = list(type = "ora", database = "hallmark")),
    diff = new_analysis_bundle("run_diff",
             input_info = list(omics_type = "proteomics"),
             params = list(method = "limma", analysis_type = "group")),
    qc = new_analysis_bundle("run_qc",
           input_info = list(omics_type = "proteomics"),
           params = list(impute_method = "half_min"))
  )
  lines <- export_script(proj)
  # Listed out of order above; enrichment consumes `diff`, so the script
  # has to define it first or it will not run.
  expect_lt(grep("^qc <- run_qc\\($", lines)[[1L]],
            grep("^diff <- run_diff\\($", lines)[[1L]])
  expect_lt(grep("^diff <- run_diff\\($", lines)[[1L]],
            grep("^enrich <- run_enrichment\\($", lines)[[1L]])
})

test_that("a non-literal parameter becomes a note, not broken code", {
  proj <- fixture_project("raw/x.xlsx")
  proj$bundles <- list(diff = new_analysis_bundle(
    analysis_name = "run_diff",
    input_info = list(omics_type = "proteomics"),
    params = list(method = "limma", covariates = quote(f(x)))
  ))
  lines <- export_script(proj)
  # The note names the argument, so it does appear in the text -- what
  # must not appear is a `covariates = <something unrunnable>` line in
  # the call itself.
  expect_true(any(grepl("not a literal", lines, fixed = TRUE)))
  expect_true(any(grepl("covariates", lines, fixed = TRUE)))
  expect_false(any(grepl("^\\s+covariates\\s+=", lines)))
})

test_that("export_script writes to disk when given a path", {
  path <- tempfile(fileext = ".R"); on.exit(unlink(path), add = TRUE)
  lines <- export_script(fixture_project("raw/x.xlsx"), path)
  expect_true(file.exists(path))
  expect_equal(readLines(path), lines)
})

test_that("plotting calls are opt-in", {
  proj <- fixture_project("raw/x.xlsx")
  proj$bundles <- list(diff = new_analysis_bundle("run_diff",
    input_info = list(omics_type = "proteomics"),
    params = list(method = "limma")))
  expect_false(any(grepl("plot_volcano", export_script(proj))))
  expect_true(any(grepl("plot_volcano",
                        export_script(proj, include_plots = TRUE))))
})

# ---- round trip -------------------------------------------------------

test_that("an exported script reproduces the analysis it describes", {
  skip_on_cran()
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("readxl")
  skip_if_not_installed("limma")
  rscript <- file.path(R.home("bin"), "Rscript")
  skip_if_not(file.exists(rscript), "Rscript not found")

  work <- tempfile("roundtrip")
  dir.create(file.path(work, "raw"), recursive = TRUE)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)

  data_file <- file.path(work, "raw", "fixture.xlsx")
  write_tiny_workbook(data_file)

  parsed <- read_omics(data_file, omics_type = "proteomics",
                       assay_type = "intensity")
  inp <- parsed$input
  inp$source_path <- data_file

  # `method = "auto"` on purpose: the script must record the engine that
  # actually ran, or a reader cannot tell what produced the numbers.
  original <- run_diff(inp, method = "auto", analysis_type = "group",
                       group_col = "group", control_group = "G1",
                       case_group = "G2")

  proj <- omics_project("roundtrip", experiments = list(proteomics = inp))
  proj$bundles <- list(diff = original)

  script <- file.path(work, "reproduce.R")
  export_script(proj, script)

  expect_false(any(grepl('method        = "auto"', readLines(script),
                         fixed = TRUE)))
  expect_true(any(grepl('method        = "limma"', readLines(script),
                        fixed = TRUE)))

  out_rds <- file.path(work, "out.rds")
  writeLines(
    c(readLines(script),
      sprintf("saveRDS(diff$results$diff_result_df, %s)",
              deparse(out_rds))),
    script
  )

  old_wd <- setwd(work)
  on.exit(setwd(old_wd), add = TRUE)
  log <- suppressWarnings(system2(
    rscript, c("--vanilla", shQuote(script)),
    stdout = TRUE, stderr = TRUE,
    env = paste0("R_LIBS=", shQuote(paste(.libPaths(), collapse = .Platform$path.sep)))
  ))

  expect_true(file.exists(out_rds),
              info = paste(utils::tail(log, 30), collapse = "\n"))

  reproduced <- readRDS(out_rds)
  expected <- original$results$diff_result_df

  expect_equal(nrow(reproduced), nrow(expected))
  expect_equal(reproduced$feature_id, expected$feature_id)
  # The whole point: the numbers a reader would check must match.
  expect_equal(reproduced$effect, expected$effect, tolerance = 1e-12)
  expect_equal(reproduced$p_value, expected$p_value, tolerance = 1e-12)
  expect_equal(reproduced$adj_p_value, expected$adj_p_value,
               tolerance = 1e-12)
})
