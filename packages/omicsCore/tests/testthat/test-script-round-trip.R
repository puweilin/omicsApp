# The exported script has to produce the numbers it documents.
#
# It always ran. That was the problem: it re-read the archived file and told
# the analysis the values were already normalized, because the layer knew it
# was on a log scale and not that anything had put it there. So the script
# modelled linear intensities and reported them as log2 fold changes -- the
# same defect as A04, reached from the other end. The replayed effects
# differed from the recorded ones by seven million.
#
# Which means "the script runs" and "the script produces the same result"
# are different assertions, and only the second one is worth making.

skip_unless_replayable <- function() {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("readxl")
  skip_if_not_installed("withr")
  skip_if_not_installed("limma")
}

# A workbook on the linear scale, so normalization is a step that has to
# happen rather than a no-op.
rt_workbook <- function(path, n_features = 40L, n_samples = 8L, seed = 11,
                        with_na = FALSE) {
  set.seed(seed)
  ids <- paste0("s", seq_len(n_samples))
  fids <- paste0("P", sprintf("%03d", seq_len(n_features)))
  m <- matrix(rlnorm(n_features * n_samples, 12, 0.6), nrow = n_features,
              dimnames = list(fids, ids))
  half <- n_samples / 2
  m[1:8, seq_len(half)] <- m[1:8, seq_len(half)] * 4
  if (with_na) m[1:3, 1] <- NA
  mat <- data.frame(feature_id = fids, m, check.names = FALSE)
  samples <- data.frame(sample_id = ids,
                        group = rep(c("A", "B"), each = half),
                        age = seq(30, by = 3, length.out = n_samples),
                        stringsAsFactors = FALSE)
  openxlsx::write.xlsx(list(matrix = mat, samples = samples), path)
  path
}

# Import, process and analyse the way the app does, then hand back both the
# project and the results, so the replay has something to be compared with.
rt_pipeline <- function(dir, normalize = "log2", qc = NULL, seed = 11,
                        n_features = 40L) {
  raw_dir <- file.path(dir, "raw")
  dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)
  path <- rt_workbook(file.path(raw_dir, "layer.xlsx"), seed = seed,
                      n_features = n_features, with_na = !is.null(qc))

  inp <- suppressWarnings(read_omics(path, omics_type = "proteomics",
                                     assay_type = "raw_intensity")$input)
  inp$source_path <- path
  if (!is.null(normalize)) {
    inp <- suppressMessages(normalize_omics(inp, method = normalize))
  }
  proj <- omics_project("round trip", list(proteomics = inp))
  if (!is.null(qc)) {
    b <- do.call(run_qc, c(list(input = inp), qc))
    proj <- set_preprocessing(proj, "proteomics",
                              b$results$cleaned_input, qc)
  }
  eff <- effective_input(proj, "proteomics")$input
  diff <- run_diff(eff, method = "limma", analysis_type = "group",
                   group_col = "group", control_group = "A", case_group = "B")
  proj$bundles <- list(diff = diff)
  list(project = proj, diff = diff, dir = dir)
}

# Run the exported script in a directory of its own, holding only the
# archived file -- which is the situation the script is written for.
rt_replay <- function(p) {
  out <- withr::local_tempdir()
  dir.create(file.path(out, "raw"), recursive = TRUE, showWarnings = FALSE)
  file.copy(file.path(p$dir, "raw", "layer.xlsx"), file.path(out, "raw"))
  script <- file.path(out, "replay.R")
  export_script(p$project, script)
  env <- new.env(parent = globalenv())
  withr::with_dir(out, suppressWarnings(suppressMessages(
    sys.source(script, envir = env))))
  env
}

rt_replayed_diff <- function(env) {
  nm <- grep("^diff", ls(env), value = TRUE)
  if (length(nm) == 0L) return(NULL)
  env[[nm[[1L]]]]
}

# Aligned by feature, never by row position: the replay is free to order
# its output differently, and a positional comparison would either hide a
# real difference or invent one.
expect_same_effects <- function(replayed, original, label) {
  a <- replayed$results$diff_result_df
  b <- original$results$diff_result_df
  expect_setequal(a$feature_id, b$feature_id)
  i <- match(b$feature_id, a$feature_id)
  expect_equal(a$effect[i], b$effect, tolerance = 1e-8, info = label)
  expect_equal(a$statistic[i], b$statistic, tolerance = 1e-8, info = label)
  expect_equal(a$p_value[i], b$p_value, tolerance = 1e-8, info = label)
  expect_equal(a$adj_p_value[i], b$adj_p_value, tolerance = 1e-8, info = label)
}

test_that("the script emits the normalization it performed", {
  skip_unless_replayable()
  dir <- withr::local_tempdir()
  p <- rt_pipeline(dir, normalize = "log2")
  script <- withr::local_tempfile(fileext = ".R")
  export_script(p$project, script)
  txt <- paste(readLines(script), collapse = "\n")
  expect_match(txt, "normalize_omics")
  # And re-reads the file as what it is, not as what it became.
  expect_match(txt, 'assay_type\\s*=\\s*"raw_intensity"')
})

test_that("a log2-normalized layer replays to the same numbers", {
  skip_unless_replayable()
  dir <- withr::local_tempdir()
  p <- rt_pipeline(dir, normalize = "log2")
  env <- rt_replay(p)
  replayed <- rt_replayed_diff(env)
  expect_false(is.null(replayed))
  expect_same_effects(replayed, p$diff, "log2")
})

test_that("a vsn-normalized layer replays to the same numbers", {
  skip_unless_replayable()
  skip_if_not_installed("vsn")
  dir <- withr::local_tempdir()
  # vsn refuses to estimate from a handful of rows, and rightly.
  p <- rt_pipeline(dir, normalize = "vsn", n_features = 300L)
  env <- rt_replay(p)
  expect_same_effects(rt_replayed_diff(env), p$diff, "vsn")
})

test_that("an applied QC step replays too", {
  skip_unless_replayable()
  dir <- withr::local_tempdir()
  p <- rt_pipeline(dir, normalize = "log2",
                   qc = list(missing_threshold = 0.5, impute_method = "min"))
  script <- withr::local_tempfile(fileext = ".R")
  export_script(p$project, script)
  txt <- paste(readLines(script), collapse = "\n")
  expect_match(txt, "run_qc")
  expect_match(txt, "cleaned_input")

  env <- rt_replay(p)
  replayed <- rt_replayed_diff(env)
  expect_same_effects(replayed, p$diff, "qc")
  # And the layer it modelled is the cleaned one, not the import.
  expect_equal(nrow(replayed$results$diff_result_df),
               nrow(p$diff$results$diff_result_df))
})

test_that("the QC arguments that changed the values are all emitted", {
  skip_unless_replayable()
  dir <- withr::local_tempdir()
  # `shift` and `scale` reach impute_matrix() through `...` and were
  # recorded nowhere, so the same recorded parameters gave a different
  # answer on replay.
  p <- rt_pipeline(dir, normalize = "log2",
                   qc = list(missing_threshold = 0.5, impute_method = "man",
                             shift = 1.6, scale = 0.4))
  script <- withr::local_tempfile(fileext = ".R")
  export_script(p$project, script)
  txt <- paste(readLines(script), collapse = "\n")
  expect_match(txt, "shift")
  expect_match(txt, "scale")
  expect_same_effects(rt_replayed_diff(rt_replay(p)), p$diff, "man imputation")
})

test_that("a layer with nothing recorded stops instead of guessing", {
  skip_unless_replayable()
  dir <- withr::local_tempdir()
  p <- rt_pipeline(dir, normalize = "log2")
  # A project from before any of this was recorded.
  p$project$experiments$proteomics$provenance <- NULL

  script <- withr::local_tempfile(fileext = ".R")
  res <- export_script(p$project, script)
  txt <- paste(readLines(script), collapse = "\n")
  # Emitting the import alone would look complete and produce different
  # numbers, because the missing step is known to change them.
  expect_match(txt, "stop\\(")
  expect_match(txt, "no recorded processing")
  expect_match(paste(attr(res, "notes") %||% "", collapse = " "),
               "predates provenance|", fixed = FALSE)
})
