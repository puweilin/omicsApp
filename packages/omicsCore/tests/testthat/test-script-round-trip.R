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

# --- an integration, replayed --------------------------------------------

# Concordance takes two differentials and the second is computed inside
# the integration view, displayed by nothing. It used to exist only in the
# worker's scope, so the integration named an upstream the project did not
# hold, the exported script emitted one differential, and running it
# stopped on "`diff_bundles` is required". A NOTE about the gap does not
# make a script executable.
rt_two_layers <- function(dir) {
  raw <- file.path(dir, "raw")
  dir.create(raw, showWarnings = FALSE, recursive = TRUE)
  syms <- paste0("SYM", 1:40)
  mk <- function(file, seed, shift) {
    set.seed(seed)
    n <- 8; ids <- paste0("s", seq_len(n))
    m <- matrix(rlnorm(40 * n, 12, 0.5), nrow = 40, dimnames = list(syms, ids))
    m[1:8, 5:8] <- m[1:8, 5:8] * shift
    openxlsx::write.xlsx(list(
      matrix = data.frame(feature_id = syms, m, check.names = FALSE),
      samples = data.frame(sample_id = ids,
                           group = rep(c("A", "B"), each = 4),
                           stringsAsFactors = FALSE)), file)
    file
  }
  p1 <- mk(file.path(raw, "prot.xlsx"), 21, 4)
  p2 <- mk(file.path(raw, "rna.xlsx"), 22, 2)

  layer <- function(path) {
    inp <- suppressWarnings(read_omics(path, omics_type = "proteomics",
                                       assay_type = "raw_intensity")$input)
    inp$source_path <- path
    suppressMessages(normalize_omics(inp, method = "log2"))
  }
  a <- layer(p1); b <- layer(p2)
  proj <- omics_project("two layer", list(prot = a, rna = b))

  run <- function(tag) {
    run_diff(effective_input(proj, tag)$input, method = "limma",
             analysis_type = "group", group_col = "group",
             control_group = "A", case_group = "B")
  }
  d_a <- run("prot"); d_b <- run("rna")

  # As the application now files them: both differentials in the registry,
  # the integration naming both as upstream.
  bind <- function(b, tag, id) {
    b$input_info$layer_tag <- tag
    b$input_info$input_version <- effective_input(proj, tag)$version
    b$provenance_binding <- new_provenance_binding(
      analysis_id = id, project_id = proj$metadata$project_id,
      inputs = stats::setNames(list(effective_input(proj, tag)$version), tag),
      params_signature = params_signature(b$params))
    b
  }
  d_a <- bind(d_a, "prot", "ana-a")
  d_b <- bind(d_b, "rna", "ana-b")

  integ <- run_integration(project = proj, method = "concordance",
                           experiments = c("prot", "rna"),
                           diff_bundles = list(prot = d_a, rna = d_b))
  integ$provenance_binding <- new_provenance_binding(
    analysis_id = "ana-int", project_id = proj$metadata$project_id,
    inputs = list(prot = effective_input(proj, "prot")$version,
                  rna = effective_input(proj, "rna")$version),
    upstream_ids = c("ana-a", "ana-b"),
    params_signature = params_signature(integ$params))

  proj <- register_analysis(proj, d_a)
  proj <- register_analysis(proj, d_b)
  proj <- register_analysis(proj, integ)
  proj$bundles <- list(diff = d_a, integration = integ)
  list(project = proj, integration = integ, dir = dir)
}

test_that("a two-layer integration replays to the same numbers", {
  skip_unless_replayable()
  dir <- withr::local_tempdir()
  p <- rt_two_layers(dir)

  out <- withr::local_tempdir()
  dir.create(file.path(out, "raw"), recursive = TRUE, showWarnings = FALSE)
  file.copy(list.files(file.path(dir, "raw"), full.names = TRUE),
            file.path(out, "raw"))
  script <- file.path(out, "replay.R")
  export_script(p$project, script)
  txt <- paste(readLines(script), collapse = "\n")
  # Both layers, and the bundles the method needs.
  expect_match(txt, "diff_bundles")
  expect_equal(length(gregexpr("run_diff\\(", txt)[[1L]]), 3L)

  env <- new.env(parent = globalenv())
  withr::with_dir(out, suppressWarnings(suppressMessages(
    sys.source(script, envir = env))))
  replayed <- env$integration
  expect_false(is.null(replayed))

  a <- replayed$results$integration_df
  b <- p$integration$results$integration_df
  i <- match(b$feature_id, a$feature_id)
  expect_false(anyNA(i))
  for (col in c("effect", "effect_a", "effect_b", "p_value", "adj_p_value")) {
    expect_equal(a[[col]][i], b[[col]], tolerance = 1e-8, info = col)
  }
})

test_that("an integration whose upstreams are not registered says so", {
  skip_unless_replayable()
  dir <- withr::local_tempdir()
  p <- rt_two_layers(dir)
  p$project$analyses <- list()          # a project from before they were filed
  script <- withr::local_tempfile(fileext = ".R")
  res <- export_script(p$project, script)
  txt <- paste(readLines(script), collapse = "\n")
  expect_match(paste(res, collapse = "\n"), "not registered")
})
