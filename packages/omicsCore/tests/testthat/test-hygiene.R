# What a call leaves behind.
#
# A Shiny process serves every user from one R session, and a library
# call runs inside whatever script called it. Either way, an analysis
# that changes state it was not asked to change charges the cost to
# someone else: the next session's imputation draws from a stream a
# stranger's GSEA left behind, a figure jitters differently on every
# render, the report and the script that reproduces it disagree.
#
# Found this way: GSEA advanced the caller's random stream on every run
# and planted a seed in a session that had none, DESeq2 planted one
# too, and the feature-expression plot drew unseeded jitter on each
# render. The report's figure directory was written beside the
# template inside the installed package.

skip_if_not_installed("limma")

prot <- realistic_input("proteomics")
prot$meta_df$age <- seq_len(ncol(prot$expr_mat)) * 3 + 20
rna <- realistic_input("rnaseq")
two_layers <- omics_project("hygiene", experiments = list(proteomics = prot, rnaseq = rna))

group_diff <- function(input, method) {
  suppressMessages(run_diff(input, method = method, analysis_type = "group",
                            group_col = "group", control_group = "G1",
                            case_group = "G2"))
}

has <- function(pkg) requireNamespace(pkg, quietly = TRUE)

# The first use of an engine loads its namespace, and a namespace load
# is allowed to set options, environment variables and a seed of its
# own -- HDF5Array and magick do. The contract here is about the call,
# not the load, so every engine is used once before any baseline is
# taken.
warm <- local({
  done <- FALSE
  function() {
    if (done) return(invisible())
    done <<- TRUE
    suppressWarnings(suppressMessages({
      db <- group_diff(prot, "limma")
      if (has("DESeq2")) group_diff(rna, "deseq2")
      if (has("edgeR")) group_diff(rna, "edger")
      if (has("clusterProfiler")) {
        run_enrichment(db, type = "ora", database = "hallmark")
        run_enrichment(db, type = "gsea", database = "hallmark")
      }
      if (has("GSVA")) run_gsva(prot, database = "hallmark")
    }))
    invisible()
  }
})

seed_state <- function() {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  }
}

expect_stream_untouched <- function(expr, label) {
  set.seed(20260906L)
  before <- seed_state()
  force(expr)
  expect_identical(seed_state(), before, label = label)
}

expect_no_seed_planted <- function(expr, label) {
  suppressWarnings(rm(".Random.seed", envir = globalenv()))
  force(expr)
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE),
               label = label)
}

test_that("a seeded stream is where the caller left it", {
  warm()
  withr::defer(set.seed(1))
  db <- group_diff(prot, "limma")
  expect_stream_untouched(run_qc(prot), "run_qc")
  expect_stream_untouched(group_diff(prot, "limma"), "limma")
  expect_stream_untouched(group_diff(prot, "ttest"), "ttest")
  expect_stream_untouched(run_diff_continuous(prot, method = "lm", continuous_col = "age"), "lm")
  if (has("DESeq2")) expect_stream_untouched(group_diff(rna, "deseq2"), "deseq2")
  if (has("edgeR")) expect_stream_untouched(group_diff(rna, "edger"), "edger")
  if (has("clusterProfiler")) {
    expect_stream_untouched(run_enrichment(db, type = "ora", database = "hallmark"), "ora")
    expect_stream_untouched(
      suppressWarnings(run_enrichment(db, type = "gsea", database = "hallmark")), "gsea")
  }
  if (has("GSVA")) expect_stream_untouched(run_gsva(prot, database = "hallmark"), "gsva")
  if (has("imputeLCMD")) {
    expect_stream_untouched(run_qc(prot, impute_method = "MinProb"), "MinProb imputation")
  }
  if (has("edgeR")) {
    dbe <- group_diff(rna, "edger")
    expect_stream_untouched(
      run_integration(two_layers, method = "concordance",
                      experiments = c("proteomics", "rnaseq"),
                      diff_bundles = list(proteomics = db, rnaseq = dbe)),
      "integration")
  }
  expect_stream_untouched(ggplot2::ggplot_build(plot_volcano(db)), "volcano")
  expect_stream_untouched(ggplot2::ggplot_build(plot_pca(prot, color_by = "group")), "pca")
  expect_stream_untouched(
    ggplot2::ggplot_build(plot_feature_expression(
      prot, features = rownames(prot$expr_mat)[1:2], group_by = "group")),
    "feature expression")
})

test_that("no seed is planted in a session that had none", {
  warm()
  withr::defer(set.seed(1))
  db <- group_diff(prot, "limma")
  expect_no_seed_planted(group_diff(prot, "limma"), "limma")
  expect_no_seed_planted(group_diff(prot, "ttest"), "ttest")
  expect_no_seed_planted(run_qc(prot), "run_qc")
  if (has("DESeq2")) expect_no_seed_planted(group_diff(rna, "deseq2"), "deseq2")
  if (has("clusterProfiler")) {
    expect_no_seed_planted(run_enrichment(db, type = "ora", database = "hallmark"), "ora")
    expect_no_seed_planted(
      suppressWarnings(run_enrichment(db, type = "gsea", database = "hallmark")), "gsea")
  }
  expect_no_seed_planted(
    ggplot2::ggplot_build(plot_feature_expression(
      prot, features = rownames(prot$expr_mat)[1:2], group_by = "group")),
    "feature expression")
})

test_that("the engines answer the same table twice, whatever the stream", {
  warm()
  withr::defer(set.seed(1))
  db <- group_diff(prot, "limma")
  set.seed(1); a <- group_diff(prot, "limma")$results$diff_result_df
  set.seed(2); b <- group_diff(prot, "limma")$results$diff_result_df
  expect_identical(a, b)
  if (has("DESeq2")) {
    set.seed(1); a <- group_diff(rna, "deseq2")$results$diff_result_df
    set.seed(2); invisible(stats::runif(100)); b <- group_diff(rna, "deseq2")$results$diff_result_df
    expect_identical(a$p_value, b$p_value)
  }
  if (has("clusterProfiler")) {
    set.seed(1); a <- suppressWarnings(run_enrichment(db, type = "gsea", database = "hallmark"))
    set.seed(2); b <- suppressWarnings(run_enrichment(db, type = "gsea", database = "hallmark"))
    expect_identical(a$results$enrich_result_df$p_value, b$results$enrich_result_df$p_value)
  }
  if (has("imputeLCMD")) {
    set.seed(1); a <- run_qc(prot, impute_method = "MinProb")
    set.seed(2); b <- run_qc(prot, impute_method = "MinProb")
    expect_identical(a$results$imputed_matrix, b$results$imputed_matrix)
  }
})

test_that("a jittered figure is the same on every render", {
  p <- plot_feature_expression(prot, features = rownames(prot$expr_mat)[1:3],
                               group_by = "group")
  x1 <- ggplot2::ggplot_build(p)$data
  set.seed(99)
  x2 <- ggplot2::ggplot_build(p)$data
  expect_identical(x1, x2)
})

test_that("nothing else is left behind by a full pipeline", {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("readxl")
  warm()
  withr::defer(set.seed(1))
  mine <- tempfile("hygiene-")
  dir.create(mine)
  withr::defer(unlink(mine, recursive = TRUE))
  xlsx <- file.path(mine, "prot.xlsx")
  openxlsx::write.xlsx(list(
    expression = data.frame(feature_id = rownames(prot$expr_mat), prot$expr_mat,
                            check.names = FALSE),
    metadata = data.frame(sample_id = rownames(prot$meta_df), prot$meta_df),
    features = prot$feature_df), xlsx)

  can_report <- has("rmarkdown") && rmarkdown::pandoc_available()

  before <- list(
    conns = nrow(showConnections()), devs = length(grDevices::dev.list()),
    wd = getwd(), locale = Sys.getlocale(), sink = sink.number(),
    globals = ls(globalenv(), all.names = TRUE),
    opts = options(), env = Sys.getenv(),
    tmp = list.files(tempdir(), recursive = TRUE, all.files = TRUE)
  )

  parsed <- read_omics(xlsx, omics_type = "proteomics", assay_type = "normalized_intensity")
  input <- parsed$input
  qc <- run_qc(input)
  db <- group_diff(input, "limma")
  proj <- omics_project("hygiene", experiments = list(proteomics = input))
  proj$bundles <- list(qc = qc, diff = db)
  if (has("clusterProfiler")) {
    proj$bundles$enrich <- run_enrichment(db, type = "ora", database = "hallmark")
    ggplot2::ggplot_build(plot_enrichment(proj$bundles$enrich))
  }
  ggplot2::ggplot_build(plot_qc(qc, view = "pca"))
  ggplot2::ggplot_build(plot_volcano(db))
  export_bundle(db, file.path(mine, "bundle"), formats = c("tsv", "pdf"))
  if (can_report) export_report(proj, file.path(mine, "report.html"))
  export_script(proj, file.path(mine, "script.R"))
  save_project(proj, file.path(mine, "p.omp"))
  load_project(file.path(mine, "p.omp"))

  after <- list(
    conns = nrow(showConnections()), devs = length(grDevices::dev.list()),
    wd = getwd(), locale = Sys.getlocale(), sink = sink.number(),
    globals = ls(globalenv(), all.names = TRUE),
    opts = options(), env = Sys.getenv(),
    tmp = list.files(tempdir(), recursive = TRUE, all.files = TRUE)
  )
  expect_identical(after$conns, before$conns, label = "open connections")
  expect_identical(after$devs, before$devs, label = "graphics devices")
  expect_identical(after$wd, before$wd)
  expect_identical(after$locale, before$locale)
  expect_identical(after$sink, before$sink)
  # clusterProfiler keeps an annotation cache in the global environment;
  # that is its design, not a leak of ours.
  expect_setequal(setdiff(after$globals, c(before$globals, ".Anno_clusterProfiler_Env", ".Random.seed")),
                  character(0))
  # Options and variables that existed before must read the same after.
  # New ones are what a namespace load adds, and are allowed.
  shared_opts <- intersect(names(before$opts), names(after$opts))
  changed <- shared_opts[!vapply(shared_opts, function(n) identical(before$opts[[n]], after$opts[[n]]), logical(1))]
  expect_identical(changed, character(0), label = "changed options")
  shared_env <- intersect(names(before$env), names(after$env))
  changed_env <- shared_env[!vapply(shared_env, function(n) identical(unname(before$env[n]), unname(after$env[n])), logical(1))]
  expect_identical(changed_env, character(0), label = "changed environment variables")
  # The only new files in tempdir() are the ones this test asked for,
  # under its own directory. A report used to leave nothing here but
  # rendered its intermediates beside the template instead.
  new_tmp <- setdiff(after$tmp, before$tmp)
  foreign <- new_tmp[!startsWith(new_tmp, basename(mine))]
  expect_identical(foreign, character(0), label = "temp files left by the pipeline")
})

test_that("the report renders nothing beside its template", {
  skip_if_not_installed("rmarkdown")
  skip_if_not(rmarkdown::pandoc_available(), "pandoc is not available")
  template_dir <- dirname(system.file("rmd", "default-report.Rmd", package = "omicsCore"))
  before <- list.files(template_dir, recursive = TRUE, all.files = TRUE)
  proj <- omics_project("p", experiments = list(proteomics = prot))
  proj$bundles <- list(diff = group_diff(prot, "limma"))
  out <- tempfile(fileext = ".html")
  withr::defer(unlink(out))
  export_report(proj, out)
  expect_true(file.exists(out))
  expect_identical(list.files(template_dir, recursive = TRUE, all.files = TRUE), before)
})

test_that("package code seeds the stream only where it puts it back", {
  r_dir <- if (dir.exists(file.path("..", "..", "R"))) file.path("..", "..", "R") else NULL
  skip_if(is.null(r_dir), "package source is not beside the tests")
  files <- list.files(r_dir, pattern = "\\.R$", full.names = TRUE)
  seeders <- Filter(function(f) {
    code <- grep("^\\s*#", readLines(f, warn = FALSE), value = TRUE, invert = TRUE)
    any(grepl("set\\.seed\\(", code))
  }, files)
  # with_fixed_seed() restores the caller's stream on exit; nothing else
  # may call set.seed().
  expect_identical(basename(seeders), "qc-imputation.R")
})
