# The items carried over from the round-4 review: readable plots with
# long labels, one vocabulary for effects and p-values, and the steps a
# long computation reports while it runs.

bl_diff <- function(groups, n_per = 4L, n_feat = 120L) {
  set.seed(3)
  ids <- paste0("S", seq_len(length(groups) * n_per))
  meta <- data.frame(group = rep(groups, each = n_per), row.names = ids)
  m <- matrix(stats::rnorm(n_feat * length(ids), 20, 0.5), n_feat,
              dimnames = list(paste0("G", seq_len(n_feat)), ids))
  for (k in seq_along(groups)[-1L]) {
    m[seq_len(10L * k), meta$group == groups[k]] <- m[seq_len(10L * k), meta$group == groups[k]] + 1.5
  }
  omics_input(m, meta, data.frame(feature_id = rownames(m), feature_symbol = rownames(m)),
              omics_type = "proteomics", assay_type = "normalized_intensity")
}

test_that("long labels wrap at spaces, inside long words, and by display width", {
  expect_identical(wrap_label("Short"), "Short")
  expect_identical(wrap_label(NA_character_), NA_character_)
  w <- wrap_label("Treatment_high_dose_week_12_vs_Control", width = 20)
  expect_identical(gsub("\n", "", w), "Treatment_high_dose_week_12_vs_Control")
  expect_true(all(nchar(strsplit(w, "\n")[[1]]) <= 20))
  # No break character at all: cut by length.
  expect_identical(strsplit(wrap_label(strrep("A", 30), width = 12), "\n")[[1]],
                   c(strrep("A", 12), strrep("A", 12), strrep("A", 6)))
  # A Chinese label is twice as wide per character.
  zh <- wrap_label("处理组甲高剂量第十二周与对照组相比较", width = 12)
  expect_true(all(nchar(strsplit(zh, "\n")[[1]], type = "width") <= 12))
  long <- wrap_label(paste(rep("word", 40), collapse = " "), width = 10, max_lines = 2)
  expect_length(strsplit(long, "\n")[[1]], 2L)
  expect_match(long, "…$")
})

test_that("a comparison label keeps its control however long the case name", {
  lbl <- wrap_comparison("Treatment with a very long descriptive name number 1 at high dose_vs_Control")
  expect_match(lbl, "\nvs Control$")
  expect_identical(wrap_comparison("A_vs_B"), "A vs B")
  expect_identical(comparison_label_lines(c("A_vs_B", "A very long case name indeed, longer still_vs_Control")),
                   c(1L, 3L))
})

test_that("the contrast and overlap plots wrap their comparison names", {
  skip_if_not_installed("limma")
  groups <- c("Control", paste("Treatment with a very long descriptive name", 1:2), "Short")
  d <- run_diff(bl_diff(groups), method = "limma", group_col = "group",
                control_group = "Control", case_group = groups[-1])
  p <- plot_diff_contrasts(d)
  b <- ggplot2::ggplot_build(p)
  labs <- b$layout$panel_params[[1]]$y$get_labels()
  expect_true(any(grepl("\n", labs, fixed = TRUE)))
  expect_true(all(grepl("Control$", labs)))
  expect_identical(p$theme$plot.title.position, "plot")
  expect_identical(p$labels$x, "features, down | up")
  ov <- plot_diff_overlap(d)
  expect_s3_class(ov, "patchwork")
  dots <- ggplot2::ggplot_build(ov[[2]])
  expect_true(all(grepl("\\(\\d+\\)$", dots$layout$panel_params[[1]]$y$get_labels())))
})

test_that("a PCA of up to six groups draws them as shapes as well as hues", {
  inp <- bl_diff(c("Control", "A", "B"))
  q <- run_qc(inp, outlier_method = "none", impute_method = "none")
  p <- plot_qc(q, view = "pca", color_by = "group")
  expect_true("shape" %in% names(p$mapping))
  inp8 <- bl_diff(LETTERS[1:8], n_per = 2L)
  q8 <- run_qc(inp8, outlier_method = "none", impute_method = "none")
  expect_false("shape" %in% names(plot_qc(q8, view = "pca", color_by = "group")$mapping))
  expect_true("shape" %in% names(plot_pca(inp, color_by = "group")$mapping))
})

test_that("effects and p-values are named the same way everywhere", {
  expect_identical(effect_label("log2FC"), "log2FC")
  expect_identical(effect_label("log2FC_per_unit"), "slope")
  expect_identical(effect_label("correlation"), "rho")
  expect_identical(effect_label(character(0)), "log2FC")
  skip_if_not_installed("limma")
  d <- run_diff(bl_diff(c("Control", "A")), method = "limma", group_col = "group",
                control_group = "Control", case_group = "A")
  v <- plot_volcano(d, effect_threshold = 1)
  expect_match(v$labels$caption, "adjusted p < 0.05, |log2FC| >= 1", fixed = TRUE)
  expect_identical(v$labels$y, "-log10(adjusted p)")
  expect_match(plot_diff_contrasts(d, effect_cutoff = 0.5)$labels$subtitle, "|log2FC|", fixed = TRUE)
})

test_that("run_qc reports its steps to whoever listens, and to no one by default", {
  steps <- character(0)
  withr::local_options(omicsCore.progress = function(s) steps <<- c(steps, s))
  inp <- bl_diff(c("Control", "A"))
  inp$expr_mat[1:5, 1] <- NA
  run_qc(inp, outlier_method = "pca", impute_method = "MinDet")
  expect_identical(steps, c("Counting missing values", "Checking for outlier samples",
                            "Imputing missing values"))
  withr::local_options(omicsCore.progress = NULL)
  expect_silent(run_qc(inp, outlier_method = "none", impute_method = "none"))
})

bl_counts <- function() {
  set.seed(11)
  ids <- paste0("S", 1:8)
  meta <- data.frame(group = rep(c("ctrl", "trt"), each = 4), row.names = ids)
  mu <- matrix(rep(stats::rgamma(400, 2, 0.02), 8), 400)
  mu[1:40, 5:8] <- mu[1:40, 5:8] * 4
  m <- matrix(stats::rnbinom(length(mu), mu = mu, size = 10), 400,
              dimnames = list(paste0("ENSG", sprintf("%011d", 1:400)), ids))
  omics_input(m, meta, data.frame(feature_id = rownames(m), feature_symbol = rownames(m)),
              omics_type = "rnaseq", assay_type = "raw_count")
}

test_that("edgeR reports its four steps", {
  skip_if_not_installed("edgeR")
  steps <- character(0)
  withr::local_options(omicsCore.progress = function(s) steps <<- c(steps, s))
  run_diff(bl_counts(), method = "edger", group_col = "group",
           control_group = "ctrl", case_group = "trt")
  expect_identical(steps, c("Normalising library sizes (1 of 4)", "Estimating dispersions (2 of 4)",
                            "Fitting the model (3 of 4)", "Testing (4 of 4)"))
})

test_that("DESeq2 reports its steps from its own messages and still prints none", {
  skip_if_not_installed("DESeq2")
  inp <- bl_counts()
  steps <- character(0)
  withr::local_options(omicsCore.progress = function(s) steps <<- c(steps, s))
  msgs <- character(0)
  d <- withCallingHandlers(
    run_diff(inp, method = "deseq2", group_col = "group",
             control_group = "ctrl", case_group = "trt"),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  # The fit's own messages are turned into steps, not printed.
  expect_false(any(grepl("estimating|fitting model|dispersion", msgs)))
  expect_true(all(c("Estimating size factors (1 of 4)", "Estimating dispersions (2 of 4)",
                    "Fitting the model and testing (3 of 4)") %in% steps))
  # The same fit as the quiet one it replaced.
  cd <- inp$meta_df
  cd$group <- factor(cd$group)
  dds <- suppressMessages(DESeq2::DESeqDataSetFromMatrix(round(inp$expr_mat), cd, ~ group))
  ref <- with_fixed_seed(1L, DESeq2::DESeq(dds, quiet = TRUE))
  ref <- as.data.frame(DESeq2::results(ref, contrast = c("group", "trt", "ctrl")))
  got <- d$results$diff_result_df
  expect_equal(got$p_value[match(rownames(ref), got$feature_id)], ref$pvalue, tolerance = 1e-10)
})

test_that("DESeq2 in parallel, when asked for, gives the serial result", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("BiocParallel")
  skip_on_os("windows")
  inp <- bl_counts()
  serial <- run_diff(inp, method = "deseq2", group_col = "group",
                     control_group = "ctrl", case_group = "trt")
  expect_null(deseq2_bpparam())
  withr::local_options(omicsCore.deseq2_workers = 2L)
  expect_s4_class(deseq2_bpparam(), "BiocParallelParam")
  par <- run_diff(inp, method = "deseq2", group_col = "group",
                  control_group = "ctrl", case_group = "trt")
  expect_equal(par$results$diff_result_df$p_value, serial$results$diff_result_df$p_value,
               tolerance = 1e-8)
})

test_that("a result is matched to the layer it records, then by omics type", {
  inp <- bl_diff(c("A", "B"))
  proj <- omics_project("P", list(batch1 = inp, batch2 = inp))
  b <- new_analysis_bundle("run_qc", input_info = list(omics_type = "proteomics"))
  expect_identical(bundle_layer(proj, b), "batch1")       # legacy: first of its type
  b$input_info$layer <- "batch2"
  expect_identical(bundle_layer(proj, b), "batch2")
  b$input_info$layer <- "gone"
  expect_true(is.na(bundle_layer(proj, b)))
  expect_true(is.na(bundle_layer(proj, new_analysis_bundle("x", input_info = list(omics_type = "rnaseq")))))
  # Enrichment inherits the layer from the differential result it reads.
  skip_if_not_installed("limma")
  d <- run_diff(inp, method = "limma", group_col = "group", control_group = "A", case_group = "B")
  d$input_info$layer <- "batch2"
  expect_identical(resolve_tag(proj, d), "batch2")
})

test_that("an uploaded .rds holding code is refused, data is read", {
  inp <- bl_diff(c("A", "B"))
  ok <- withr::local_tempfile(fileext = ".rds")
  saveRDS(inp, ok)
  expect_false(is.null(read_omics(ok)$input))
  bad <- withr::local_tempfile(fileext = ".rds")
  evil <- inp
  evil$misc <- list(hook = function() stop("ran"))
  saveRDS(evil, bad)
  r <- read_omics(bad)
  expect_null(r$input)
  expect_match(paste(import_report_warnings(r$report), collapse = " "), "a function", fixed = TRUE)
})

test_that("a layer merged from quantification files is scripted with read_quant_files", {
  dir <- withr::local_tempdir()
  raw <- file.path(dir, "raw"); dir.create(raw)
  ids <- paste0("ENSG", sprintf("%011d", 1:30))
  set.seed(5)
  mk <- function(f) {
    reads <- stats::rpois(30, 200); eff <- stats::runif(30, 500, 3000)
    tpm <- reads / eff; tpm <- tpm / sum(tpm) * 1e6
    utils::write.table(data.frame(Name = ids, Length = round(eff) + 150, EffectiveLength = eff,
                                  TPM = tpm, NumReads = reads),
                       file.path(raw, f), sep = "\t", quote = FALSE, row.names = FALSE)
    file.path(raw, f)
  }
  # Archived under store names; uploaded as these.
  paths <- vapply(c("a__1.sf", "b__2.sf", "c__3.sf", "d__4.sf"), mk, character(1))
  shown <- c("ctrl_1.quant.sf", "ctrl_2.quant.sf", "ko_1.quant.sf", "ko_2.quant.sf")
  inp <- read_quant_files(unname(paths), file_names = shown)$input
  inp$quant_source <- list(paths = unname(paths), names = shown)
  proj <- omics_project("q", list(rnaseq = inp))
  out <- withr::local_tempfile(fileext = ".R")
  export_script(proj, out)
  txt <- readLines(out)
  expect_true(any(grepl("read_quant_files(", txt, fixed = TRUE)))
  expect_false(any(grepl("not archived", txt, fixed = TRUE)))
  # The read line runs from the folder the script sits in, and gives
  # back the same samples the app had.
  i <- grep("^[^#]*read_quant_files\\(", txt)[1L]
  j <- i
  while (!grepl("\\$input$", txt[j])) j <- j + 1L
  env <- new.env()
  withr::with_dir(dir, eval(parse(text = txt[i:j]), envir = env))
  got <- get(ls(env)[1], env)
  expect_identical(colnames(got$expr_mat), colnames(inp$expr_mat))
  expect_equal(got$expr_mat, inp$expr_mat)
  # And it survives a QC exclusion.
  sub <- subset_omics(inp, samples = colnames(inp$expr_mat)[-1])
  expect_identical(sub$quant_source, inp$quant_source)
})
