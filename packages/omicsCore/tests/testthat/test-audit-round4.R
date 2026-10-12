# Defects found in the round-4 accuracy, UX and UI audits.

r4_input <- function(m, meta, type = "proteomics", assay = "normalized_intensity") {
  omics_input(m, meta, data.frame(feature_id = rownames(m), feature_symbol = rownames(m)),
              omics_type = type, assay_type = assay)
}

test_that("paired limma fits the pairs as a fixed block, as limma recommends", {
  skip_if_not_installed("limma")
  set.seed(1)
  n <- 4
  m <- matrix(stats::rnorm(300 * 2 * n, 10), 300) +
    matrix(rep(stats::rnorm(300 * n, 0, 1), 2), 300)
  rownames(m) <- paste0("P", 1:300)
  colnames(m) <- paste0("s", 1:(2 * n))
  meta <- data.frame(group = rep(c("A", "B"), each = n), pair = rep(paste0("d", 1:n), 2),
                     row.names = colnames(m))
  b <- run_diff(r4_input(m, meta), method = "limma", group_col = "group",
                control_group = "A", case_group = "B", paired_col = "pair")
  ref <- limma::topTable(limma::eBayes(limma::lmFit(
    m, stats::model.matrix(~ factor(pair) + factor(group), meta))),
    coef = n + 1, number = Inf, sort.by = "none")
  df <- b$results$diff_result_df
  expect_equal(df$p_value, ref$P.Value, tolerance = 1e-10)
  expect_equal(df$effect, ref$logFC, tolerance = 1e-10)

  # A covariate that is constant within pairs can only be estimated with
  # pairs as a random effect; that fit remains, and says so.
  meta$sex <- rep(c("F", "M"), length.out = n)[as.integer(factor(meta$pair))]
  b2 <- run_diff(r4_input(m, meta), method = "limma", group_col = "group",
                 control_group = "A", case_group = "B", paired_col = "pair",
                 covariates = "sex")
  expect_match(b2$warnings, "random effect", all = FALSE)
})

test_that("a custom contrast may leave groups out", {
  set.seed(2)
  m <- matrix(stats::rnorm(40 * 9), 40, dimnames = list(paste0("P", 1:40), paste0("s", 1:9)))
  meta <- data.frame(group = rep(c("Control", "TreatA", "TreatB"), 3), row.names = colnames(m))
  inp <- r4_input(m, meta)
  b <- run_diff(inp, method = "limma", group_col = "group", contrasts = "TreatA - Control")
  expect_identical(diff_comparisons(b), "TreatA_vs_Control")
  b2 <- run_diff(inp, method = "limma", group_col = "group", contrasts = "`TreatB` - `TreatA`")
  expect_identical(diff_comparisons(b2), "TreatB_vs_TreatA")
})

test_that("the limma global test works with two groups", {
  set.seed(3)
  m <- matrix(stats::rnorm(30 * 6), 30, dimnames = list(paste0("P", 1:30), paste0("s", 1:6)))
  meta <- data.frame(group = rep(c("A", "B"), 3), row.names = colnames(m))
  b <- run_diff(r4_input(m, meta), method = "limma", analysis_type = "anova", group_col = "group")
  two <- run_diff(r4_input(m, meta), method = "limma", group_col = "group",
                  control_group = "A", case_group = "B")
  expect_equal(b$results$diff_result_df$p_value, two$results$diff_result_df$p_value,
               tolerance = 1e-10)
})

test_that("the paired t-test effect is the mean difference over complete pairs", {
  m <- rbind(X1 = c(1, 2, 3, 10, 1.6, 2.5, 3.6, NA))
  colnames(m) <- paste0("s", 1:8)
  meta <- data.frame(group = rep(c("A", "B"), each = 4), pair = rep(paste0("p", 1:4), 2),
                     row.names = colnames(m))
  inp <- r4_input(rbind(m, X2 = m[1, ] + 1), meta)
  b <- suppressWarnings(run_diff(inp, method = "ttest", group_col = "group",
                                 control_group = "A", case_group = "B", paired_col = "pair"))
  df <- b$results$diff_result_df
  expect_equal(df$effect[1], mean(c(1.6, 2.5, 3.6) - c(1, 2, 3)))
  expect_identical(sign(df$effect[1]), sign(df$statistic[1]))
})

test_that("GSVA leaves out features with missing values instead of returning NA pathways", {
  skip_if_not_installed("GSVA")
  set.seed(4)
  genes <- paste0("GENE", 1:40)
  m <- matrix(stats::rnorm(40 * 8, 20), 40, dimnames = list(paste0("F", 1:40), paste0("s", 1:8)))
  m[1, 2] <- NA
  inp <- omics_input(m, data.frame(group = rep(c("A", "B"), 4), row.names = colnames(m)),
                     data.frame(feature_id = rownames(m), feature_symbol = genes),
                     omics_type = "proteomics", assay_type = "normalized_intensity")
  g <- suppressWarnings(run_gsva(inp, gene_sets = list(SET1 = genes[1:20], SET2 = genes[21:40]),
                                 min_size = 10))
  expect_false(anyNA(g$results$gsva_matrix))
  expect_match(g$warnings, "missing values", all = FALSE)
})

test_that("enrichment refuses a direction on a result that has none", {
  set.seed(5)
  m <- matrix(stats::rnorm(30 * 9), 30, dimnames = list(paste0("P", 1:30), paste0("s", 1:9)))
  meta <- data.frame(group = rep(c("A", "B", "C"), 3), row.names = colnames(m))
  an <- run_diff(r4_input(m, meta), method = "limma", analysis_type = "anova", group_col = "group")
  expect_error(run_enrichment(an, type = "gsea"), "no direction")
  expect_error(run_enrichment(an, type = "ora", direction = "up"), "no direction")
})

test_that("GSEA ranks by the statistic's own sign", {
  df <- data.frame(feature_symbol = c("A", "B"), effect = c(-0.2, 0.3),
                   statistic = c(4, -5), statistic_type = "t", p_value = c(0.01, 0.001))
  r <- gsea_rank_vector(df, "feature_symbol")
  expect_equal(unname(r[c("A", "B")]), c(4, -5))
})

test_that("raw counts become voom-style log-CPM, which barely shrinks low counts", {
  m <- matrix(c(10L, 10L, 20L, 20L), 1, dimnames = list("g", paste0("s", 1:4)))
  m <- rbind(m, matrix(1000000L, 1, 4, dimnames = list("big", NULL)))
  inp <- omics_input(m, data.frame(group = c("A", "A", "B", "B"), row.names = colnames(m)),
                     data.frame(feature_id = rownames(m)), omics_type = "rnaseq",
                     assay_type = "raw_count")
  sc <- prepare_diff_scale(inp, "limma")$input$expr_mat
  lib <- colSums(m)
  # (TMM moves the library size by a hair; the old log2(CPM + 0.5) was 0.06 lower.)
  expect_equal(unname(sc["g", 1]), log2((10 + 0.5) / (lib[[1]] + 1) * 1e6), tolerance = 1e-4)
})

test_that("edgeR length offsets keep TMM normalisation", {
  skip_if_not_installed("edgeR")
  set.seed(6)
  g <- 300
  mu <- matrix(exp(stats::rnorm(g, log(200), 1)), g, 6)
  up <- 1:30
  mu[up, 4:6] <- mu[up, 4:6] * 8
  m <- matrix(stats::rnbinom(g * 6, mu = mu, size = 20), g,
              dimnames = list(paste0("G", 1:g), paste0("s", 1:6)))
  inp <- omics_input(m, data.frame(group = rep(c("A", "B"), each = 3), row.names = colnames(m)),
                     data.frame(feature_id = rownames(m)), omics_type = "rnaseq",
                     assay_type = "raw_count")
  inp$misc <- list(tximport = list(length = matrix(1000, g, 6, dimnames = dimnames(m)),
                                   counts_from_abundance = "no"))
  b <- run_diff(inp, method = "edger", group_col = "group", control_group = "A", case_group = "B")
  df <- b$results$diff_result_df
  unchanged <- df[-up, ]
  expect_lt(abs(stats::median(unchanged$effect, na.rm = TRUE)), 0.25)
})

test_that("BAM-tool suffixes come off sample names", {
  expect_identical(strip_vendor_decoration(c("/data/ctrl_1.sorted.bam",
                                             "x/ko_2Aligned.sortedByCoord.out.bam",
                                             "wt_3.dedup.bam")),
                   c("ctrl_1", "ko_2", "wt_3"))
})

test_that("a samples-in-rows table keeps its text columns as sample information", {
  skip_if_not_installed("writexl")
  set.seed(7)
  prot <- c("IL6", "TNF", "CXCL8", "IL10", "CCL2", "IL1B", "IFNG", "VEGFA", "MMP9", "CD40")
  samp <- sprintf("S%02d", 1:16)
  npx <- data.frame(SampleID = samp, Disease = rep(c("ctrl", "case"), 8),
                    matrix(round(stats::rnorm(16 * 10, 5), 3), 16, dimnames = list(NULL, prot)),
                    check.names = FALSE)
  path <- withr::local_tempfile(fileext = ".xlsx")
  writexl::write_xlsx(list(npx = npx), path)
  r <- read_omics(path, omics_type = "proteomics", assay_type = "normalized_intensity")
  expect_identical(colnames(r$input$expr_mat), samp)
  expect_identical(r$input$meta_df$Disease, npx$Disease)
})

test_that("a sample sheet with unrecognised headings is found by its contents", {
  skip_if_not_installed("writexl")
  set.seed(8)
  samp <- paste0("S", 1:6)
  mat <- data.frame(Gene = paste0("G", 1:30),
                    matrix(round(stats::rnorm(180, 20), 2), 30, dimnames = list(NULL, samp)),
                    check.names = FALSE)
  sheet <- data.frame(a = samp, b = rep(c("x", "y"), 3))
  names(sheet) <- c("样本", "分组")
  path <- withr::local_tempfile(fileext = ".xlsx")
  writexl::write_xlsx(list(expr = mat, info = sheet), path)
  r <- read_omics(path, omics_type = "proteomics", assay_type = "normalized_intensity")
  expect_true("分组" %in% names(r$input$meta_df))
})

test_that("a separate sample sheet is matched to a matrix file", {
  set.seed(9)
  samp <- paste0("S", 1:6)
  mat <- data.frame(Gene = paste0("G", 1:30),
                    matrix(round(stats::rnorm(180, 20), 2), 30, dimnames = list(NULL, samp)),
                    check.names = FALSE)
  mp <- withr::local_tempfile(fileext = ".csv")
  sp <- withr::local_tempfile(fileext = ".csv")
  utils::write.csv(mat, mp, row.names = FALSE)
  utils::write.csv(data.frame(sample = rev(samp), group = rep(c("B", "A"), each = 3)), sp,
                   row.names = FALSE)
  r <- read_omics(mp, omics_type = "proteomics", assay_type = "normalized_intensity",
                  sample_sheet = sp)
  expect_identical(r$input$meta_df$group, rep(c("A", "B"), each = 3))
  expect_identical(r$input$sample_sheet_path, sp)
  # The script reads it the same way.
  r$input$source_path <- mp
  lines <- export_script(omics_project("p", list(proteomics = r$input)))
  expect_true(any(grepl("sample_sheet", lines, fixed = TRUE)))
})

test_that("an ordinary matrix with an annotation column is not flagged as a guessed orientation", {
  skip_if_not_installed("writexl")
  set.seed(10)
  samp <- paste0("S", 1:8)
  mat <- data.frame(Protein = paste0("P", 100 + 1:400), Genes = paste0("GENE", 1:400),
                    matrix(round(stats::rnorm(3200, 20), 2), 400, dimnames = list(NULL, samp)),
                    check.names = FALSE)
  path <- withr::local_tempfile(fileext = ".xlsx")
  writexl::write_xlsx(list(x = mat), path)
  r <- read_omics(path, omics_type = "proteomics", assay_type = "normalized_intensity")
  expect_false(any(grepl("was a guess", r$report$warnings)))
})

test_that("an exported script repeats a QC exclusion", {
  inp <- realistic_input(n_per_group = 4L)
  inp$source_path <- "upload.xlsx"
  sub <- subset_omics(inp, samples = colnames(inp$expr_mat)[-1])
  sub$excluded_samples <- colnames(inp$expr_mat)[1]
  lines <- export_script(omics_project("p", list(proteomics = sub)))
  expect_true(any(grepl("subset_omics(", lines, fixed = TRUE)))
  expect_silent(parse(text = lines))
})

test_that("the limma trend fit survives an all-missing feature, and its extras match lm()", {
  set.seed(12)
  ids <- paste0("s", 1:12)
  m <- matrix(stats::rnorm(40 * 12, 10), 40, dimnames = list(paste0("P", 1:40), ids))
  m[1, ] <- NA
  m[2, 3] <- NA
  meta <- data.frame(age = seq(20, 75, by = 5), row.names = ids)
  b <- run_diff(r4_input(m, meta), method = "limma", analysis_type = "continuous",
                continuous_col = "age")
  df <- b$results$diff_result_df
  expect_true(is.na(df$p_value[1]))
  y <- m[2, ]
  ref <- summary(stats::lm(y ~ meta$age))$adj.r.squared
  expect_equal(df$model_fit[2], ref, tolerance = 1e-10)
})

test_that("the batched lm and t-test backends match lm() and t.test()", {
  set.seed(13)
  ids <- paste0("s", 1:10)
  m <- matrix(stats::rnorm(20 * 10, 5), 20, dimnames = list(paste0("P", 1:20), ids))
  m[3, c(1, 6)] <- NA
  meta <- data.frame(group = rep(c("A", "B"), each = 5), sex = rep(c("F", "M"), 5),
                     row.names = ids)
  inp <- r4_input(m, meta)
  lm_b <- run_diff(inp, method = "lm", group_col = "group", control_group = "A",
                   case_group = "B", covariates = "sex")$results$diff_result_df
  y <- m[3, ]
  ref <- summary(stats::lm(y ~ factor(meta$group) + meta$sex))$coefficients[2, ]
  expect_equal(lm_b$effect[3], unname(ref[1]), tolerance = 1e-10)
  expect_equal(lm_b$p_value[3], unname(ref[4]), tolerance = 1e-10)
  tt_b <- run_diff(inp, method = "ttest", group_col = "group", control_group = "A",
                   case_group = "B")$results$diff_result_df
  tref <- stats::t.test(m[3, 6:10], m[3, 1:5])
  expect_equal(tt_b$p_value[3], tref$p.value, tolerance = 1e-10)
})

test_that("count-model objects are not kept in the bundle", {
  skip_if_not_installed("edgeR")
  set.seed(14)
  m <- matrix(stats::rnbinom(200 * 6, mu = 100, size = 10), 200,
              dimnames = list(paste0("G", 1:200), paste0("s", 1:6)))
  inp <- omics_input(m, data.frame(group = rep(c("A", "B"), each = 3), row.names = colnames(m)),
                     data.frame(feature_id = rownames(m)), omics_type = "rnaseq",
                     assay_type = "raw_count")
  b <- run_diff(inp, method = "edger", group_col = "group", control_group = "A", case_group = "B")
  expect_null(b$results$diff_object)
  withr::local_options(omicsCore.keep_count_models = TRUE)
  b2 <- run_diff(inp, method = "edger", group_col = "group", control_group = "A", case_group = "B")
  expect_false(is.null(b2$results$diff_object))
})

test_that("pairwise-complete correlation by products equals cor()", {
  set.seed(15)
  x <- matrix(stats::rnorm(500 * 8), 500)
  x[sample(length(x), 400)] <- NA
  expect_equal(unname(pairwise_cor(x)), unname(stats::cor(x, use = "pairwise.complete.obs")),
               tolerance = 1e-12)
})
