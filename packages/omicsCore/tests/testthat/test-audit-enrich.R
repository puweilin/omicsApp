# Enrichment, integration and export defects found in the 2026-10 audit.

ae_diff_df <- function(effect, statistic, p, symbol, type = "t") {
  data.frame(feature_id = paste0("f", seq_along(effect)), feature_symbol = symbol,
             effect = effect, statistic = statistic, statistic_type = type,
             p_value = p, stringsAsFactors = FALSE)
}

test_that("GSEA ranks by the signed statistic, one value per gene, tested genes only", {
  df <- ae_diff_df(effect = c(3, 0.5, -2, 1, 4),
                   statistic = c(1, 9, -8, 2, NA),
                   p = c(0.3, 1e-6, 1e-5, 0.05, NA),
                   symbol = c("A", "B", "C", "C", "D"))
  r <- gsea_rank_vector(df, "feature_symbol")
  # B has the small effect and the strong evidence; D was never tested.
  expect_identical(names(r), c("B", "A", "C"))
  # C measured twice: the stronger value (-8) wins, not the more positive.
  expect_equal(unname(r[["C"]]), -8)
  expect_identical(attr(r, "metric"), "signed test statistic")

  # edgeR's QL F is a squared t.
  f <- ae_diff_df(effect = c(1, -1), statistic = c(16, 4), p = c(0.01, 0.1),
                  symbol = c("X", "Y"), type = "F")
  expect_equal(as.numeric(gsea_rank_vector(f, "feature_symbol")), c(4, -2))
})

test_that("ORA passes the gene-set size limits on, and its universe is what was tested", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("msigdbr")
  b <- realistic_diff_bundle()
  enr <- run_enrichment(b, type = "ora", database = "hallmark",
                        output_p_cutoff = 1, min_size = 30L, max_size = 60L)
  obj <- enr$results$enrich_object[[1L]]
  skip_if(is.null(obj), "no ORA result to inspect")
  bg <- as.integer(sub("/.*", "", obj@result$BgRatio))
  expect_true(all(bg >= 30L & bg <= 60L))
  expect_identical(enr$params$output_p_cutoff, 1)

  # A gene with no p-value is not in the background.
  b2 <- b
  df <- b2$results$diff_result_df
  untested <- utils::tail(order(df$p_value), 50L)   # the least significant
  df$p_value[untested] <- NA
  df$adj_p_value[untested] <- NA
  b2$results$diff_result_df <- df
  enr2 <- run_enrichment(b2, type = "ora", database = "hallmark",
                         output_p_cutoff = 1)
  skip_if(length(enr2$results$enrich_object) == 0L, "no ORA result to inspect")
  obj2 <- enr2$results$enrich_object[[1L]]
  n_bg <- as.integer(sub(".*/", "", obj2@result$BgRatio[1L]))
  n_bg_full <- as.integer(sub(".*/", "", obj@result$BgRatio[1L]))
  expect_lt(n_bg, n_bg_full)
})

test_that("concordance quadrant counts are over the features significant in both", {
  set.seed(11)
  sym <- c(paste0("S", 1:10), paste0("G", 1:30))
  groups <- rep(c("ctrl", "case"), each = 4)
  mk <- function(prefix, omics, assay) {
    ids <- paste0(prefix, seq_along(sym))
    samp <- paste0(prefix, "_S", 1:8)
    m <- matrix(stats::rnorm(length(sym) * 8, 8, 0.3), length(sym),
                dimnames = list(ids, samp))
    m[1:10, groups == "case"] <- m[1:10, groups == "case"] + 3
    omics_input(m, data.frame(group = groups, row.names = samp),
                data.frame(feature_id = ids, feature_symbol = sym, row.names = ids),
                omics_type = omics, assay_type = assay)
  }
  p <- omics_project("q", list(prot = mk("p", "proteomics", "normalized_intensity"),
                               rna = mk("r", "rnaseq", "logcpm")))
  d <- lapply(p$experiments, run_diff, method = "ttest", group_col = "group",
              control_group = "ctrl", case_group = "case")
  b <- run_integration(p, "concordance", c("prot", "rna"), diff_bundles = d)
  info <- b$params$method_info
  expect_gt(info$n_significant_both, 0L)
  expect_equal(sum(unlist(info$quadrant_counts)), info$n_significant_both)
  expect_gte(sum(unlist(info$quadrant_counts_all)), info$n_significant_both)
  g <- plot_integration(b, view = "quadrant")
  expect_equal(sum(g$data$n), sum(b$results$integration_df$significant_a &
                                    b$results$integration_df$significant_b &
                                    !is.na(b$results$integration_df$quadrant)))
})

test_that("tables with list columns or tabs export as plain cells", {
  df <- data.frame(id = c("a", "b"), stringsAsFactors = FALSE)
  df$genes <- list(c("TP53", "MYC"), character(0))
  df$note <- c("tab\there", "line\nbreak")
  out <- flatten_for_export(df)
  expect_identical(out$genes, c("TP53;MYC", NA))
  expect_identical(out$note, c("tab here", "line break"))
  path <- withr::local_tempfile(fileext = ".tsv")
  utils::write.table(out, path, sep = "\t", quote = FALSE, row.names = FALSE)
  expect_identical(nrow(utils::read.delim(path)), 2L)
})

test_that("the script reproduces an in-app normalization and the study design", {
  inp <- realistic_input()
  inp$expr_mat <- 2^inp$expr_mat
  inp$assay_type <- "raw_intensity"
  norm <- suppressMessages(normalize_omics(inp, method = "log2"))
  norm <- set_study_design(norm, "group", "G1")
  norm$source_path <- "upload.xlsx"
  proj <- omics_project("n", list(proteomics = norm))
  lines <- export_script(proj)
  read_line <- grep("read_omics(", lines, fixed = TRUE, value = TRUE)
  expect_true(any(grepl("raw_intensity", lines, fixed = TRUE)))
  expect_false(any(grepl("normalized_intensity", read_line, fixed = TRUE)))
  expect_true(any(grepl("normalize_omics(", lines, fixed = TRUE)))
  expect_true(any(grepl("set_study_design(", lines, fixed = TRUE)))
  # The emitted lines parse.
  expect_silent(parse(text = lines))
})

test_that("subsetting keeps where the data came from and what was done to it", {
  inp <- realistic_input()
  inp$source_path <- "upload.xlsx"
  inp$normalization <- list(method = "log2", offset = 1, from_assay_type = "raw_intensity")
  sub <- subset_omics(inp, samples = colnames(inp$expr_mat)[1:6])
  expect_identical(sub$source_path, "upload.xlsx")
  expect_identical(sub$normalization$method, "log2")
})

ae_int_project <- function() {
  set.seed(21)
  sym <- c(paste0("S", 1:10), paste0("G", 1:30))
  groups <- rep(c("ctrl", "trtA", "trtB"), each = 4)
  mk <- function(prefix, omics, assay, shift) {
    ids <- paste0(prefix, seq_along(sym))
    samp <- paste0(prefix, "_S", 1:12)
    m <- matrix(stats::rnorm(length(sym) * 12, 8, 0.3), length(sym),
                dimnames = list(ids, samp))
    m[1:10, groups == "trtA"] <- m[1:10, groups == "trtA"] + shift
    omics_input(m, data.frame(group = groups, donor = paste0("D", 1:12), row.names = samp),
                data.frame(feature_id = ids, feature_symbol = sym, row.names = ids),
                omics_type = omics, assay_type = assay)
  }
  omics_project("int", list(proteomics = mk("p", "proteomics", "normalized_intensity", 2),
                            rnaseq = mk("r", "rnaseq", "logcpm", 1.5)))
}

# Runs an exported script with read_omics() answering from memory.
ae_run_script <- function(lines, inputs) {
  env <- new.env(parent = asNamespace("omicsCore"))
  env$read_omics <- function(path, ...) {
    list(input = inputs[[sub("^<path-to-(.*)-file>$", "\\1", path)]])
  }
  env$library <- function(...) invisible(NULL)
  env$sessionInfo <- function() invisible(NULL)
  pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  eval(parse(text = lines), envir = env)
  env
}

test_that("the script reproduces a concordance integration, partner diff and settings included", {
  p <- ae_int_project()
  full <- lapply(p$experiments, run_diff, method = "limma", group_col = "group",
                 control_group = "ctrl", case_group = c("trtA", "trtB"))
  d <- lapply(full, select_comparison, "trtA_vs_ctrl")
  p$bundles <- list(
    diff = full$proteomics,
    integration = run_integration(p, "concordance", c("proteomics", "rnaseq"),
                                  diff_bundles = d, p_cutoff = 0.01,
                                  effect_cutoff = 0.5))
  lines <- export_script(p)
  expect_true(any(grepl("p_cutoff", lines)))
  expect_true(any(grepl("effect_cutoff", lines)))
  expect_true(any(grepl("diff_rnaseq <- run_diff(", lines, fixed = TRUE)))
  expect_true(any(grepl('view = "quadrant"', lines, fixed = TRUE)))
  expect_false(any(grepl("NOTE: run_integration", lines, fixed = TRUE)))
  env <- ae_run_script(lines, p$experiments)
  expect_equal(env$integration$results$integration_df,
               p$bundles$integration$results$integration_df)
})

test_that("the script reproduces a correlation integration with its sample pairing", {
  p <- ae_int_project()
  link <- data.frame(
    tag = rep(c("proteomics", "rnaseq"), each = 12),
    sample_id = c(paste0("p_S", 1:12), paste0("r_S", 12:1)),
    donor_id = c(paste0("D", 1:12), paste0("D", 1:12)),
    stringsAsFactors = FALSE)
  p$sample_link <- link
  p$bundles <- list(integration = run_integration(
    p, "correlation", c("proteomics", "rnaseq"), cor_method = "pearson", min_samples = 5L))
  lines <- export_script(p)
  expect_true(any(grepl("sample_link = data.frame(", lines, fixed = TRUE)))
  expect_true(any(grepl('cor_method', lines)))
  expect_true(any(grepl('view = "scatter"', lines, fixed = TRUE)))
  expect_false(any(grepl("dual_volcano", lines)))
  env <- ae_run_script(lines, p$experiments)
  expect_equal(env$integration$results$integration_df,
               p$bundles$integration$results$integration_df)
})
