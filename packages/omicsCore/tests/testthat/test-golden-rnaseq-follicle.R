# Golden test: the follicle RNA-seq counts, raw file to result, against
# the legacy pipeline's stored output.
#
# The proteomics golden (in omicsApp) closed the gap that let a missing
# normalization step through: it starts at the raw workbook. The RNA-seq
# path had no such test, and the DESeq2/edgeR backends were the least
# covered code in the package. This one starts at all.counts.txt --
# 63,241 genes by 129 follicles, as the pipeline wrote it -- and walks
# the legacy G2-vs-G1 recipe with omicsCore: drop the novel-transcript
# genes, keep the grouped samples, filter with edgeR's filterByExpr,
# winsorize at k = 20, then edgeR QL with age as a covariate. The
# result table it compares against is the one the legacy report saved.
#
# Needs the parent project's data/ and results/, which are not in this
# repository. Set OMICSAPP_FOLLICLE_ROOT to the project root, or keep
# the conventional layout; the test skips when either is absent.

follicle_paths <- function() {
  root <- Sys.getenv("OMICSAPP_FOLLICLE_ROOT", "")
  if (!nzchar(root)) {
    dir <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
    repeat {
      if (file.exists(file.path(dir, "data", "Transcriptomic-Follicle", "all.counts.txt"))) {
        root <- dir
        break
      }
      parent <- dirname(dir)
      if (identical(parent, dir)) break
      dir <- parent
    }
  }
  if (!nzchar(root)) return(NULL)
  results <- file.path(root, "results", "rnaseq", "follicle")
  p <- list(
    counts   = file.path(root, "data", "Transcriptomic-Follicle", "all.counts.txt"),
    novel    = file.path(results, "Tables", "removed_novel_description_genes.tsv"),
    manifest = file.path(results, "G2_vs_G1", "Tables", "sample_manifest_G1_G2.tsv"),
    filter   = file.path(results, "G2_vs_G1", "Tables", "feature_filter_summary.tsv"),
    winsor   = file.path(results, "G2_vs_G1", "Tables", "qc_winsorize_stats.tsv"),
    deg      = file.path(results, "G2_vs_G1", "Tables", "diff.summary.xlsx")
  )
  if (!all(vapply(p, file.exists, logical(1)))) return(NULL)
  p
}

skip_unless_follicle <- function() {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("edgeR")
  testthat::skip_if_not_installed("readxl")
  p <- follicle_paths()
  testthat::skip_if(is.null(p),
    "Follicle counts and legacy results not found; set OMICSAPP_FOLLICLE_ROOT")
  p
}

# Read once per session: the file is 22 MB and edgeR takes its time.
follicle_cache <- new.env(parent = emptyenv())

follicle_import <- function(p) {
  if (!is.null(follicle_cache$import)) return(follicle_cache$import)
  out <- suppressWarnings(read_omics(p$counts, omics_type = "rnaseq",
                                     assay_type = "raw_count"))
  follicle_cache$import <- out
  out
}

follicle_prepared <- function(p) {
  if (!is.null(follicle_cache$prepared)) return(follicle_cache$prepared)
  counts <- follicle_import(p)$input$expr_mat
  novel <- utils::read.delim(p$novel, stringsAsFactors = FALSE)$gene_id
  manifest <- utils::read.delim(p$manifest, stringsAsFactors = FALSE)
  manifest <- manifest[manifest$in_analysis %in% c(TRUE, "TRUE"), , drop = FALSE]

  counts <- counts[!rownames(counts) %in% novel, manifest$sample_id, drop = FALSE]
  meta <- data.frame(Group = manifest$Group, AGE = manifest$AGE,
                     row.names = manifest$sample_id, stringsAsFactors = FALSE)

  keep_nonzero <- rowSums(counts) > 0
  design <- stats::model.matrix(~ Group, data = meta)
  keep_expr <- edgeR::filterByExpr(edgeR::DGEList(counts = counts), design = design)
  counts <- counts[keep_nonzero & keep_expr, , drop = FALSE]

  wins <- winsorize_counts(counts, k = 20)
  clipped <- round(wins$count_mat)
  storage.mode(clipped) <- "integer"

  input <- omics_input(clipped, meta, data.frame(feature_id = rownames(clipped)),
                       omics_type = "rnaseq", assay_type = "raw_count")
  follicle_cache$prepared <- list(input = input, winsor = wins$stats,
                                  n_nonzero = sum(keep_nonzero),
                                  n_kept = nrow(clipped))
  follicle_cache$prepared
}

# ---- the reader, on the real file -------------------------------------------

test_that("the pipeline's counts file imports whole, ids and samples intact", {
  p <- skip_unless_follicle()
  out <- follicle_import(p)
  m <- out$input$expr_mat
  expect_identical(dim(m), c(63241L, 129L))
  expect_true(all(startsWith(rownames(m), "ENSG")))
  expect_true(all(grepl("^RD[0-9]{3}_Folli$", colnames(m))))
  expect_true(all(m == round(m)))
  expect_false(anyNA(m))
  # The novel-transcript genes the legacy run removed are all present
  novel <- utils::read.delim(p$novel, stringsAsFactors = FALSE)$gene_id
  expect_true(all(novel %in% rownames(m)))
})

# ---- the preprocessing, step by step against the legacy tables ------------------

test_that("the feature filter keeps exactly the genes the legacy run kept", {
  p <- skip_unless_follicle()
  prep <- follicle_prepared(p)
  summary <- utils::read.delim(p$filter, stringsAsFactors = FALSE)
  expect_identical(prep$n_nonzero,
                   summary$n_features[summary$step == "nonzero_total_count"])
  expect_identical(prep$n_kept,
                   summary$n_features[summary$step == "final_kept_features"])
})

test_that("winsorize_counts reproduces the legacy clipping, gene by gene", {
  p <- skip_unless_follicle()
  prep <- follicle_prepared(p)
  legacy <- utils::read.delim(p$winsor, stringsAsFactors = FALSE)
  ours <- prep$winsor[match(legacy$feature_id, prep$winsor$feature_id), ]
  expect_false(anyNA(ours$feature_id))
  expect_equal(ours$threshold, legacy$threshold)
  expect_identical(as.integer(ours$n_clipped), as.integer(legacy$n_clipped))
  expect_equal(ours$max_original, legacy$max_original)
  # And nothing outside the legacy table was clipped
  clipped_elsewhere <- prep$winsor$n_clipped[!prep$winsor$feature_id %in% legacy$feature_id]
  expect_true(all(clipped_elsewhere == 0L))
})

# ---- the conclusions -------------------------------------------------------------

test_that("edgeR with age as covariate reproduces the legacy table", {
  p <- skip_unless_follicle()
  prep <- follicle_prepared(p)
  bundle <- run_diff(prep$input, method = "edger", analysis_type = "group",
                     group_col = "Group", control_group = "G1", case_group = "G2",
                     covariates = "AGE")
  got <- bundle$results$diff_result_df
  legacy <- as.data.frame(readxl::read_excel(p$deg, sheet = "Raw"))
  expect_setequal(got$feature_id, legacy$feature_id)
  got <- got[match(legacy$feature_id, got$feature_id), ]

  # Two implementations of the same recipe: an equality, not a tolerance
  # in any meaningful sense.
  expect_equal(got$effect, legacy$logFC, tolerance = 1e-8)
  expect_equal(got$statistic, legacy$F, tolerance = 1e-8)
  expect_equal(got$p_value, legacy$PValue, tolerance = 1e-8)
  expect_equal(got$adj_p_value, legacy$FDR, tolerance = 1e-8)
  expect_equal(got$base_mean, legacy$logCPM, tolerance = 1e-8)

  # What actually matters: the same genes come out, at the cut the
  # legacy report used (raw p < 0.05, |log2FC| > 0.2)
  sig_ref <- legacy$feature_id[legacy$PValue < 0.05 & abs(legacy$logFC) > 0.2]
  sig_new <- got$feature_id[got$p_value < 0.05 & abs(got$effect) > 0.2]
  expect_gt(length(sig_ref), 1000L)
  expect_setequal(sig_ref, sig_new)
})

test_that("the count backends survive the real matrix, constant genes and all", {
  p <- skip_unless_follicle()
  prep <- follicle_prepared(p)
  # PCA on the full 63k-gene import, where a fifth of the genes never
  # vary -- the case that once stopped the QC scatter.
  full <- follicle_import(p)$input$expr_mat
  pca <- pca_over_samples(log2(full + 1))
  expect_gt(attr(pca, "n_dropped"), 5000L)
  expect_identical(nrow(pca$x), 129L)
  expect_true(all(is.finite(pca$sdev)))
})
