# =============================================================================
# Golden test: raw workbook -> omicsApp -> results that match the legacy pipeline
# =============================================================================
# This exists because of how the normalization gap got in. Phase 1's exit
# criterion was "all current CHISSS analyses can be reproduced via
# library(omicsCore)", and the validation script that signed it off did:
#
#     manager  <- readRDS(legacy_rds)      # already vsn-normalized and imputed
#     expr_mat <- manager$get_imputed()
#     input    <- omics_input(..., assay_type = "normalized_intensity")
#
# It started downstream of normalization, so it proved the differential layer
# was faithful and said nothing about the preprocessing layer -- which had not
# been ported. r = 1.000 against the legacy diff, and the gap sailed through.
#
# So this test starts at the raw .xlsx. Anything that only holds for
# pre-processed input cannot pass here.
#
# It needs the real SkinProteomics workbook, which does not live in this repo.
# Point CHISSS_PROTEOMICS_XLSX at it, or keep the conventional CHISSS layout.

golden_workbook <- function() {
  from_env <- Sys.getenv("CHISSS_PROTEOMICS_XLSX", "")
  if (nzchar(from_env) && file.exists(from_env)) return(from_env)

  # Walk up rather than counting "..": how deep the working directory sits
  # depends on whether the suite was started by test_file(), test_dir(),
  # R CMD check or the IDE, and getting that wrong silently skips the test.
  relative <- file.path("data", "SkinProteomics", "Proteomics_Data.xlsx")
  dir <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  repeat {
    candidate <- file.path(dir, relative)
    if (file.exists(candidate)) return(candidate)
    parent <- dirname(dir)
    if (identical(parent, dir)) break
    dir <- parent
  }
  NA_character_
}

skip_unless_golden <- function() {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("readxl")
  testthat::skip_if_not_installed("vsn")
  testthat::skip_if_not_installed("limma")
  path <- golden_workbook()
  testthat::skip_if(
    is.na(path),
    "Real SkinProteomics workbook not found; set CHISSS_PROTEOMICS_XLSX"
  )
  path
}

# What the legacy CHISSS framework does to get from raw intensities to
# normalized_intensity, transcribed from
# scripts/frameworks/omics_core/R/data_input/proteomics/build_proteomics_se.R.
#
# The log2 and the 2^ cancel on paper. They are kept here because they are
# kept there: 2^log2(x) != x at ~1e-15 for most real intensities, and vsn's
# optimiser turns that into visible differences. Reproducing the reference
# means reproducing the round trip.
legacy_normalize <- function(mat) {
  m <- mat
  m[m == 0] <- NA
  m <- suppressWarnings(log2(m))
  linear <- 2^m
  vsn::predict(vsn::vsnMatrix(linear), linear)
}

read_cheek_matrix <- function(path, max_missing = 0.5) {
  mat  <- suppressMessages(readxl::read_excel(path, sheet = 1))
  meta <- suppressMessages(readxl::read_excel(path, sheet = 2))

  m <- as.matrix(mat[, -1])
  storage.mode(m) <- "double"
  rownames(m) <- mat[[1]]

  keep <- intersect(meta$label[meta$tissue == "Cheek"], colnames(m))
  m <- m[, keep, drop = FALSE]
  m <- m[rowMeans(is.na(m)) <= max_missing, , drop = FALSE]

  meta_df <- as.data.frame(meta[match(colnames(m), meta$label), ])
  rownames(meta_df) <- meta_df$label
  list(
    expr = m,
    meta = meta_df,
    feat = data.frame(feature_id = rownames(m), row.names = rownames(m))
  )
}

# ---- C1: the normalized layer -----------------------------------------------

test_that("normalize_omics reproduces the legacy pipeline's normalized values", {
  path <- skip_unless_golden()
  d <- read_cheek_matrix(path)

  reference <- legacy_normalize(d$expr)
  produced <- suppressMessages(normalize_omics(
    omics_input(d$expr, d$meta, d$feat,
                omics_type = "proteomics", assay_type = "raw_intensity"),
    method = "vsn"
  ))$expr_mat

  expect_identical(dim(produced), dim(reference))
  expect_identical(dimnames(produced), dimnames(reference))

  # Not bit-identical, and the reason is known and bounded: the legacy path
  # sends vsn `2^log2(x)` rather than `x`, which differs at ~1e-15 relative
  # precision, and vsn's optimiser (factr = 5e7) amplifies that. Measured on
  # this dataset: r = 0.9999990, max |diff| 0.22, median 2.5e-03.
  r <- stats::cor(as.vector(reference), as.vector(produced),
                  use = "complete.obs")
  expect_gt(r, 0.9999)
  expect_lt(max(abs(reference - produced), na.rm = TRUE), 1.0)
  expect_lt(stats::median(abs(reference - produced), na.rm = TRUE), 0.01)

  # Both land on the same scale; a double or missing transform would not
  expect_lt(abs(max(produced, na.rm = TRUE) - max(reference, na.rm = TRUE)), 0.5)
  expect_silent(check_assay_scale(
    omics_input(produced, d$meta, d$feat, omics_type = "proteomics",
                assay_type = "normalized_intensity")
  ))
})

# ---- C1b: the conclusions ---------------------------------------------------

test_that("the differential results agree with the legacy pipeline", {
  path <- skip_unless_golden()
  d <- read_cheek_matrix(path)

  run_de <- function(mat_norm) {
    input <- omics_input(mat_norm, d$meta, d$feat, omics_type = "proteomics",
                         assay_type = "normalized_intensity")
    bundle <- run_diff(input, method = "limma", analysis_type = "group",
                       group_col = "condition", control_group = "G1",
                       case_group = "G2", covariates = "age")
    bundle$results$diff_result_df
  }

  reference <- run_de(legacy_normalize(d$expr))
  produced <- run_de(suppressMessages(normalize_omics(
    omics_input(d$expr, d$meta, d$feat,
                omics_type = "proteomics", assay_type = "raw_intensity"),
    method = "vsn"
  ))$expr_mat)

  joined <- merge(
    reference[, c("feature_id", "effect", "p_value", "adj_p_value")],
    produced[,  c("feature_id", "effect", "p_value", "adj_p_value")],
    by = "feature_id", suffixes = c(".ref", ".new")
  )
  expect_gt(nrow(joined), 3000)

  # Measured on this dataset: logFC r = 0.9999957, max |diff| 0.0088
  expect_gt(stats::cor(joined$effect.ref, joined$effect.new), 0.9999)
  expect_lt(max(abs(joined$effect.ref - joined$effect.new)), 0.05)
  expect_gt(
    stats::cor(-log10(joined$p_value.ref), -log10(joined$p_value.new)),
    0.9999
  )

  # What actually matters: the same proteins come out significant. Measured
  # 91 / 91 with Jaccard 1.0; the bound leaves room for an FDR-boundary flip.
  sig_ref <- joined$feature_id[joined$adj_p_value.ref < 0.05]
  sig_new <- joined$feature_id[joined$adj_p_value.new < 0.05]
  expect_gt(length(sig_ref), 50)
  jaccard <- length(intersect(sig_ref, sig_new)) / length(union(sig_ref, sig_new))
  expect_gt(jaccard, 0.95)
})

# ---- C2: the failure this was built to prevent -------------------------------

test_that("importing the raw workbook does not leave linear values in the project", {
  path <- skip_unless_golden()
  skip_if_not_installed("openxlsx")

  shiny::testServer(import_view_server, {
    session$setInputs(
      omics_type = "proteomics",
      file = list(datapath = path, name = basename(path),
                  size = file.info(path)$size)
    )
    expect_true(parse_ok())
    # The scale is read off the data, not assumed
    expect_identical(parsed()$input$assay_type, "raw_intensity")

    session$setInputs(normalize = TRUE, normalize_method = "vsn")
    session$setInputs(confirm = 1)

    inp <- session$returned()
    # This is the assertion the old validation could not make, because it
    # never started from a raw file: what reaches the project is on the scale
    # limma needs, and the untouched values are still there for QC.
    expect_identical(inp$assay_type, "normalized_intensity")
    expect_silent(check_assay_scale(inp))
    expect_gt(max(inp$raw_mat, na.rm = TRUE), 1e5)
    expect_lt(max(inp$expr_mat, na.rm = TRUE), 100)
  })
})

test_that("an already-normalized workbook is imported without a second transform", {
  path <- skip_unless_golden()
  skip_if_not_installed("openxlsx")

  # Same data, pre-normalized, which is what a collaborator would send
  d <- read_cheek_matrix(path)
  pre <- legacy_normalize(d$expr)

  xlsx <- tempfile(fileext = ".xlsx")
  on.exit(unlink(xlsx), add = TRUE)
  openxlsx::write.xlsx(
    list(
      expression = data.frame(feature_id = rownames(pre), pre,
                              check.names = FALSE),
      metadata = data.frame(sample_id = rownames(d$meta),
                            group = d$meta$condition, age = d$meta$age),
      feature_annot = data.frame(feature_id = rownames(pre))
    ),
    file = xlsx, overwrite = TRUE
  )

  shiny::testServer(import_view_server, {
    session$setInputs(
      omics_type = "proteomics",
      file = list(datapath = xlsx, name = "prenormalized.xlsx",
                  size = file.info(xlsx)$size)
    )
    expect_identical(parsed()$input$assay_type, "normalized_intensity")
    # No normalization offered, so no way to compress the range a second time
    expect_false(normalizable())

    session$setInputs(confirm = 1)
    inp <- session$returned()
    expect_identical(inp$assay_type, "normalized_intensity")
    expect_null(inp$raw_mat)
    # Values survive intact rather than being log-transformed again
    expect_equal(max(inp$expr_mat, na.rm = TRUE),
                 max(pre, na.rm = TRUE), tolerance = 1e-6)
  })
})

# ---- C3: projects saved before the vocabulary existed ------------------------

test_that("a layer labelled with the old 'intensity' value still analyses", {
  path <- skip_unless_golden()
  d <- read_cheek_matrix(path)

  # What the import view stamped on every proteomics upload before Phase A.
  # Such a layer is loaded from disk, so it never passes back through
  # omics_input() and keeps the superseded spelling.
  legacy_layer <- suppressWarnings(
    omics_input(legacy_normalize(d$expr), d$meta, d$feat,
                omics_type = "proteomics", assay_type = "normalized_intensity")
  )
  legacy_layer$assay_type <- "intensity"

  # Recognised as superseded rather than rejected or called unknown
  w <- tryCatch(validate_omics_input(legacy_layer),
                warning = function(w) conditionMessage(w))
  expect_match(w, "superseded")
  expect_false(grepl("Unrecognised", w))

  # And it still runs: an old project must not become unopenable
  bundle <- suppressWarnings(run_diff(
    legacy_layer, method = "limma", analysis_type = "group",
    group_col = "condition", control_group = "G1", case_group = "G2"
  ))
  expect_gt(nrow(bundle$results$diff_result_df), 3000)
})
