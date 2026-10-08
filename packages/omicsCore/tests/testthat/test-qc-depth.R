# The missingness panels answer the proteomics question: a peptide that
# was not detected is a hole in the matrix. A counts matrix has no holes
# -- every gene has a number for every sample and most are zero -- so
# the panel reported "63,241 features, all at 0%", which is true, carries
# no information, and occupied the space that should have been showing
# whether a library was under-sequenced.

depth_input <- function(lib_scale = rep(1, 6), n_feat = 40L) {
  ids <- paste0("S", seq_along(lib_scale))
  set.seed(7)
  base <- matrix(as.numeric(stats::rpois(n_feat * length(ids), 100)),
                 nrow = n_feat, dimnames = list(paste0("G", seq_len(n_feat)), ids))
  mat <- sweep(base, 2L, lib_scale, "*")
  omics_input(mat,
              data.frame(sample_id = ids, condition = "G1", row.names = ids,
                         stringsAsFactors = FALSE),
              data.frame(feature_id = rownames(mat), row.names = rownames(mat),
                         stringsAsFactors = FALSE),
              omics_type = "rnaseq", assay_type = "raw_count")
}

test_that("qc_depth reports library size and detection per sample", {
  d <- qc_depth(depth_input())
  expect_equal(nrow(d), 6L)
  expect_setequal(names(d), c("sample_id", "library_size", "n_detected",
                              "detection_rate", "library_size_ratio"))
  expect_true(all(d$library_size > 0))
  expect_true(all(d$detection_rate <= 1))
})

test_that("a shallow library is flagged against the median, not a constant", {
  # What counts as shallow depends entirely on the experiment; a fixed
  # count would be wrong for every study but one.
  d <- qc_depth(depth_input(c(1, 1, 1, 1, 1, 0.1)))
  expect_identical(qc_depth_outliers(d), "S6")

  # The same matrix scaled up is not suddenly healthy.
  d2 <- qc_depth(depth_input(c(1, 1, 1, 1, 1, 0.1) * 1000))
  expect_identical(qc_depth_outliers(d2), "S6")
})

test_that("an even cohort flags nothing", {
  expect_length(qc_depth_outliers(qc_depth(depth_input())), 0L)
})

test_that("detection counts signal, whether absence is written 0 or NA", {
  # A counts matrix says "nothing seen" with 0 and an intensity matrix
  # with NA. The question is the same one.
  inp <- depth_input()
  inp$expr_mat[1:10, "S1"] <- 0
  inp$expr_mat[1:10, "S2"] <- NA_real_
  d <- qc_depth(inp)
  n <- nrow(inp$expr_mat)
  expect_equal(d$n_detected[d$sample_id == "S1"], n - 10L)
  expect_equal(d$n_detected[d$sample_id == "S2"], n - 10L)
})

test_that("run_qc carries a depth summary for both modalities", {
  for (type in c("rnaseq", "proteomics")) {
    inp <- depth_input()
    inp$omics_type <- type
    inp$assay_type <- if (type == "rnaseq") "raw_count" else "raw_intensity"
    b <- run_qc(inp)
    expect_equal(nrow(b$results$qc_summary$depth), 6L, info = type)
  }
})

test_that("the depth view draws, and refuses a bundle that has none", {
  b <- run_qc(depth_input())
  expect_s3_class(plot_qc(b, view = "depth"), c("patchwork", "ggplot"))

  b$results$qc_summary$depth <- NULL
  expect_error(plot_qc(b, view = "depth"), "no depth summary")
})

# ---- the depth view ----------------------------------------------------
# Same rules as the missingness panel's samples: it named up to 30
# samples in the app's 360 px card, and past ~10 the names overlapped.

depth_df <- function(n, low = integer(0), ids = sprintf("S%02d", seq_len(n))) {
  lib <- seq(1.3, 0.7, length.out = n) * 1e6
  lib[low] <- 0.1e6
  det <- round(seq(15000, 14000, length.out = n))
  det[low] <- 9000
  data.frame(sample_id = ids, library_size = lib, n_detected = as.integer(det),
             detection_rate = det / 20000,
             library_size_ratio = lib / stats::median(lib),
             stringsAsFactors = FALSE)
}
depth_bundle <- function(d) {
  b <- run_qc(depth_input())
  b$results$qc_summary$depth <- d
  b
}
depth_panels <- function(d) {
  p <- plot_qc(depth_bundle(d), view = "depth")
  list(library = p[[1]], detection = p[[2]])
}

test_that("depth names samples on bars up to the shared cap, then ranks them", {
  bars <- depth_panels(depth_df(SAMPLE_MAX_NAMED_BARS))
  expect_false("rank" %in% names(bars$library$data))
  expect_false("rank" %in% names(bars$detection$data))
  curve <- depth_panels(depth_df(SAMPLE_MAX_NAMED_BARS + 1L))
  expect_true("rank" %in% names(curve$library$data))
  expect_true("rank" %in% names(curve$detection$data))
  # Every sample stays on the curve.
  expect_equal(nrow(depth_panels(depth_df(60L))$library$data), 60L)
})

test_that("depth bar names shrink with the missingness panel's, worst at the top", {
  p <- depth_panels(depth_df(SAMPLE_MAX_NAMED_BARS, low = 3L))$library
  expect_equal(p$theme$axis.text.y$size, 7)
  # A discrete y axis draws bottom-up: the last level is the top row.
  lv <- levels(p$data$sample_id)
  expect_identical(lv[length(lv)], "S03")
  expect_null(depth_panels(depth_df(4L))$library$theme$axis.text.y$size)
})

test_that("shallow libraries are amber in both panels, with the cutoff dashed", {
  d <- depth_df(8L, low = c(2L, 5L))
  panels <- depth_panels(d)
  lib <- panels$library
  expect_setequal(as.character(lib$data$sample_id[lib$data$.low]), c("S02", "S05"))
  expect_setequal(as.character(lib$data$sample_id[lib$data$.low]),
                  qc_depth_outliers(d))
  # The dashed line, over a white halo that keeps it visible on the bars.
  vl <- Filter(function(l) inherits(l$geom, "GeomVline"), lib$layers)
  dashed <- Filter(function(l) identical(l$aes_params$linetype, "dashed"), vl)
  expect_length(dashed, 1L)
  expect_equal(dashed[[1]]$data$xintercept %||% dashed[[1]]$aes_params$xintercept,
               DEPTH_LOW_RATIO * stats::median(d$library_size))
  expect_identical(vl[[1]]$aes_params$colour, "white")
  expect_match(lib$labels$subtitle, "2 below 30% of median", fixed = TRUE)
  fills <- unique(ggplot2::ggplot_build(lib)$data[[1]]$fill)
  expect_setequal(fills, c(MISSING_FILL, MISSING_OVER_FILL))
  # Amber, not the "up" red it used to be.
  expect_false(omics_colors$up %in% fills)
  det <- panels$detection
  expect_setequal(as.character(det$data$sample_id[det$data$.low]), c("S02", "S05"))
  expect_identical(det$labels$subtitle, "of 20,000 \u00b7 amber: shallow library")
})

test_that("the ranked depth curve names the shallow or lowest samples", {
  ids <- sprintf("Patient_%03d_PBMC_RNA_rep1", 1:24)
  p <- depth_panels(depth_df(24L, low = c(4L, 9L), ids = ids))
  expect_match(p$library$labels$subtitle,
               "24 samples \u00b7 2 below 30% of median\nShallow: ",
               fixed = TRUE)
  expect_match(p$library$labels$subtitle, ids[4], fixed = TRUE)
  expect_match(p$detection$labels$subtitle, "\nFewest: ", fixed = TRUE)
  # Rank 1 is the lowest.
  lib <- p$library$data
  expect_equal(lib$library_size[lib$rank == 1L], min(lib$library_size))
  even <- depth_panels(depth_df(24L))$library
  expect_match(even$labels$subtitle, "none below 30% of median", fixed = TRUE)
  expect_match(even$labels$subtitle, "\nLowest: S24", fixed = TRUE)
})

test_that("each depth panel starts its own axis at zero", {
  p <- depth_panels(depth_df(8L))
  r1 <- ggplot2::ggplot_build(p$library)$layout$panel_params[[1]]$x.range
  r2 <- ggplot2::ggplot_build(p$detection)$layout$panel_params[[1]]$x.range
  expect_equal(r1[1], 0)
  expect_equal(r2[1], 0)
  expect_gt(r1[2], 1e6)
  expect_lt(r2[2], 2e4)
})

test_that("depth axis labels format each value on its own", {
  # format() over the whole vector wrote the zero as "0e+00".
  expect_identical(depth_axis_labels(c(0, 500, 50000, 2.5e6, NA)),
                   c("0", "500", "50k", "2.5M", ""))
})

test_that("depth subtitles fit a phone-width line", {
  # About 42 characters fit the app's 293 px panel; "of the median
  # (dashed)" ran off it.
  ids <- sprintf("Patient_%03d_PBMC_RNA_rep1", 1:60)
  p <- depth_panels(depth_df(60L, low = c(2L, 9L, 30L), ids = ids))
  lines <- unlist(strsplit(c(p$library$labels$subtitle,
                             p$detection$labels$subtitle), "\n"))
  expect_true(all(nchar(lines) <= 42L), info = paste(lines, collapse = " | "))
})

test_that("the depth view survives the app's phone-width theme", {
  for (n in c(8L, 24L, 60L)) {
    ids <- sprintf("Patient_%03d_PBMC_RNA_rep1", seq_len(n))
    p <- plot_qc(depth_bundle(depth_df(n, low = 2L, ids = ids)), view = "depth") &
      ggplot2::theme(text = ggplot2::element_text(size = 9))
    f <- tempfile(fileext = ".png")
    grDevices::png(f, width = 293, height = 360, res = 96)
    expect_no_error(print(p))
    grDevices::dev.off()
    unlink(f)
  }
})
