# The volcano's cloud, thinned.
#
# A 63k-feature volcano is 6.5 MB of JSON once plotly has serialised it,
# and the browser then draws 63,000 marks. The R side takes half a second,
# which is why this never appeared in a server-side timing and reached the
# user as "the differential analysis is very slow".
#
# What must survive thinning is everything that carries information. What
# may go is the grey cloud, which at that density is solid.

thin_df <- function(n_ns = 5000L, n_sig = 40L, n_lab = 5L) {
  df <- data.frame(
    feature_id = paste0("f", seq_len(n_ns + n_sig)),
    effect = c(stats::rnorm(n_ns), stats::rnorm(n_sig, 3)),
    stringsAsFactors = FALSE)
  df$.sig <- factor(c(rep("ns", n_ns), rep("significant", n_sig)),
                    levels = c("ns", "significant"))
  df$.label <- NA_character_
  df$.label[seq_len(n_lab)] <- paste0("L", seq_len(n_lab))  # labelled but ns
  df
}

test_that("nothing significant or labelled is ever dropped", {
  df <- thin_df()
  out <- thin_volcano_points(df, max_points = 500L)

  expect_lte(nrow(out), nrow(df))
  expect_equal(sum(out$.sig == "significant"), sum(df$.sig == "significant"))
  expect_equal(sum(!is.na(out$.label)), sum(!is.na(df$.label)))
  # a labelled non-significant point is the case a naive filter loses
  expect_true(all(df$feature_id[!is.na(df$.label)] %in% out$feature_id))
})

test_that("the cap is respected, and only the cloud pays for it", {
  df <- thin_df(n_ns = 5000L, n_sig = 40L, n_lab = 5L)
  out <- thin_volcano_points(df, max_points = 500L)
  expect_lte(nrow(out), 500L)
  expect_gt(sum(out$.sig == "ns"), 0L)
})

test_that("a cap below the number of protected points keeps them all", {
  # Better an over-full plot than one missing its findings.
  df <- thin_df(n_ns = 100L, n_sig = 300L, n_lab = 0L)
  out <- thin_volcano_points(df, max_points = 50L)
  expect_equal(sum(out$.sig == "significant"), 300L)
})

test_that("below the cap nothing changes at all", {
  df <- thin_df(n_ns = 100L, n_sig = 10L)
  expect_identical(thin_volcano_points(df, max_points = 10000L), df)
  expect_identical(thin_volcano_points(df, max_points = NULL), df)
  expect_identical(thin_volcano_points(df, max_points = 0L), df)
})

# A figure that reshuffles between two renders of one result is a figure
# nobody can cite, and one that consumes random numbers changes whatever
# the caller does next.
test_that("thinning is deterministic and does not touch the RNG", {
  df <- thin_df()
  expect_identical(thin_volcano_points(df, 500L),
                   thin_volcano_points(df, 500L))

  set.seed(1); before <- .Random.seed
  invisible(thin_volcano_points(df, 500L))
  expect_identical(.Random.seed, before)
})

make_thin_bundle <- function(n = 1200L) {
  set.seed(3)
  ids <- paste0("f", seq_len(n))
  m <- matrix(stats::rnorm(n * 8, 10), nrow = n,
              dimnames = list(ids, paste0("s", 1:8)))
  m[1:30, 5:8] <- m[1:30, 5:8] + 4
  x <- omics_input(m,
                   data.frame(group = rep(c("a", "b"), each = 4),
                              row.names = colnames(m)),
                   data.frame(feature_id = ids, feature_symbol = ids,
                              row.names = ids),
                   omics_type = "proteomics",
                   assay_type = "normalized_intensity")
  run_diff(x, method = "ttest", analysis_type = "group", group_col = "group",
           control_group = "a", case_group = "b")
}

test_that("plot_volcano says so when it has thinned", {
  skip_if_not_installed("ggplot2")
  b <- make_thin_bundle()
  cap_all <- ggplot2::ggplot_build(plot_volcano(b, max_points = NULL))$plot$labels$caption
  cap_cut <- ggplot2::ggplot_build(plot_volcano(b, max_points = 200L))$plot$labels$caption
  expect_false(grepl("thinned", cap_all %||% ""))
  expect_match(cap_cut, "thinned")
  # and it says how many of how many, since the subtitle counts them all
  expect_match(cap_cut, "of [0-9,]+ features")
})
