# Relationships between runs.
#
# The most dangerous bug class in this project's history is a complete,
# plausible, wrong table: is_significant answered FALSE for everything,
# GSEA obeyed a feature threshold it must ignore, enrichment ran on a
# different gene list from the one on screen, PCA died on a constant
# gene. Edge tests do not catch that class, because the wrong table has
# the right shape. Relationships do: swap the contrast and the effect
# must flip, shuffle the samples and nothing may change, tighten a
# threshold and the hit list may only shrink.

diff_on <- function(input, method, control = "G1", case = "G2", ...) {
  suppressMessages(run_diff(
    input, method = method, analysis_type = "group",
    group_col = "group", control_group = control, case_group = case, ...
  ))$results$diff_result_df
}

aligned <- function(df, ids) df[match(ids, df$feature_id), , drop = FALSE]

continuous_engines <- function() {
  c("ttest", "lm", if (requireNamespace("limma", quietly = TRUE)) "limma")
}
count_engines <- function() {
  c(if (requireNamespace("DESeq2", quietly = TRUE)) "deseq2",
    if (requireNamespace("edgeR", quietly = TRUE)) "edger")
}
tol_for <- function(method) if (method %in% c("deseq2", "edger")) 1e-8 else 1e-12

input_for <- function(method, ...) {
  if (method %in% c("deseq2", "edger")) realistic_input("rnaseq", n_per_group = 4L, ...)
  else realistic_input(n_per_group = 5L, ...)
}

permute_samples <- function(input, seed = 3) {
  set.seed(seed)
  idx <- sample(ncol(input$expr_mat))
  out <- input
  out$expr_mat <- input$expr_mat[, idx, drop = FALSE]
  out$meta_df <- input$meta_df[idx, , drop = FALSE]
  out
}

rename_samples <- function(input) {
  out <- input
  new <- paste0("X_", colnames(input$expr_mat))
  colnames(out$expr_mat) <- new
  rownames(out$meta_df) <- new
  out
}

# ---- every engine ---------------------------------------------------------

for (method in c(continuous_engines(), count_engines())) {
  tol <- tol_for(method)

  test_that(paste(method, ": swapping control and case negates the effect and keeps p"), {
    inp <- input_for(method)
    a <- diff_on(inp, method)
    z <- aligned(diff_on(inp, method, control = "G2", case = "G1"), a$feature_id)
    # DESeq2 re-parameterises the GLM when the reference level changes
    # and its IRLS fit converges to a slightly different point; the two
    # answers agree to about five digits, not to machine precision.
    swap_tol <- if (method == "deseq2") 1e-4 else tol
    expect_equal(z$effect, -a$effect, tolerance = swap_tol)
    expect_equal(z$p_value, a$p_value, tolerance = swap_tol)
    expect_equal(z$adj_p_value, a$adj_p_value, tolerance = swap_tol)
    flipped <- c(up = "down", down = "up", ns = "ns")
    expect_identical(z$direction, unname(flipped[a$direction]))
  })

  test_that(paste(method, ": the sample order does not matter"), {
    inp <- input_for(method)
    a <- diff_on(inp, method)
    z <- aligned(diff_on(permute_samples(inp), method), a$feature_id)
    expect_equal(z$effect, a$effect, tolerance = tol)
    expect_equal(z$p_value, a$p_value, tolerance = tol)
  })

  test_that(paste(method, ": the feature order does not matter"), {
    inp <- input_for(method)
    set.seed(4)
    idx <- sample(nrow(inp$expr_mat))
    shuffled <- inp
    shuffled$expr_mat <- inp$expr_mat[idx, , drop = FALSE]
    shuffled$feature_df <- inp$feature_df[idx, , drop = FALSE]
    a <- diff_on(inp, method)
    z <- aligned(diff_on(shuffled, method), a$feature_id)
    expect_equal(z$effect, a$effect, tolerance = tol)
    expect_equal(z$p_value, a$p_value, tolerance = tol)
    expect_equal(z$adj_p_value, a$adj_p_value, tolerance = tol)
  })

  test_that(paste(method, ": sample ids are labels, not data"), {
    inp <- input_for(method)
    a <- diff_on(inp, method)
    z <- aligned(diff_on(rename_samples(inp), method), a$feature_id)
    expect_equal(z$effect, a$effect, tolerance = tol)
    expect_equal(z$p_value, a$p_value, tolerance = tol)
  })

  test_that(paste(method, ": a duplicated feature gets its twin's statistics"), {
    inp <- input_for(method)
    twin <- inp
    extra <- inp$expr_mat[1, , drop = FALSE]
    rownames(extra) <- "TWIN"
    twin$expr_mat <- rbind(inp$expr_mat, extra)
    twin$feature_df <- rbind(inp$feature_df,
                             data.frame(feature_id = "TWIN", feature_symbol = "TWIN"))
    df <- diff_on(twin, method)
    one <- df[df$feature_id == rownames(inp$expr_mat)[1], ]
    two <- df[df$feature_id == "TWIN", ]
    expect_equal(two$effect, one$effect, tolerance = tol)
    expect_equal(two$p_value, one$p_value, tolerance = tol)
  })

  test_that(paste(method, ": run_diff leaves significance unanswered"), {
    df <- diff_on(input_for(method), method)
    expect_true(all(is.na(df$is_significant)))
  })
}

# ---- the continuous engines only ----------------------------------------

for (method in continuous_engines()) {
  test_that(paste(method, ": a constant added to every value changes nothing"), {
    inp <- input_for(method)
    shifted <- inp
    shifted$expr_mat <- inp$expr_mat + 3
    a <- diff_on(inp, method)
    z <- aligned(diff_on(shifted, method), a$feature_id)
    expect_equal(z$effect, a$effect, tolerance = 1e-10)
    expect_equal(z$p_value, a$p_value, tolerance = 1e-10)
  })

  test_that(paste(method, ": scaling every value scales the effect and keeps p"), {
    inp <- input_for(method)
    scaled <- inp
    scaled$expr_mat <- inp$expr_mat * 2
    a <- diff_on(inp, method)
    z <- aligned(diff_on(scaled, method), a$feature_id)
    expect_equal(z$effect, 2 * a$effect, tolerance = 1e-10)
    expect_equal(z$p_value, a$p_value, tolerance = 1e-10)
  })
}

# ---- filtering ------------------------------------------------------------

test_that("filter_diff_results keeps exactly the rows inside the cutoffs", {
  df <- diff_on(realistic_input(), "ttest")
  for (pref in c("adjusted", "raw")) {
    col <- if (pref == "adjusted") "adj_p_value" else "p_value"
    for (p in c(0.001, 0.05, 0.5)) {
      kept <- filter_diff_results(df, p_cutoff = p, p_preference = pref)
      expect_setequal(kept$feature_id,
                      df$feature_id[!is.na(df[[col]]) & df[[col]] < p])
      expect_true(all(kept$is_significant))
      with_fc <- filter_diff_results(df, p_cutoff = p, p_preference = pref,
                                     effect_cutoff = 1)
      expect_setequal(with_fc$feature_id,
                      df$feature_id[!is.na(df[[col]]) & df[[col]] < p &
                                      !is.na(df$effect) & abs(df$effect) >= 1])
    }
  }
})

test_that("a looser cutoff admits a superset", {
  df <- diff_on(realistic_input(), "ttest")
  cuts <- c(0.001, 0.01, 0.05, 0.2, 1)
  sets <- lapply(cuts, function(p) filter_diff_results(df, p_cutoff = p)$feature_id)
  for (i in seq_along(cuts)[-1]) {
    expect_true(all(sets[[i - 1]] %in% sets[[i]]))
  }
  expect_setequal(sets[[length(sets)]], df$feature_id[!is.na(df$adj_p_value)])
})

# ---- figures --------------------------------------------------------------

# diff_significance() colours a point when p < threshold and, if an
# effect threshold is given, |effect| > it (strict on both sides; see
# plot-diff.R). filter_diff_results() uses >= for the effect. The
# figure and the table therefore agree exactly when no point sits on
# the effect boundary, which a continuous fixture guarantees.
test_that("the volcano colours exactly the points the thresholds admit", {
  b <- suppressMessages(run_diff(realistic_input(), method = "ttest",
                                 analysis_type = "group", group_col = "group",
                                 control_group = "G1", case_group = "G2"))
  df <- b$results$diff_result_df
  for (p in c(0.01, 0.05)) {
    for (e in list(NULL, 0.5)) {
      fig <- plot_volcano(b, p_threshold = p, effect_threshold = e)
      coloured <- sum(fig$data$.sig == "significant")
      expected <- nrow(filter_diff_results(df, p_cutoff = p,
                                           p_preference = "adjusted",
                                           effect_cutoff = e))
      expect_identical(coloured, expected, info = sprintf("p=%s e=%s", p, e))
    }
  }
})

test_that("with no threshold the volcano colours nothing", {
  b <- suppressMessages(run_diff(realistic_input(), method = "ttest",
                                 analysis_type = "group", group_col = "group",
                                 control_group = "G1", case_group = "G2"))
  fig <- plot_volcano(b, p_threshold = NULL)
  expect_identical(sum(fig$data$.sig == "significant"), 0L)
})

test_that("the MA plot agrees with the volcano about who is significant", {
  b <- suppressMessages(run_diff(realistic_input(), method = "ttest",
                                 analysis_type = "group", group_col = "group",
                                 control_group = "G1", case_group = "G2"))
  v <- plot_volcano(b, p_threshold = 0.05, effect_threshold = 0.5)
  m <- plot_ma(b, p_threshold = 0.05, effect_threshold = 0.5)
  expect_identical(v$data$feature_id[v$data$.sig == "significant"],
                   m$data$feature_id[m$data$.sig == "significant"])
})

# ---- enrichment ------------------------------------------------------------

test_that("GSEA ignores the feature thresholds ORA obeys", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("msigdbr")
  skip_if_not_installed("fgsea")
  b <- realistic_diff_bundle()
  gsea <- function(...) {
    set.seed(11)
    suppressWarnings(run_enrichment(b, type = "gsea", database = "hallmark",
                                    output_p_cutoff = 1, ...))$results$enrich_result_df
  }
  plain <- gsea()
  strict <- gsea(p_cutoff = 1e-9, effect_cutoff = 5, p_preference = "raw")
  strict <- strict[match(plain$pathway_id, strict$pathway_id), ]
  expect_equal(strict$effect, plain$effect)
  expect_equal(strict$p_value, plain$p_value)

  ora <- function(...) {
    run_enrichment(b, type = "ora", database = "hallmark", ...)$results$enrich_object[["both__hallmark"]]@gene
  }
  loose <- ora(p_cutoff = 0.5, effect_cutoff = 0)
  tight <- ora(p_cutoff = 0.5, effect_cutoff = 1)
  expect_true(all(tight %in% loose))
  expect_lt(length(tight), length(loose))
})

# ---- PCA ---------------------------------------------------------------------

test_that("constant features neither break the PCA nor change it", {
  inp <- realistic_input()
  mat <- inp$expr_mat
  padded <- rbind(mat,
                  ZERO = rep(0, ncol(mat)),
                  FLAT = rep(7.5, ncol(mat)))
  a <- pca_over_samples(mat)
  z <- pca_over_samples(padded)
  expect_identical(attr(a, "n_dropped"), 0L)
  expect_identical(attr(z, "n_dropped"), 2L)
  expect_equal(abs(z$x), abs(a$x), tolerance = 1e-10)
  expect_equal(z$sdev, a$sdev, tolerance = 1e-10)
})

test_that("too few varying features is an error that says how many", {
  mat <- matrix(1, 5, 4, dimnames = list(paste0("f", 1:5), paste0("s", 1:4)))
  mat[1, ] <- c(1, 2, 3, 4)
  expect_error(pca_over_samples(mat), "1 of 5")
})

test_that("applicable_diff_methods is decided by the assay, not the matrix", {
  counts <- realistic_input("rnaseq", n_per_group = 3L)
  offered <- applicable_diff_methods(counts)
  expect_true(all(c("deseq2", "edger") %in% offered))
  expect_false("limma" %in% offered)

  relabelled <- counts
  relabelled$expr_mat <- log2(counts$expr_mat + 1)
  relabelled$assay_type <- "logcpm"
  offered <- applicable_diff_methods(relabelled)
  expect_false(any(c("deseq2", "edger") %in% offered))
  expect_true("limma" %in% offered)
})
