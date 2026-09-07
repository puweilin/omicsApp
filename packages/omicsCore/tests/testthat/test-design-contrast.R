# Metadata column names and factor levels come from an uploaded file, so both
# are attacker-controlled. These tests cover the two ways that went wrong:
# the name reaching the R parser, and the label being matched against a
# coefficient by string.
#
# The security assertions deliberately require the call to SUCCEED as well as
# to leave no trace. Asserting only "an error came back" would have passed
# against the vulnerable code, which executed the payload and then failed on a
# design-dimension mismatch -- the marker file was already there.

# A fixture whose group column is named `col`, on the scale each backend wants.
fixture_named <- function(col, method = "limma", n_per_group = 6,
                          n_features = 24, levels = c("A", "B"),
                          extra_numeric = NULL) {
  set.seed(11)
  n <- n_per_group * length(levels)
  ids <- paste0("s", seq_len(n))
  fids <- paste0("g", seq_len(n_features))
  counts <- method %in% c("deseq2", "edger")
  expr <- if (counts) {
    matrix(rpois(n_features * n, 200), nrow = n_features, dimnames = list(fids, ids))
  } else {
    matrix(rnorm(n_features * n, 8, 1), nrow = n_features, dimnames = list(fids, ids))
  }
  # Signal, so a "did it actually fit" assertion has something to see.
  if (counts) expr[1:4, seq_len(n_per_group)] <- expr[1:4, seq_len(n_per_group)] * 3L
  else expr[1:4, seq_len(n_per_group)] <- expr[1:4, seq_len(n_per_group)] + 2

  meta <- data.frame(row.names = ids)
  meta[[col]] <- rep(levels, each = n_per_group)
  if (!is.null(extra_numeric)) {
    for (nm in extra_numeric) meta[[nm]] <- rnorm(n, 50, 8)
  }
  feat <- data.frame(feature_id = fids, feature_symbol = fids,
                     row.names = fids, stringsAsFactors = FALSE)
  omics_input(
    expr, meta, feat,
    omics_type = if (counts) "rnaseq" else "proteomics",
    assay_type = if (counts) "raw_count" else "normalized_intensity"
  )
}

payloads <- function(marker) {
  c(
    sprintf('file.create("%s")', marker),
    sprintf('I(file.create("%s"))+0', marker),
    sprintf('grp + (function(){file.create("%s"); 0})()', marker)
  )
}

test_that("a payload column name reaches the model and has no effect", {
  skip_if_not_installed("limma")
  for (method in c("limma", "lm", "ttest")) {
    for (payload in payloads("placeholder")) {
      marker <- tempfile("marker-")
      col <- sub("placeholder", marker, payload, fixed = TRUE)
      inp <- fixture_named(col, method = method)
      res <- run_diff(inp, method = method, analysis_type = "group",
                      group_col = col, control_group = "A", case_group = "B")
      # The model was actually fitted -- not refused before reaching it.
      expect_s3_class(res, "analysis_bundle")
      expect_equal(nrow(res$results$diff_result_df), 24L)
      expect_false(file.exists(marker),
                   info = paste(method, "executed:", col))
    }
  }
})

test_that("a payload covariate name has no effect either", {
  skip_if_not_installed("limma")
  marker <- tempfile("marker-")
  cov <- sprintf('file.create("%s")', marker)
  inp <- fixture_named("grp", extra_numeric = cov)
  res <- run_diff(inp, method = "limma", analysis_type = "group",
                  group_col = "grp", control_group = "A", case_group = "B",
                  covariates = cov)
  expect_s3_class(res, "analysis_bundle")
  expect_false(file.exists(marker))
})

test_that("a payload continuous column has no effect", {
  skip_if_not_installed("limma")
  marker <- tempfile("marker-")
  col <- sprintf('file.create("%s")', marker)
  inp <- fixture_named("grp", extra_numeric = col)
  for (m in c("limma", "lm")) {
    res <- run_diff(inp, method = m, analysis_type = "continuous",
                    continuous_col = col)
    expect_s3_class(res, "analysis_bundle")
  }
  expect_false(file.exists(marker))
})

# --- correctness -------------------------------------------------------

# limma fitted through a design built by hand, with the contrast given as
# positions rather than as a parsed string. Any disagreement with run_diff()
# means the production path resolved the contrast to a different coefficient.
reference_limma <- function(inp, group_col, control, case, covariates = NULL) {
  meta <- inp$meta_df
  g <- factor(meta[[group_col]], levels = c(control, case))
  d <- data.frame(g = g)
  if (!is.null(covariates)) for (i in seq_along(covariates)) {
    d[[paste0("cov", i)]] <- meta[[covariates[[i]]]]
  }
  design <- stats::model.matrix(~ 0 + ., data = d)
  v <- numeric(ncol(design))
  v[[match("gA", colnames(design))]] <- -1
  v[[match("gB", colnames(design))]] <- 1
  # The hand-built design names the two group columns gA / gB; the labels in
  # the fixture are always A and B, so this stays a fixed mapping.
  fit <- limma::eBayes(limma::contrasts.fit(
    limma::lmFit(inp$expr_mat, design),
    matrix(v, ncol = 1L, dimnames = list(colnames(design), "c"))))
  tt <- limma::topTable(fit, number = Inf, sort.by = "none")
  list(effect = tt$logFC, p = tt$P.Value)
}

test_that("group labels that collide with internal identifiers still contrast correctly", {
  skip_if_not_installed("limma")
  # `t001` and `L001` are what build_design() names its own term and level;
  # `v2` is what the first attempt at this fix named a covariate. A group is
  # allowed to be called any of them.
  for (case_label in c("B", "v2", "t001", "L001", "t001L002")) {
    inp <- fixture_named("grp", levels = c("A", case_label),
                         extra_numeric = "age")
    got <- run_diff(inp, method = "limma", analysis_type = "group",
                    group_col = "grp", control_group = "A",
                    case_group = case_label, covariates = "age")
    # Reference fit, with the same collision present in the data.
    meta <- inp$meta_df
    g <- factor(meta$grp, levels = c("A", case_label))
    d <- data.frame(g = g, cov1 = meta$age)
    design <- stats::model.matrix(~ 0 + ., data = d)
    v <- numeric(ncol(design))
    v[[match(paste0("g", "A"), colnames(design))]] <- -1
    v[[match(paste0("g", case_label), colnames(design))]] <- 1
    fit <- limma::eBayes(limma::contrasts.fit(
      limma::lmFit(inp$expr_mat, design),
      matrix(v, ncol = 1L, dimnames = list(colnames(design), "c"))))
    ref <- limma::topTable(fit, number = Inf, sort.by = "none")

    expect_equal(got$results$diff_result_df$effect, ref$logFC,
                 tolerance = 1e-10, info = case_label)
    expect_equal(got$results$diff_result_df$p_value, ref$P.Value,
                 tolerance = 1e-10, info = case_label)
  }
})

test_that("a covariate named after the group column does not steal the contrast", {
  skip_if_not_installed("limma")
  inp <- fixture_named("grp", extra_numeric = "grp_extra")
  got <- run_diff(inp, method = "limma", analysis_type = "group",
                  group_col = "grp", control_group = "A", case_group = "B",
                  covariates = "grp_extra")
  ref <- reference_limma(inp, "grp", "A", "B", covariates = "grp_extra")
  expect_equal(got$results$diff_result_df$effect, ref$effect, tolerance = 1e-10)
  expect_equal(got$results$diff_result_df$p_value, ref$p, tolerance = 1e-10)
})

test_that("ten or more model terms do not truncate a design column", {
  skip_if_not_installed("limma")
  # `t001` prefixed a `t0010` when the identifiers were v1..v10 and the
  # prefix was stripped with startsWith(): the tenth column became "0",
  # which makeContrasts() then rejected as a non-syntactic name.
  covs <- paste0("c", 1:12)
  inp <- fixture_named("grp", n_per_group = 20, extra_numeric = covs)
  design <- build_design(inp$meta_df, c("grp", covs),
                         primary = "grp", intercept = FALSE)
  expect_equal(anyDuplicated(colnames(design)), 0L)
  expect_true(all(nzchar(colnames(design))))
  expect_false(any(grepl("^[0-9]+$", colnames(design))))
  res <- run_diff(inp, method = "limma", analysis_type = "group",
                  group_col = "grp", control_group = "A", case_group = "B",
                  covariates = covs)
  expect_s3_class(res, "analysis_bundle")
})

test_that("spaces, CJK and punctuation in column names model correctly", {
  skip_if_not_installed("limma")
  # These are ordinary names in a real submission. Before the design was
  # built from generated identifiers they broke the formula outright, and in
  # lm they broke differently again: data.frame() ran make.names() over the
  # column while the pasted formula still asked for the original.
  for (col in c("年龄 (岁)", "batch-1", "treatment group", "a.b", "a|b")) {
    inp <- fixture_named(col)
    for (m in c("limma", "lm", "ttest")) {
      res <- run_diff(inp, method = m, analysis_type = "group",
                      group_col = col, control_group = "A", case_group = "B")
      expect_s3_class(res, "analysis_bundle")
      expect_false(anyNA(res$results$diff_result_df$p_value), info = paste(m, col))
    }
  }
})

test_that("group labels with spaces and CJK model correctly", {
  skip_if_not_installed("limma")
  for (lv in c("对照组", "G 1", "G-2")) {
    inp <- fixture_named("grp", levels = c("ctrl", lv))
    res <- run_diff(inp, method = "limma", analysis_type = "group",
                    group_col = "grp", control_group = "ctrl", case_group = lv)
    expect_s3_class(res, "analysis_bundle")
    expect_false(anyNA(res$results$diff_result_df$p_value), info = lv)
  }
})

test_that("build_design refuses to emit duplicate column names", {
  meta <- data.frame(a = rep(c("x", "y"), 4), b = rnorm(8))
  d <- build_design(meta, c("a", "b"), primary = "a", intercept = FALSE)
  expect_equal(anyDuplicated(colnames(d)), 0L)
  expect_identical(attr(d, "primary_id"), "t001")
  expect_identical(attr(d, "level_labels")$t001, c("x", "y"))
})

test_that("design_contrast places -1 and +1 on the requested levels only", {
  meta <- data.frame(a = rep(c("x", "y", "z"), each = 3), b = rnorm(9))
  d <- build_design(meta, c("a", "b"), primary = "a", intercept = FALSE)
  cm <- design_contrast(d, "x", "z")
  expect_equal(sum(cm), 0)
  expect_equal(cm[match("t001L001", rownames(cm)), 1L], -1)
  expect_equal(cm[match("t001L003", rownames(cm)), 1L], 1)
  expect_equal(cm[match("t002", rownames(cm)), 1L], 0)   # the covariate
})
