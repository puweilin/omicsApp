# Directional ActivePathways.
#
# Brown's method merged the two layers' p-values with no regard for sign,
# so a pathway came back "significant" with no word on which way it went,
# and genes up in one layer and down in the other counted as much as
# genes up in both. The default is now the directional merge (DPM), with
# the layers expected to agree, and every pathway carries its direction.

skip_if_no_ap <- function() {
  skip_if_not_installed("ActivePathways")
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("msigdbr")
}

# One layer of the realistic fixture, with each named block shifted in
# the case group by the given amount.
ap_layer <- function(shifts, seed, prefix) {
  set.seed(seed)
  symbols <- realistic_symbols()
  n <- length(symbols)
  grp <- rep(c("G1", "G2"), each = 6L)
  m <- matrix(stats::rnorm(n * 12L, sd = 0.6), n, 12L) + stats::rnorm(n, 20, 2)
  for (blk in names(shifts)) {
    hit <- symbols %in% REAL_GENE_SETS[[blk]]
    m[hit, grp == "G2"] <- m[hit, grp == "G2"] + shifts[[blk]]
  }
  ids <- paste0(prefix, "_", symbols)
  sids <- sprintf("%s%02d", prefix, 1:12)
  dimnames(m) <- list(ids, sids)
  omics_input(m, data.frame(group = grp, row.names = sids),
              data.frame(feature_id = ids, feature_symbol = symbols),
              omics_type = "proteomics", assay_type = "normalized_intensity")
}

# G2M goes up in both layers; ADIPO goes up in the first and by
# `adipo_b` in the second.
ap_run <- function(adipo_b, ...) {
  a <- ap_layer(list(G2M = 1.5, ADIPO = 1.5, OXPHOS = -1.5), 1, "A")
  b <- ap_layer(list(G2M = 1.5, ADIPO = adipo_b, OXPHOS = -1.5), 2, "B")
  proj <- omics_project("x", list(a = a, b = b))
  diffs <- lapply(list(a = a, b = b), function(x) {
    run_diff(x, method = "limma", group_col = "group",
             control_group = "G1", case_group = "G2")
  })
  suppressMessages(run_integration(proj, method = "active_pathways",
                                   experiments = c("a", "b"), diff_bundles = diffs,
                                   database = "hallmark", ...))
}

row_of <- function(bundle, pathway) {
  df <- bundle$results$integration_df
  df[df$feature_symbol == pathway, , drop = FALSE]
}

test_that("a pathway driven by genes up in both layers is reported as up", {
  skip_if_no_ap()
  r <- ap_run(1.5)
  expect_identical(r$params$merge_method, "DPM")
  expect_equal(r$params$constraints_vector, c(1, 1))
  g2m <- row_of(r, "HALLMARK_G2M_CHECKPOINT")
  expect_equal(nrow(g2m), 1L)
  expect_identical(g2m$direction, "up")
  expect_identical(g2m$direction_a, "up")
  expect_identical(g2m$direction_b, "up")
  expect_true(g2m$layers_agree)
  expect_gt(g2m$n_genes_agree, 10L)
  expect_identical(g2m$n_genes_disagree, 0L)
  expect_identical(g2m$evidence, "shared")
  ox <- row_of(r, "HALLMARK_OXIDATIVE_PHOSPHORYLATION")
  if (nrow(ox)) expect_identical(ox$direction, "down")
  check_integration_result_schema(r$results$integration_df)
})

test_that("layers moving opposite ways are penalised against layers that agree", {
  skip_if_no_ap()
  agree <- ap_run(1.5)
  conflict <- ap_run(-1.5)
  adipo <- "HALLMARK_ADIPOGENESIS"
  p_agree <- row_of(agree, adipo)$adj_p_value
  p_conflict <- row_of(conflict, adipo)$adj_p_value
  expect_length(p_agree, 1L)
  expect_lt(p_agree, 0.05)
  # Gone from the significant pathways, or at least much weaker.
  expect_true(length(p_conflict) == 0L || p_conflict > p_agree)

  # It is the direction doing it: Brown's method, blind to sign, finds
  # the conflicting pathway as strongly -- and still says it is mixed.
  brown <- ap_run(-1.5, merge_method = "Brown")
  expect_identical(brown$params$merge_method, "Brown")
  expect_null(brown$params$constraints_vector)
  b_row <- row_of(brown, adipo)
  expect_lt(b_row$adj_p_value, 0.05)
  expect_identical(b_row$direction, "mixed")
  expect_false(b_row$layers_agree)
  expect_gt(b_row$n_genes_disagree, 10L)
})

test_that("an ActivePathways without the directional method falls back to Brown's, and says so", {
  skip_if_no_ap()
  testthat::local_mocked_bindings(active_pathways_directional = function() FALSE)
  r <- ap_run(1.5)
  expect_identical(r$params$merge_method, "Brown")
  expect_match(r$warnings, "Brown's method", all = FALSE)
  expect_identical(row_of(r, "HALLMARK_G2M_CHECKPOINT")$direction, "up")
})

test_that("a bad constraints vector is refused in plain words", {
  skip_if_no_ap()
  expect_error(ap_run(1.5, constraints_vector = c(2, 1)), "constraints_vector")
})

test_that("the dot plot colours by direction; an older result still plots", {
  df <- data.frame(
    feature_id = c("P1", "P2", "P3"),
    feature_symbol = c("HALLMARK_ONE", "HALLMARK_TWO", "HALLMARK_THREE"),
    result_type = "active_pathways", experiments = "a vs b", comparison = "x | y",
    effect = c(10, 6, 3), effect_type = "neg_log10_padj",
    statistic = c(1e-10, 1e-6, 1e-3), statistic_type = "adjusted_p_val",
    p_value = c(1e-10, 1e-6, 1e-3), adj_p_value = c(1e-10, 1e-6, 1e-3),
    direction = c("up", "down", "mixed"), quadrant = "a,b",
    is_significant = TRUE, source_label = "ap",
    evidence = c("shared", "unique", "combined"), stringsAsFactors = FALSE)
  b <- new_analysis_bundle("run_integration",
                           params = list(method = "active_pathways",
                                         experiments = c("a", "b"),
                                         merge_method = "DPM"),
                           results = list(integration_df = df))
  p <- plot_integration(b, view = "dotplot")
  built <- ggplot2::ggplot_build(p)
  expect_setequal(built$data[[1L]]$colour,
                  c(omics_colors$up, omics_colors$down, omics_colors$conc_up_down))
  expect_identical(p$scales$get_scales("colour")$name, "direction")
  expect_match(p$labels$subtitle, "agree in direction", fixed = TRUE)

  # Before the directional rework: evidence in `direction`, no direction.
  old <- df
  old$evidence <- NULL
  old$direction <- c("shared", "unique", "combined")
  b$results$integration_df <- old
  b$params$merge_method <- "Brown"
  p_old <- plot_integration(b, view = "dotplot")
  expect_identical(p_old$scales$get_scales("colour")$name, "adj p")
  expect_s3_class(ggplot2::ggplot_build(p_old), "ggplot_built")
})

test_that("the script repeats the directional method", {
  inp <- ap_layer(list(G2M = 1), 1, "A")
  proj <- omics_project("p", list(a = inp, b = inp))
  ib <- new_analysis_bundle(
    "run_integration",
    params = list(method = "active_pathways", experiments = c("a", "b"),
                  by = "feature_symbol", database = "hallmark",
                  organism = "Homo sapiens", p_preference = "raw",
                  significant = 0.05, geneset_filter = c(5L, 1000L),
                  merge_method = "DPM", constraints_vector = c(1, 1)))
  proj$bundles <- list(integration = ib)
  txt <- paste(export_script(proj, include_plots = FALSE), collapse = "\n")
  expect_match(txt, 'merge_method       = "DPM"', fixed = TRUE)
  expect_match(txt, "constraints_vector = c(1, 1)", fixed = TRUE)
  # A result made with Brown's method is repeated with it.
  ib$params$merge_method <- "Brown"
  ib$params$constraints_vector <- NULL
  proj$bundles$integration <- ib
  txt <- paste(export_script(proj, include_plots = FALSE), collapse = "\n")
  expect_match(txt, 'merge_method   = "Brown"', fixed = TRUE)
  expect_false(grepl("constraints_vector", txt, fixed = TRUE))
})
