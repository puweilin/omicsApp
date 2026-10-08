# The remaining P1 items on the app side.

test_that("a result whose settings were changed afterwards says so until re-run", {
  proj <- shiny::reactiveVal(tutorial_project())
  shiny::testServer(diff_view_server, args = list(current_project = proj), {
    stale_note <- function() grepl("settings have changed",
                                   paste(as.character(output$notices), collapse = ""), fixed = TRUE)
    session$setInputs(layer = "proteomics", group_col = "group", control = "Control",
                      case = "TreatA", method = "limma", rerun = 1)
    expect_false(stale_note())
    session$setInputs(case = c("TreatA", "TreatB"))
    expect_true(stale_note())
    # Back to what was run: no longer out of date.
    session$setInputs(case = "TreatA")
    expect_false(stale_note())
    session$setInputs(covariates = "batch")
    expect_true(stale_note())
    session$setInputs(rerun = 2)
    expect_false(stale_note())
  })
})

test_that("a restored result is not flagged until a control is changed", {
  p <- tutorial_project()
  p$bundles <- list(diff = omicsCore::run_diff(p$experiments$proteomics, method = "limma",
                                               group_col = "group", control_group = "Control",
                                               case_group = "TreatA"))
  proj <- shiny::reactiveVal(NULL)
  shiny::testServer(diff_view_server, args = list(current_project = proj), {
    proj(p)
    session$flushReact()
    expect_false(is.null(diff_bundle()))
    html <- paste(as.character(output$notices), collapse = "")
    expect_false(grepl("settings have changed", html, fixed = TRUE))
    session$setInputs(method = "ttest")
    expect_match(paste(as.character(output$notices), collapse = ""), "settings have changed",
                 fixed = TRUE)
  })
})

test_that("an enrichment whose controls or thresholds moved says so until re-run", {
  skip_if_not(has_pkg("clusterProfiler"))
  p <- tutorial_project()
  d <- omicsCore::run_diff(p$experiments$proteomics, method = "limma", group_col = "group",
                           control_group = "Control", case_group = "TreatA")
  thr <- shiny::reactiveVal(list(p_cutoff = 0.05, p_preference = "adjusted", effect_cutoff = NULL))
  shiny::testServer(enrich_view_server,
                    args = list(diff_bundle = shiny::reactiveVal(d),
                                diff_thresholds = shiny::reactive(thr()),
                                current_project = shiny::reactiveVal(p)), {
    stale <- function() grepl("settings have changed",
                              paste(as.character(output$notices), collapse = ""), fixed = TRUE)
    session$setInputs(type = "ora", database = "hallmark", direction = "separate",
                      organism = "Hs", rerun = 1)
    expect_false(is.null(enrich_bundle()))
    expect_false(stale())
    session$setInputs(direction = "both")
    expect_true(stale())
    session$setInputs(direction = "separate")
    expect_false(stale())
    thr(list(p_cutoff = 0.01, p_preference = "adjusted", effect_cutoff = NULL))
    session$flushReact()
    expect_true(stale())
    session$setInputs(rerun = 2)
    expect_false(stale())
  })
})

# ---- several layers of one omics type ------------------------------------

p1_input <- function(fingerprint, shift = 0) {
  set.seed(7)
  mat <- matrix(stats::rnorm(40 * 6, 20, 0.3), 40,
                dimnames = list(paste0("G", 1:40), paste0("s", 1:6)))
  mat[1:8, 4:6] <- mat[1:8, 4:6] + 2 + shift
  meta <- data.frame(group = rep(c("A", "B"), each = 3), row.names = paste0("s", 1:6))
  omicsCore::omics_input(mat, meta, data.frame(feature_id = rownames(mat),
                                               feature_symbol = rownames(mat)),
                         omics_type = "proteomics", assay_type = "normalized_intensity",
                         source_fingerprint = fingerprint)
}

test_that("a second file of the same type can be kept beside the first", {
  proj <- shiny::reactiveVal(omicsCore::omics_project("P", list(proteomics = p1_input("fp-A"))))
  shiny::testServer(import_view_server, args = list(current_project = proj), {
    parsed(list(input = p1_input("fp-B"), report = NULL))
    session$setInputs(confirm = 1)
    expect_null(confirmed_input())          # asked first
    session$setInputs(confirm_keep_both = 1)
    expect_identical(confirmed_input()$layer_tag, "proteomics_2")
  })
})

test_that("a layer named in the Import view goes in under that name", {
  proj <- shiny::reactiveVal(omicsCore::omics_project("P", list(proteomics = p1_input("fp-A"))))
  shiny::testServer(import_view_server, args = list(current_project = proj), {
    parsed(list(input = p1_input("fp-B"), report = NULL))
    session$setInputs(layer_name = "Batch 2 (May)", confirm = 1)
    expect_identical(confirmed_input()$layer_tag, "Batch_2_May")
  })
  expect_identical(clean_layer_name("  ", "rnaseq"), "rnaseq")
  expect_identical(next_free_tag("proteomics", c("proteomics", "proteomics_2")), "proteomics_3")
})

test_that("results remember their layer when two layers share an omics type", {
  p <- omicsCore::omics_project("P", list(batch1 = p1_input("fp-A"),
                                          batch2 = p1_input("fp-B", shift = 1)))
  proj <- shiny::reactiveVal(p)
  shiny::testServer(diff_view_server, args = list(current_project = proj), {
    session$setInputs(layer = "batch2", group_col = "group", control = "A", case = "B",
                      method = "limma", rerun = 1)
    expect_identical(diff_bundle()$input_info$layer, "batch2")
    expect_identical(omicsCore::bundle_layer(p, diff_bundle()), "batch2")
  })
  d2 <- omicsCore::run_diff(p$experiments$batch2, method = "limma", group_col = "group",
                            control_group = "A", case_group = "B")
  d2$input_info$layer <- "batch2"
  p$bundles <- list(diff = d2)
  # Removing the other layer keeps this result; removing its own drops it.
  expect_identical(names(drop_layer_bundles(p$bundles, "batch1", p)), "diff")
  expect_length(drop_layer_bundles(p$bundles, "batch2", p), 0L)
  # Restored, it goes back on its own layer, not the first of its type.
  proj2 <- shiny::reactiveVal(NULL)
  shiny::testServer(diff_view_server, args = list(current_project = proj2), {
    proj2(p)
    session$flushReact()
    expect_identical(active()$tag, "batch2")
    expect_false(is.null(diff_bundle()))
  })
})

test_that("a saved QC result goes back to its own layer of two of the same type", {
  a <- p1_input("fp-A"); b <- p1_input("fp-B", shift = 1)
  q <- omicsCore::run_qc(b, outlier_method = "pca", impute_method = "none")
  q$input_info$layer <- "batch2"
  exps <- list(batch1 = a, batch2 = b)
  expect_identical(qc_bundle_layer(exps, q), "batch2")
  q$input_info$layer <- NULL
  expect_identical(qc_bundle_layer(exps, q), "batch1")   # older result: first that fits
})

# ---- figures --------------------------------------------------------------

test_that("the volcano labels its top features, and falls back to SVG without WebGL", {
  proj <- shiny::reactiveVal(tutorial_project())
  shiny::testServer(diff_view_server, args = list(current_project = proj), {
    session$setInputs(layer = "proteomics", group_col = "group", control = "Control",
                      case = "TreatA", method = "limma", rerun = 1)
    plain <- jsonlite::fromJSON(output$volcano, simplifyVector = FALSE)$x
    expect_length(plain$layout$annotations %||% list(), 0L)
    session$setInputs(label_top = TRUE)
    fig <- jsonlite::fromJSON(output$volcano, simplifyVector = FALSE)$x
    ann <- fig$layout$annotations
    expect_length(ann, 20L)
    top <- diff_bundle()$results$diff_result_df
    top <- top[order(top$adj_p_value), ][1, ]
    expect_identical(ann[[1]]$text, top$feature_symbol)
    expect_true(any(vapply(fig$data, function(t) identical(t$type, "scattergl"), logical(1))))
    session$rootScope()$setInputs(omics_webgl = FALSE)
    svg <- jsonlite::fromJSON(output$volcano, simplifyVector = FALSE)$x
    expect_false(any(vapply(svg$data, function(t) identical(t$type, "scattergl"), logical(1))))
    expect_match(svg$data[[1]]$text[[1]], "log2FC: ", fixed = TRUE)
  })
  expect_match(as.character(webgl_probe()), "omics_webgl", fixed = TRUE)
})
