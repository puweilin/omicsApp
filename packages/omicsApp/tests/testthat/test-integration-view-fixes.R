# The Integration view's defects from the 2026-09 review.

iv_rendered <- function(x) !is.null(x) && nzchar(paste(as.character(x), collapse = ""))

iv_input <- function(prefix, omics, groups = rep(c("ctrl", "A", "B"), each = 4),
                     donor = TRUE, symbols = NULL, shift = 2) {
  set.seed(if (omics == "proteomics") 1 else 2)
  n <- 30L
  sym <- symbols %||% c("TP53", "EGFR", "MYC", "AKT1", "PTEN", paste0("G", 1:25))
  ids <- paste0(prefix, seq_len(n))
  samp <- paste0(prefix, "_", seq_along(groups))
  m <- matrix(stats::rnorm(n * length(groups), 8, 0.3), n,
              dimnames = list(ids, samp))
  m[1:5, groups == "A"] <- m[1:5, groups == "A"] + shift
  m[6:10, groups == "B"] <- m[6:10, groups == "B"] - shift
  meta <- data.frame(group = groups, age = seq_along(groups) + 30,
                     row.names = samp, stringsAsFactors = FALSE)
  if (donor) meta$donor <- paste0("D", seq_along(groups))
  omicsCore::omics_input(
    m, meta,
    data.frame(feature_id = ids, feature_symbol = sym, row.names = ids,
               stringsAsFactors = FALSE),
    omics_type = omics,
    assay_type = if (omics == "proteomics") "normalized_intensity" else "logcpm")
}

iv_project <- function(...) {
  omicsCore::omics_project("iv", list(proteomics = iv_input("p", "proteomics", ...),
                                      rnaseq = iv_input("r", "rnaseq", ...)))
}

iv_diff <- function(proj, tag = "rnaseq", cases = "A") {
  b <- omicsCore::run_diff(proj$experiments[[tag]], method = "limma",
                           group_col = "group", control_group = "ctrl",
                           case_group = cases)
  if (length(cases) > 1L) b <- omicsCore::select_comparison(b, paste0(cases[1], "_vs_ctrl"))
  b
}

test_that("the primary layer is the one the diff ran on, not a guess from its omics type", {
  proj <- iv_project()
  shiny::testServer(integration_view_server, args = list(
    current_project = shiny::reactiveVal(proj),
    diff_bundle = shiny::reactiveVal(iv_diff(proj, "rnaseq")),
    diff_layer = shiny::reactiveVal("rnaseq")), {
    session$flushReact()
    expect_identical(layers()$primary, "rnaseq")
    expect_identical(layers()$partner, "proteomics")
    b <- integration_bundle()
    expect_s3_class(b, "analysis_bundle")
    expect_identical(b$params$experiments, c("rnaseq", "proteomics"))
  })
})

test_that("the partner is repeated with the primary's covariates and the Differential thresholds", {
  proj <- iv_project()
  primary <- omicsCore::run_diff(proj$experiments$proteomics, method = "limma",
                                 group_col = "group", control_group = "ctrl",
                                 case_group = "A", covariates = "age")
  shiny::testServer(integration_view_server, args = list(
    current_project = shiny::reactiveVal(proj),
    diff_bundle = shiny::reactiveVal(primary),
    diff_layer = shiny::reactiveVal("proteomics"),
    diff_thresholds = shiny::reactiveVal(list(p_cutoff = 0.1, p_preference = "raw",
                                              effect_cutoff = 1))), {
    session$flushReact()
    b <- integration_bundle()
    expect_identical(b$params$p_cutoff, 0.1)
    expect_identical(b$params$p_preference, "raw")
    expect_identical(b$params$effect_cutoff, 1)
    expect_identical(sec_cache()$bundle$params$covariates, "age")
  })
})

test_that("a multi-contrast diff is repeated with the same groups in the partner's model", {
  proj <- iv_project()
  primary <- iv_diff(proj, "proteomics", cases = c("A", "B"))
  shiny::testServer(integration_view_server, args = list(
    current_project = shiny::reactiveVal(proj),
    diff_bundle = shiny::reactiveVal(primary),
    diff_layer = shiny::reactiveVal("proteomics")), {
    session$flushReact()
    expect_identical(can_run()$all_cases, c("A", "B"))
    sec <- sec_cache()$bundle
    expect_identical(sec$params$comparison, "A_vs_ctrl")
    expect_identical(sec$params$all_case_groups, c("A", "B"))
    df <- integration_bundle()$results$integration_df
    expect_match(unique(df$comparison), "A_vs_ctrl | A_vs_ctrl", fixed = TRUE)
  })
})

test_that("filing the result into the project does not start a second run", {
  p0 <- iv_project()
  proj <- shiny::reactiveVal(p0)
  shiny::testServer(integration_view_server, args = list(
    current_project = proj,
    diff_bundle = shiny::reactiveVal(iv_diff(p0)),
    diff_layer = shiny::reactiveVal("rnaseq")), {
    session$flushReact()
    expect_identical(runs$token, 1L)
    # What app_server() does with the result.
    p <- current_project()
    p$bundles <- list(integration = integration_bundle())
    current_project(p)
    session$flushReact()
    expect_identical(runs$token, 1L)
    # Re-run is still a request to compute.
    session$setInputs(rerun = 1)
    expect_identical(runs$token, 2L)
  })
})

test_that("changing only the thresholds reuses the partner's diff", {
  proj <- iv_project()
  th <- shiny::reactiveVal(list(p_cutoff = 0.05, p_preference = "adjusted", effect_cutoff = 0))
  shiny::testServer(integration_view_server, args = list(
    current_project = shiny::reactiveVal(proj),
    diff_bundle = shiny::reactiveVal(iv_diff(proj)),
    diff_layer = shiny::reactiveVal("rnaseq"),
    diff_thresholds = th), {
    session$flushReact()
    first <- sec_cache()$bundle
    th(list(p_cutoff = 0.2, p_preference = "adjusted", effect_cutoff = 0))
    session$flushReact()
    expect_identical(runs$token, 2L)
    expect_identical(sec_cache()$bundle, first)
    expect_identical(integration_bundle()$params$p_cutoff, 0.2)
  })
})

test_that("sample-level correlation runs from the view on a donor pairing", {
  proj <- iv_project()
  shiny::testServer(integration_view_server, args = list(
    current_project = shiny::reactiveVal(proj)), {
    session$setInputs(method = "correlation")
    session$flushReact()
    expect_true(isTRUE(can_run()$ok))
    b <- integration_bundle()
    expect_identical(b$params$method, "correlation")
    expect_identical(b$params$method_info$pairing_source, "donor")
    expect_true(iv_rendered(output$stats))
  })
})

test_that("a failed run shows the error, not the demo", {
  # No shared symbols: concordance has nothing to join on.
  proj <- omicsCore::omics_project("x", list(
    proteomics = iv_input("p", "proteomics", symbols = paste0("PX", 1:30)),
    rnaseq = iv_input("r", "rnaseq", symbols = paste0("RX", 1:30))))
  shiny::testServer(integration_view_server, args = list(
    current_project = shiny::reactiveVal(proj),
    diff_bundle = shiny::reactiveVal(iv_diff(proj)),
    diff_layer = shiny::reactiveVal("rnaseq"),
    diff_thresholds = shiny::reactiveVal(NULL)), {
    session$setInputs(method = "correlation")
    session$flushReact()
    expect_false(isTRUE(is_demo()))
    expect_null(integration_bundle())
    expect_false(is.null(integration_error()))
    html <- paste(as.character(output$notices), collapse = "")
    expect_match(html, "notice-error", fixed = TRUE)
    expect_match(html, "gene symbols in common", fixed = TRUE)
  })
})

test_that("a partner without the contrast's levels names what is missing", {
  proj <- iv_project()
  sec <- proj$experiments$proteomics
  sec$meta_df$group[sec$meta_df$group == "A"] <- "ctrl"
  proj$experiments$proteomics <- sec
  shiny::testServer(integration_view_server, args = list(
    current_project = shiny::reactiveVal(proj),
    diff_bundle = shiny::reactiveVal(iv_diff(proj)),
    diff_layer = shiny::reactiveVal("rnaseq")), {
    session$flushReact()
    expect_false(isTRUE(can_run()$ok))
    expect_match(can_run()$detail, "'A'", fixed = TRUE)
    expect_match(paste(as.character(output$notices), collapse = ""), "no samples in")
  })
})

test_that("accepting a pairing keeps the pairings saved for other layers", {
  mk <- function(ids) {
    m <- matrix(seq_len(4 * length(ids)) * 1000, 4, dimnames = list(paste0("F", 1:4), ids))
    omicsCore::omics_input(m, data.frame(g = rep("a", length(ids)), row.names = ids),
                           data.frame(feature_id = paste0("F", 1:4), row.names = paste0("F", 1:4)),
                           omics_type = "proteomics", assay_type = "raw_intensity")
  }
  proj <- omicsCore::omics_project("p", list(
    prot = mk(c("RD001-C", "RD002-C")), rna = mk(c("RD001_F", "RD002_F")),
    metab = mk(c("M1", "M2"))))
  proj$sample_link <- data.frame(tag = c("metab", "metab"), sample_id = c("M1", "M2"),
                                 donor_id = c("RD001", "RD002"), stringsAsFactors = FALSE)
  # The saved link does not cover prot/rna, so their pairing is still a guess.
  shiny::testServer(integration_view_server, args = list(
    current_project = shiny::reactiveVal(proj)), {
    session$flushReact()
    expect_identical(pairing()$source, "suggested")
    session$setInputs(accept_pairing = 1)
    link <- current_project()$sample_link
    expect_setequal(unique(link$tag), c("metab", "prot", "rna"))
    expect_equal(nrow(link), 6L)
  })
})

test_that("the header buttons that did nothing are gone", {
  html <- paste(as.character(integration_view_ui("i")), collapse = "")
  expect_false(grepl(">Mapping<", html, fixed = TRUE))
  expect_false(grepl("Export integration bundle", html, fixed = TRUE))
})
