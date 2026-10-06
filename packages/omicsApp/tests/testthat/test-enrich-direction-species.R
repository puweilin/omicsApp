# The Enrichment view's Direction and Species controls.
#
# ORA used to pool up- and down-regulated hits by default; the view now
# tests them as separate lists unless told otherwise, says so in plain
# words, and shows which list each pathway came from. The species menu
# offers every species msigdbr carries the gene sets to.

ora_frame <- function() {
  data.frame(
    database = "hallmark", result_type = "ora", comparison = "G2_vs_G1",
    pathway_id = c("HALLMARK_G2M_CHECKPOINT", "HALLMARK_ADIPOGENESIS"),
    pathway_name = c("G2M CHECKPOINT", "ADIPOGENESIS"),
    effect = NA_real_, effect_type = NA_character_,
    direction = c("up", "down"),
    p_value = c(1e-8, 1e-3), adj_p_value = c(1e-7, 1e-2), q_value = NA_real_,
    gene_set_size = 40, overlap_size = c(30, 8),
    overlap_features = "A", leading_features = "A",
    source_label = "ora", stringsAsFactors = FALSE)
}

ora_bundle <- function(direction = "separate", df = ora_frame(), warnings = character(0)) {
  omicsCore::new_analysis_bundle(
    "run_enrichment", input_info = list(omics_type = "proteomics"),
    params = list(type = "ora", database = "hallmark", organism = "Homo sapiens",
                  direction = direction, comparison = "G2_vs_G1"),
    results = list(enrich_result_df = df),
    warnings = warnings)
}

test_that("the Direction control defaults to separate lists, in plain words", {
  html <- as.character(enrich_params_card(shiny::NS("e")))
  expect_match(html, "Up and down separately", fixed = TRUE)
  expect_match(html, "Up only", fixed = TRUE)
  expect_match(html, "Down only", fixed = TRUE)
  expect_match(html, "Up and down pooled", fixed = TRUE)
  expect_match(html, 'value="separate" checked', fixed = TRUE)
  # "both" keeps its old meaning, so a saved project restores to it.
  expect_identical(unname(ENRICH_DIRECTION_CHOICES[["Up and down pooled"]]), "both")
})

test_that("the table says which gene list each ORA pathway was found in", {
  b <- ora_bundle()
  html <- NULL
  shiny::testServer(enrich_view_server, args = list(diff_bundle = shiny::reactiveVal(NULL)), {
    enrich_bundle(b); is_demo(FALSE)
    session$setInputs(show_p = "adjusted", show_cutoff = 0.05)
    session$flushReact()
    html <<- paste(unlist(output$hits), collapse = " ")
  })
  expect_match(html, "Found among", fixed = TRUE)
  # The cells travel to the browser separately; the rows are built here.
  tab <- enrich_hits_table(ora_frame(), "adjusted", 0.05)
  expect_identical(tab$Pathway, c("G2M CHECKPOINT", "ADIPOGENESIS"))
  expect_identical(tab$`Found among`, c("up-regulated genes", "down-regulated genes"))
  expect_false("NES" %in% names(tab))
})

test_that("a pooled result has no direction column", {
  df <- ora_frame()
  df$direction <- NA_character_
  b <- ora_bundle("both", df)
  html <- NULL
  shiny::testServer(enrich_view_server, args = list(diff_bundle = shiny::reactiveVal(NULL)), {
    enrich_bundle(b); is_demo(FALSE)
    session$setInputs(show_p = "adjusted", show_cutoff = 0.05)
    session$flushReact()
    html <<- paste(unlist(output$hits), collapse = " ")
  })
  expect_false(grepl("Found among", html, fixed = TRUE))
  tab <- enrich_hits_table(df, "adjusted", 0.05)
  expect_false(any(c("Found among", "Direction") %in% names(tab)))
  expect_identical(nrow(tab), 2L)
  # GSEA keeps its sign as Direction.
  g <- ora_frame()
  g$effect <- c(2.1, -1.4)
  g$effect_type <- "nes"
  expect_identical(enrich_hits_table(g, "adjusted", 0.05)$Direction, c("up", "down"))
})

test_that("a note from the run is shown with the result", {
  b <- ora_bundle(warnings = "The gene names are written differently from the Mus musculus gene sets.")
  html <- NULL
  shiny::testServer(enrich_view_server, args = list(diff_bundle = shiny::reactiveVal(NULL)), {
    enrich_bundle(b); is_demo(FALSE)
    session$flushReact()
    html <<- paste(unlist(output$notices), collapse = " ")
  })
  expect_match(html, "written differently", fixed = TRUE)
})

test_that("a live ORA run tests up and down separately unless asked to pool", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("msigdbr")
  diff_bundle <- shiny::reactiveVal(omicsApp:::example_diff_bundle())
  shiny::testServer(enrich_view_server, args = list(diff_bundle = diff_bundle), {
    session$setInputs(type = "ora", database = "hallmark", rerun = 1)
    expect_identical(enrich_bundle()$params$direction, "separate")
    session$setInputs(direction = "both", rerun = 2)
    expect_identical(enrich_bundle()$params$direction, "both")
  })
})

test_that("a result with no direction is enriched pooled whatever the control says", {
  b <- list(results = list(diff_result_df = data.frame(analysis_type = "continuous_spline")))
  expect_true(diff_undirected(b))
  b$results$diff_result_df$analysis_type <- "group"
  expect_false(diff_undirected(b))
  expect_false(diff_undirected(NULL))
})

test_that("the species menu offers every species omicsCore lists", {
  ch <- species_choices()
  expect_true(all(c("Hs", "Mm", "Rn", "Dr", "Dm", "Sc") %in% ch))
  expect_identical(names(ch)[ch == "Rn"], "Rat")
  expect_identical(species_code("Rattus norvegicus"), "Rn")
  expect_identical(species_code("Mm"), "Mm")
  expect_null(species_code("Bos taurus"))
  expect_null(species_code(NULL))
  expect_identical(species_label("Dr"), "Zebrafish")
  expect_identical(species_label("Sc"), "Yeast")
})

test_that("the species guess reads zebrafish, fly, worm and yeast names", {
  fish <- c("tp53", "egfra", "myca", "il6", "tnfa", "cd4-1", "actb1", "gapdh",
            "stat3", "jun", "fosab")
  expect_identical(guess_organism(fish), "Dr")
  fly <- c(paste0("CG", c(1824, 2341, 9876, 4512, 7788)), "Mcad", "yip2",
           "Arc42", "ATPCL", "mAcon1", "p53")
  expect_identical(guess_organism(fly), "Dm")
  worm <- c("unc-54", "cep-1", "daf-2", "daf-16", "lin-4", "let-7", "abt-2",
            "haf-6", "acdh-7", "aco-2", "egl-1")
  expect_identical(guess_organism(worm), "Ce")
  yeast <- c("YAL001C", "YAL002W", "YBR123C", "YCL004W", "YDR101C", "YER005W",
             "ACT1", "CDC28", "RAD9", "YGL123W", "YHR001W")
  expect_identical(guess_organism(yeast), "Sc")
  # Unchanged for the two it already told apart.
  expect_identical(guess_organism(c("Trp53", "Egfr", "Myc", "Il6", "Tnf", "Cd4", "Actb",
                                    "Gapdh", "Stat3", "Jun", "Fos")), "Mm")
  expect_identical(guess_organism(c("TP53", "EGFR", "MYC", "IL6", "TNF", "CD4", "ACTB",
                                    "GAPDH", "STAT3", "JUN", "FOS")), "Hs")
})

# ---- integration: directional ActivePathways ----------------------------

ap_frame <- function(directional = TRUE) {
  df <- data.frame(
    feature_id = c("P1", "P2", "P3"),
    feature_symbol = c("HALLMARK_ONE", "HALLMARK_TWO", "HALLMARK_THREE"),
    result_type = "active_pathways", experiments = "a vs b", comparison = "x | y",
    effect = c(10, 6, 3), effect_type = "neg_log10_padj",
    statistic = c(1e-10, 1e-6, 1e-3), statistic_type = "adjusted_p_val",
    p_value = c(1e-10, 1e-6, 1e-3), adj_p_value = c(1e-10, 1e-6, 1e-3),
    direction = c("up", "down", "mixed"), quadrant = "a,b",
    is_significant = TRUE, source_label = "ap",
    direction_a = c("up", "down", "up"), direction_b = c("up", "down", "down"),
    layers_agree = c(TRUE, TRUE, FALSE),
    evidence = c("shared", "unique", "combined"), stringsAsFactors = FALSE)
  if (!directional) {
    df$direction <- df$evidence
    df[c("evidence", "direction_a", "direction_b", "layers_agree")] <- NULL
  }
  df
}

ap_bundle <- function(df) {
  omicsCore::new_analysis_bundle(
    "run_integration",
    params = list(method = "active_pathways", experiments = c("prot", "rna"),
                  merge_method = "DPM"),
    results = list(integration_df = df))
}

test_that("ActivePathways cards and table say which way each pathway went", {
  b <- ap_bundle(ap_frame())
  html <- as.character(integration_stat_cards(b$results$integration_df, b))
  expect_match(html, "Up in both layers", fixed = TRUE)
  expect_match(html, "Down in both layers", fixed = TRUE)
  expect_match(html, "Layers disagree", fixed = TRUE)
  tab <- integration_result_table(b$results$integration_df, "active_pathways",
                                  c("prot", "rna"))
  expect_identical(tab$Direction, c("up in both layers", "down in both layers",
                                    "mixed / layers disagree"))
  expect_true(all(c("In prot", "In rna", "Found by") %in% names(tab)))
  expect_identical(tab$`Found by`, c("both layers", "one layer", "only combined"))
})

test_that("an ActivePathways result saved before directions still shows", {
  b <- ap_bundle(ap_frame(directional = FALSE))
  html <- as.character(integration_stat_cards(b$results$integration_df, b))
  expect_match(html, "Shared by both", fixed = TRUE)
  tab <- integration_result_table(b$results$integration_df, "active_pathways",
                                  c("prot", "rna"))
  expect_identical(tab$Direction, c("shared", "unique", "combined"))
})
