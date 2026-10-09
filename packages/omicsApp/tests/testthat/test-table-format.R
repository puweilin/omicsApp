# p-values in the on-screen tables (R/table_format.R): written as text,
# three significant digits and scientific below 0.001, and sorted by a
# hidden copy of the number. The browser printed signif(1.8e-248, 3) as
# "1.799999999999999e-248", which also pushed the Enriched sets table's
# last column out of its card.

test_that("p-values are written with three significant digits", {
  expect_identical(
    format_p_value(c(1.804461e-248, 4.41504e-16, 0.0123456, 0.5, 1, NA, 0.001, 0.00099996, 0)),
    c("1.80e-248", "4.42e-16", "0.0123", "0.5", "1", "\u2013", "0.001", "0.001", "0"))
  expect_identical(format_p_value(numeric(0)), character(0))
})

test_that("text columns sort by their numbers, hidden at the end of the table", {
  out <- data.frame(Name = c("a", "b"), p = c(1e-5, 0.02), Overlap = c("3/10", "12/40"),
                    stringsAsFactors = FALSE)
  tab <- dt_sortable_text(out, p_cols = "p", sort_keys = list(Overlap = c(3, 12)))
  expect_identical(names(tab$data), c("Name", "p", "Overlap", ".sort_1", ".sort_2"))
  expect_identical(tab$data$p, c("1.00e-05", "0.02"))
  expect_identical(tab$data$.sort_1, c(1e-5, 0.02))
  expect_identical(tab$data$.sort_2, c(3, 12))
  defs <- tab$column_defs
  # p (column 1, from 0) sorts by column 3, Overlap by column 4.
  expect_identical(defs[[1]][c("targets", "orderData")], list(targets = 1L, orderData = 3L))
  expect_identical(defs[[2]][c("targets", "orderData")], list(targets = 2L, orderData = 4L))
  expect_identical(defs[[1]]$className, "dt-right")
  expect_identical(defs[[3]], list(targets = 3:4, visible = FALSE, searchable = FALSE))
  # Nothing to do: the frame unchanged, no definitions.
  expect_identical(dt_sortable_text(out[1], p_cols = "p"), list(data = out[1], column_defs = list()))
})

enrich_frame <- function(type = c("ora", "gsea")) {
  type <- match.arg(type)
  df <- data.frame(
    database = "hallmark", result_type = type, comparison = "B_vs_A",
    pathway_id = paste0("H", 1:3),
    pathway_name = c("INFLAMMATORY RESPONSE", "TNFA SIGNALING VIA NFKB", "COMPLEMENT"),
    effect = NA_real_, effect_type = NA_character_, direction = c("up", "up", "down"),
    p_value = c(4.4e-250, 2.2e-17, 0.004), adj_p_value = c(1.804461e-248, 4.41504e-16, 0.0123456),
    q_value = NA_real_, gene_set_size = c(152, 156, 150), overlap_size = c(151, 37, 16),
    overlap_features = "A", leading_features = NA_character_, source_label = type,
    stringsAsFactors = FALSE)
  if (type == "gsea") {
    df$effect <- c(2.61, 1.9, -1.75)
    df$effect_type <- "nes"
    df$overlap_size <- NA_real_
  }
  df
}

test_that("the Enriched sets table writes its numbers out and keeps what each method has", {
  tab <- enrich_hits_table(enrich_frame("ora"), "adjusted", 0.05)
  d <- tab$data
  shown <- names(d)[!startsWith(names(d), ".sort_")]
  expect_identical(shown, c("Pathway", "Gene list", "adjusted p", "Overlap"))
  expect_identical(d$`adjusted p`, c("1.80e-248", "4.42e-16", "0.0123"))
  expect_identical(d$`Gene list`, c("up", "up", "down"))
  expect_identical(d$Overlap, c("151/152", "37/156", "16/150"))
  # Sorted by the numbers, not the text: "16/150" would sort before "37/156".
  expect_identical(d$.sort_1, c(1.804461e-248, 4.41504e-16, 0.0123456))
  expect_identical(d$.sort_2, c(151, 37, 16))

  # GSEA: its NES, and the pathway's size (it counts no overlap); its
  # direction is the NES's sign.
  g <- enrich_hits_table(enrich_frame("gsea"), "raw", 0.05)$data
  expect_identical(names(g)[!startsWith(names(g), ".sort_")],
                   c("Pathway", "NES", "p", "Set size"))
  expect_identical(g$NES, c("+2.61", "+1.90", "-1.75"))
  expect_identical(g$p, c("4.40e-250", "2.20e-17", "0.004"))
  expect_identical(g$`Set size`, c(152L, 156L, 150L))
  # The raw p column is the one that sorts by its number.
  expect_identical(g$.sort_1, c(4.4e-250, 2.2e-17, 0.004))
})

test_that("the Enriched sets widget fits its card and holds the pathway column on a phone", {
  w <- enrich_hits_datatable(enrich_hits_table(enrich_frame("ora"), "adjusted", 0.05))
  opts <- w$x$options
  expect_true(opts$scrollX)
  expect_identical(opts$fixedColumns, list(left = 1))
  expect_true("FixedColumns" %in% unlist(w$x$extensions))
  defs <- opts$columnDefs
  right <- Filter(function(d) identical(d$className, "dt-right") && is.null(d$orderData), defs)[[1]]
  expect_identical(right$targets, 2:3)  # adjusted p, Overlap
  expect_true(any(vapply(defs, function(d) identical(d$visible, FALSE), logical(1))))
  # The gene list in the app's colours: up red, down blue.
  expect_match(as.character(opts$rowCallback), omics_colors$up, fixed = TRUE)
  expect_match(as.character(opts$rowCallback), omics_colors$down, fixed = TRUE)
  # The card's tighter cells, scoped to this table.
  css <- as.character(enrich_hits_card(shiny::NS("enrich")))
  expect_match(css, "#enrich-hits table.dataTable>thead>tr>th{padding:6px 14px 6px 5px!important",
               fixed = TRUE)
})

# The data each table sends to the browser, caught on its way to
# DT::datatable(): the rows are fetched by the browser from the server,
# which a test session does not do.
catch_tables <- function(env = parent.frame()) {
  caught <- new.env()
  real <- DT::datatable
  testthat::local_mocked_bindings(
    datatable = function(data, ...) {
      caught$tables <- c(caught$tables, list(data))
      real(data, ...)
    },
    .package = "DT", .env = env)
  caught
}

test_that("the enrichment view's table sends its p-values as text", {
  caught <- catch_tables()
  b <- omicsCore::new_analysis_bundle(
    "run_enrichment", input_info = list(omics_type = "proteomics"),
    params = list(type = "ora", database = "hallmark", organism = "Homo sapiens",
                  direction = "separate", comparison = "B_vs_A"),
    results = list(enrich_result_df = enrich_frame("ora")))
  shiny::testServer(enrich_view_server, args = list(diff_bundle = shiny::reactiveVal(NULL)), {
    enrich_bundle(b); is_demo(FALSE)
    session$setInputs(show_p = "adjusted", show_cutoff = 0.05)
    session$flushReact()
    invisible(output$hits)
  })
  d <- caught$tables[[length(caught$tables)]]
  expect_identical(d$`adjusted p`, c("1.80e-248", "4.42e-16", "0.0123"))
})

test_that("the Differential view's Top hits and the integration table do too", {
  caught <- catch_tables()
  set.seed(3)
  groups <- rep(c("Ctrl", "Treat"), each = 4)
  samp <- paste0("S", seq_along(groups))
  ids <- paste0("P", 1:40)
  m <- matrix(stats::rnorm(40 * 8, 12, 0.2), 40, dimnames = list(ids, samp))
  m[1:10, 5:8] <- m[1:10, 5:8] + 3
  inp <- omicsCore::omics_input(
    m, data.frame(group = groups, row.names = samp, stringsAsFactors = FALSE),
    data.frame(feature_id = ids, feature_symbol = paste0("G", 1:40), row.names = ids,
               stringsAsFactors = FALSE),
    omics_type = "proteomics", assay_type = "normalized_intensity")
  proj <- omicsCore::omics_project("p", list(proteomics = inp))
  shiny::testServer(diff_view_server, args = list(current_project = shiny::reactiveVal(proj)), {
    session$setInputs(group_col = "group", control = "Ctrl", case = "Treat",
                      method = "limma", fdr_cut = 0.05, fc_cut = 0.263, p_kind = "adj",
                      rerun = 1)
    session$elapse(300)
    invisible(output$hits)
    sig <- results$hits_df()
    d <- Filter(function(x) "Direction" %in% names(x), caught$tables)
    d <- d[[length(d)]]
    expect_identical(d$`adjusted p`, format_p_value(sig$adj_p_value))
    expect_identical(d$.sort_1, sig$adj_p_value)
    expect_type(d$`adjusted p`, "character")
  })

  tab <- integration_result_table(
    data.frame(feature_symbol = c("A", "B"), feature_id = c("A", "B"),
               effect = c(1, -1), adj_p_value = c(2e-300, 0.5), direction = c("up", "down"),
               stringsAsFactors = FALSE),
    "correlation", c("x", "y"))
  out <- dt_sortable_text(tab, p_cols = "Adj. p")$data
  expect_identical(out$`Adj. p`, c("2.00e-300", "0.5"))
})
