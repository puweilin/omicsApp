# The import matches the gene-symbol column by heading, and a heading
# can be absent or wrong. Neither shows up in the outcome: a wrong match
# relabels every feature with a plausible-looking name, and no match
# leaves feature_symbol as the accession, after which enrichment returns
# an empty result rather than an error -- clusterProfiler answers "No
# gene can be mapped" and hands back NULL. So the import says which
# column it took, before anything is run.

parsed_with <- function(symbol_col = "Gene names") {
  a <- data.frame(Accession = c("A0A024RBG1", "A0A075B6I0", "A0A075B6I9"),
                  Coverage = c(1, 2, 3), stringsAsFactors = FALSE)
  if (!is.null(symbol_col)) {
    a[[symbol_col]] <- c("NUDT4B", "IGLV4-69", "IGLV8-61")
  }
  feat <- omicsCore:::materialize_feature_annot(a, a$Accession)
  m <- matrix(rnorm(3 * 4, 20, 1), 3, 4,
              dimnames = list(a$Accession, sprintf("S%d", 1:4)))
  meta <- data.frame(condition = rep(c("G1", "G2"), each = 2),
                     row.names = colnames(m))
  list(input = omicsCore::omics_input(
         m, meta, feat, omics_type = "proteomics",
         assay_type = "normalized_intensity"),
       report = omicsCore::new_import_report(source = "x.xlsx"))
}

render_source <- function(p) {
  out <- NULL
  shiny::testServer(import_view_server, args = list(), {
    parsed(p)
    session$flushReact()
    out <<- paste(unlist(output$confirm_symbol_source), collapse = " ")
  })
  out
}

test_that("the column the symbol came from is named", {
  html <- render_source(parsed_with("Gene names"))
  expect_true(grepl("Gene symbol", html, fixed = TRUE))
  expect_true(grepl("Gene names", html, fixed = TRUE))
})

test_that("a few real symbols are shown, so a wrong column is visible", {
  # Naming the column is not enough on its own: the check a reader can
  # actually perform is whether the values look like genes.
  html <- render_source(parsed_with("Gene names"))
  expect_true(grepl("NUDT4B", html, fixed = TRUE))
})

test_that("the heading is recognised whatever the vendor called it", {
  for (nm in c("Genes", "PG.Genes", "gene_symbol")) {
    html <- render_source(parsed_with(nm))
    expect_true(grepl(nm, html, fixed = TRUE), info = nm)
    expect_true(grepl("NUDT4B", html, fixed = TRUE), info = nm)
  }
})

test_that("no gene column says so, and says what it costs", {
  # This is the case worth catching: everything downstream still runs,
  # and enrichment comes back empty with nothing to explain it.
  html <- render_source(parsed_with(NULL))
  expect_true(grepl("no gene column found", html, fixed = TRUE))
  expect_true(grepl("Enrichment", html, fixed = TRUE))
  expect_false(grepl("NUDT4B", html, fixed = TRUE))
})

test_that("nothing is claimed before a file is parsed", {
  shiny::testServer(import_view_server, args = list(), {
    expect_null(output$confirm_symbol_source)
  })
})
