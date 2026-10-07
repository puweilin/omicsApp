# The Integration view's feature matching: the note saying how features
# meet, and an uploaded mapping table that is archived, saved with the
# project, used by the run and read again by the exported script.

fla_project <- function() {
  set.seed(5)
  groups <- rep(c("ctrl", "A"), each = 4)
  mk_meta <- function(samp) {
    data.frame(group = groups, donor = paste0("D", 1:8), row.names = samp,
               stringsAsFactors = FALSE)
  }
  genes <- c("TP53", "BRCA1", "MDM2", paste0("G", 1:9))
  prot_ids <- c("P04637", "P04637-2", "P38398", "Q00987", paste0("X", 1:9))
  sp <- paste0("p_", 1:8)
  sr <- paste0("r_", 1:8)
  prot <- matrix(stats::rnorm(13 * 8, 20, 0.3), 13, dimnames = list(prot_ids, sp))
  rna <- matrix(stats::rnorm(12 * 8, 8, 0.3), 12,
                dimnames = list(paste0("ENSG", 1:12), sr))
  prot[1:3, groups == "A"] <- prot[1:3, groups == "A"] + 2
  rna[1:2, groups == "A"] <- rna[1:2, groups == "A"] + 2
  omicsCore::omics_project("fla", list(
    proteomics = omicsCore::omics_input(
      prot, mk_meta(sp),
      data.frame(feature_id = prot_ids, feature_symbol = c("TP53", "TP53", "BRCA1", "MDM2",
                                                           paste0("G", 1:9)),
                 stringsAsFactors = FALSE),
      omics_type = "proteomics", assay_type = "normalized_intensity"),
    rnaseq = omicsCore::omics_input(
      rna, mk_meta(sr),
      data.frame(feature_id = rownames(rna), feature_symbol = genes,
                 stringsAsFactors = FALSE),
      omics_type = "rnaseq", assay_type = "logcpm")))
}

fla_diff <- function(proj) {
  omicsCore::run_diff(proj$experiments$proteomics, method = "limma",
                      group_col = "group", control_group = "ctrl", case_group = "A")
}

fla_html <- function(x) paste(as.character(x), collapse = "")

test_that("the note says how many features matched and that isoforms of a gene are each kept", {
  proj <- fla_project()
  shiny::testServer(integration_view_server, args = list(
    current_project = shiny::reactiveVal(proj),
    diff_bundle = shiny::reactiveVal(fla_diff(proj)),
    diff_layer = shiny::reactiveVal("proteomics")), {
    session$flushReact()
    html <- fla_html(output$feature_note)
    expect_match(html, "Matched by gene symbol.", fixed = TRUE)
    expect_match(html, "13 of 13 features of proteomics are matched with 12 of 12 features of rnaseq (13 pairs).",
                 fixed = TRUE)
    expect_match(html, "2 features of proteomics share a match in rnaseq", fixed = TRUE)
    df <- integration_bundle()$results$integration_df
    expect_equal(nrow(df), 13L)
    expect_setequal(df$feature_id_a[df$feature_symbol == "TP53"], c("P04637", "P04637-2"))
    tab <- integration_result_table(df, "concordance", c("proteomics", "rnaseq"))
    expect_true(all(c("TP53 (P04637)", "TP53 (P04637-2)") %in% tab$Feature))
    expect_true("BRCA1" %in% tab$Feature)
  })
})

test_that("an uploaded mapping table is archived, saved with the project, and used by the run", {
  store <- withr::local_tempdir()
  withr::local_envvar(OMICSAPP_DATA_DIR = store)
  proj <- fla_project()
  # Gene symbols on the protein side that match nothing: only the table
  # can pair the layers.
  proj$experiments$proteomics$feature_df$feature_symbol <- "?"
  map <- file.path(withr::local_tempdir(), "uniprot_map.csv")
  utils::write.csv(data.frame(Description = c("p53", "BRCA1"),
                              Gene = c("TP53", "BRCA1"),
                              Accession = c("P04637", "P38398")),
                   map, row.names = FALSE)
  cp <- shiny::reactiveVal(proj)
  shiny::testServer(integration_view_server, args = list(
    current_project = cp,
    diff_bundle = shiny::reactiveVal(fla_diff(proj)),
    diff_layer = shiny::reactiveVal("proteomics")), {
    session$flushReact()
    expect_match(fla_html(output$feature_note), "No features matched", fixed = TRUE)

    session$setInputs(link_file = data.frame(name = "uniprot_map.csv", size = file.size(map),
                                             type = "text/csv", datapath = map))
    # The columns are guessed from which identifiers each layer has.
    up <- link_upload()
    expect_identical(up$guess, c("Accession", "Gene"))
    session$setInputs(link_col_a = "Accession", link_col_b = "Gene")
    expect_match(fla_html(output$link_preview),
                 "With this table: 3 of 13 features of proteomics are matched with 2 of 12",
                 fixed = TRUE)
    session$setInputs(use_link = 1)
    session$flushReact()

    fl <- cp()$feature_link
    expect_identical(names(fl), c("proteomics", "rnaseq"))
    src <- attr(fl, "source")
    expect_true(startsWith(normalizePath(src$path), normalizePath(file.path(store, "raw"))))
    expect_true(file.exists(src$path))
    expect_identical(src$columns, c(proteomics = "Accession", rnaseq = "Gene"))
    expect_match(fla_html(output$feature_note),
                 "Matched by the mapping table saved with this project.", fixed = TRUE)

    b <- integration_bundle()
    expect_identical(b$params$feature_link_source$source, "project")
    expect_identical(b$params$feature_link_source$path, src$path)
    expect_setequal(b$results$integration_df$feature_id_a, c("P04637", "P04637-2", "P38398"))

    # The exported script reads the archived table.
    p <- cp()
    p$bundles <- list(diff = fla_diff(p), integration = b)
    txt <- paste(omicsCore::export_script(p, include_plots = FALSE), collapse = "\n")
    expect_match(txt, "feature_link <- read_feature_link(", fixed = TRUE)
    expect_match(txt, sprintf('"raw/%s"', basename(src$path)), fixed = TRUE)

    # And can be put aside again.
    session$setInputs(drop_link = 1)
    session$flushReact()
    expect_null(cp()$feature_link)
    expect_match(fla_html(output$feature_note), "Matched by gene symbol.", fixed = TRUE)
  })
})

test_that("removing a layer removes a mapping table that names it", {
  proj <- fla_project()
  proj$feature_link <- data.frame(proteomics = "P04637", rnaseq = "TP53")
  cp <- shiny::reactiveVal(proj)
  shiny::testServer(project_view_server, args = list(current_project = cp), {
    session$flushReact()
    session$setInputs(drop_layer_2 = 1)
    session$flushReact()
    expect_identical(pending_drop(), "rnaseq")
    session$setInputs(confirm_drop_layer = 1)
    session$flushReact()
    expect_identical(names(cp()$experiments), "proteomics")
    expect_null(cp()$feature_link)
  })
})

test_that("the script download carries the mapping table it reads", {
  withr::local_envvar(OMICSAPP_DATA_DIR = withr::local_tempdir())
  p <- fla_project()
  map <- file.path(withr::local_tempdir(), "map__0123456789ab.csv")
  utils::write.csv(data.frame(Accession = c("P04637", "P38398"), Gene = c("TP53", "BRCA1")),
                   map, row.names = FALSE)
  p$feature_link <- omicsCore::read_feature_link(map, c(proteomics = "Accession",
                                                         rnaseq = "Gene"))
  d <- fla_diff(p)
  partner <- omicsCore::run_diff(p$experiments$rnaseq, method = "limma", group_col = "group",
                                 control_group = "ctrl", case_group = "A")
  p$bundles <- list(diff = d, integration = omicsCore::run_integration(
    p, "concordance", c("proteomics", "rnaseq"),
    diff_bundles = list(proteomics = d, rnaseq = partner)))
  shiny::testServer(report_view_server, args = list(current_project = shiny::reactiveVal(p)), {
    f <- output$download_bundle
    listing <- if (grepl("zip$", f)) utils::unzip(f, list = TRUE)$Name else utils::untar(f, list = TRUE)
    expect_true(file.path("analysis", "raw", basename(map)) %in% listing)
  })
})
