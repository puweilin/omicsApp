# Feature links: which protein is which gene. Every feature is kept (two
# isoforms of one gene are two pairs), a stated link replaces symbol
# matching, and each method says what it does with several proteins of
# one gene.

fl_project <- function(prot_symbols = NULL) {
  set.seed(42)
  samp_p <- paste0("PS", 1:8)
  samp_r <- paste0("RS", 1:8)
  groups <- rep(c("ctrl", "case"), each = 4)
  donors <- paste0("D", 1:8)
  genes <- c("TP53", "BRCA1", "MDM2", "EGFR", paste0("G", 1:7))
  prot_ids <- c("P04637", "P04637-2", "P38398", "Q00987", "P00533",
                paste0("X", 1:7))
  prot_sym <- prot_symbols %||% c("TP53", "TP53", "BRCA1", "MDM2", "EGFR",
                                  paste0("G", 1:7))
  rna <- matrix(rnorm(length(genes) * 8, 8, 1), length(genes),
                dimnames = list(paste0("ENSG", seq_along(genes)), samp_r))
  rna[1:3, groups == "case"] <- rna[1:3, groups == "case"] + 3
  prot <- matrix(rnorm(length(prot_ids) * 8, 20, 0.2), length(prot_ids),
                 dimnames = list(prot_ids, samp_p))
  # The canonical TP53 follows its transcript; the isoform goes the other
  # way, so the two pairs give different answers.
  prot[1, ] <- 20 + rna[1, ] / 2 + rnorm(8, 0, 0.05)
  prot[2, ] <- 30 - rna[1, ] / 2 + rnorm(8, 0, 0.05)
  prot[3, groups == "case"] <- prot[3, groups == "case"] + 2
  mk_meta <- function(samp) {
    data.frame(group = groups, donor = donors, row.names = samp,
               stringsAsFactors = FALSE)
  }
  omics_project("fl", list(
    prot = omics_input(prot, mk_meta(samp_p),
                       data.frame(feature_id = prot_ids, feature_symbol = prot_sym,
                                  stringsAsFactors = FALSE),
                       omics_type = "proteomics", assay_type = "normalized_intensity"),
    rna = omics_input(rna, mk_meta(samp_r),
                      data.frame(feature_id = rownames(rna), feature_symbol = genes,
                                 stringsAsFactors = FALSE),
                      omics_type = "rnaseq", assay_type = "logcpm")))
}

fl_diffs <- function(p) {
  lapply(p$experiments, function(e) {
    run_diff(e, method = "ttest", group_col = "group",
             control_group = "ctrl", case_group = "case")
  })
}

# ---- pairing ------------------------------------------------------------

test_that("isoforms of one gene are both kept when matching symbols", {
  p <- fl_project()
  prev <- feature_pairing_preview(p, "prot", "rna")
  expect_identical(prev$source, "symbol")
  expect_equal(prev$n_pairs, 12L)
  expect_identical(prev$n_a_sharing, 2L)
  expect_identical(prev$n_b_shared, 1L)
  tp53 <- prev$pairs[prev$pairs$feature_b == "ENSG1", ]
  expect_setequal(tp53$feature_a, c("P04637", "P04637-2"))
  expect_setequal(tp53$feature, c("TP53 (P04637)", "TP53 (P04637-2)"))
  # A gene with one pair keeps the plain name, as results always had.
  expect_identical(prev$pairs$feature[prev$pairs$feature_a == "P38398"], "BRCA1")
})

test_that("a link given at accession level covers the isoforms, and replaces symbol matching", {
  # The protein layer's symbols are wrong for everything: only the link
  # can pair it.
  p <- fl_project(prot_symbols = rep("WRONG", 12))
  expect_error(run_integration(p, "correlation", c("prot", "rna")), "No shared")
  link <- data.frame(prot = c("P04637", "P38398"), rna = c("TP53", "BRCA1"))
  prev <- feature_pairing_preview(p, "prot", "rna", feature_link = link)
  expect_identical(prev$source, "supplied")
  expect_setequal(paste(prev$pairs$feature_a, prev$pairs$feature_b),
                  c("P04637 ENSG1", "P04637-2 ENSG1", "P38398 ENSG2"))

  # Symbols that agree are not used once a link is given: MDM2 and EGFR
  # match by name but the link leaves them out.
  p2 <- fl_project()
  prev2 <- feature_pairing_preview(p2, "prot", "rna", feature_link = link)
  expect_false(any(prev2$pairs$feature_a %in% c("Q00987", "P00533")))
  expect_equal(prev2$n_paired_a, 3L)

  # An entry for the isoform itself wins over its canonical accession.
  link3 <- rbind(link, data.frame(prot = "P04637-2", rna = "MDM2"))
  prev3 <- feature_pairing_preview(p2, "prot", "rna", feature_link = link3)
  expect_identical(prev3$pairs$feature_b[prev3$pairs$feature_a == "P04637-2"], "ENSG3")
  expect_identical(prev3$pairs$feature_b[prev3$pairs$feature_a == "P04637"], "ENSG1")

  # Protein-group ids and FASTA-style ids are matched by their members.
  expect_equal(uniprot_canonical(c("P04637-2", "NKX2-1", "Q9Y6K9-12", "TP53")),
               c("P04637", NA, "Q9Y6K9", NA))
  keys <- feature_link_keys(c("P04637-2;Q00987", "sp|P38398|BRCA1_HUMAN"), c(NA, NA))
  expect_true(all(c("Q00987", "P04637", "P38398") %in% keys$key))
})

test_that("the project's link is used, and one passed to the run overrides it", {
  p <- fl_project(prot_symbols = rep("WRONG", 12))
  d <- fl_diffs(p)
  p$feature_link <- data.frame(prot = c("P04637", "P38398", "Q00987"),
                               rna = c("TP53", "BRCA1", "MDM2"))
  b <- run_integration(p, "concordance", c("prot", "rna"), diff_bundles = d)
  expect_identical(b$params$feature_link_source$source, "project")
  expect_equal(nrow(b$results$integration_df), 4L)
  # The gene is the one the link names, not the protein table's own
  # (here wrong) symbol column.
  expect_setequal(b$results$integration_df$feature_symbol, c("TP53", "BRCA1", "MDM2"))
  expect_true("TP53 (P04637-2)" %in% b$results$integration_df$feature_id)
  only <- data.frame(prot = "P38398", rna = "ENSG2")
  b2 <- run_integration(p, "concordance", c("prot", "rna"), diff_bundles = d,
                        feature_link = only)
  expect_identical(b2$params$feature_link_source$source, "supplied")
  expect_identical(b2$results$integration_df$feature_id_a, "P38398")
  # Without a link, the run says it matched symbols.
  p0 <- fl_project()
  b0 <- run_integration(p0, "concordance", c("prot", "rna"), diff_bundles = fl_diffs(p0))
  expect_identical(b0$params$feature_link_source, list(source = "symbol"))

  expect_error(run_integration(p, "concordance", c("prot", "rna"), diff_bundles = d,
                               feature_link = data.frame(x = "P38398", rna = "TP53")),
               "missing: prot")
  expect_error(omics_project("x", p$experiments,
                             feature_link = data.frame(a = 1, b = 2)),
               "named after the layers")
})

# ---- each method with several proteins per gene ---------------------------

test_that("concordance reports each protein-gene pair as a row of its own", {
  p <- fl_project()
  d <- fl_diffs(p)
  b <- run_integration(p, "concordance", c("prot", "rna"), diff_bundles = d)
  df <- b$results$integration_df
  expect_equal(nrow(df), 12L)
  expect_false(anyDuplicated(df$feature_id) > 0L)
  tp53 <- df[df$feature_symbol == "TP53", ]
  ra <- d$prot$results$diff_result_df
  expect_equal(tp53$effect_a[tp53$feature_id_a == "P04637-2"],
               ra$effect[ra$feature_id == "P04637-2"])
  # The isoform moves against its transcript, the canonical with it.
  expect_setequal(tp53$quadrant, c("up_up", "down_up"))
  expect_equal(b$params$method_info$feature_pairing$n_pairs, 12L)
  # Corrected across pairs.
  expect_equal(df$adj_p_value, stats::p.adjust(df$p_value, "BH"))
})

test_that("correlation correlates each pair on its own data", {
  p <- fl_project()
  b <- run_integration(p, "correlation", c("prot", "rna"), cor_method = "pearson")
  df <- b$results$integration_df
  expect_equal(nrow(df), 12L)
  canon <- df[df$feature_id_a == "P04637", ]
  iso <- df[df$feature_id_a == "P04637-2", ]
  expect_gt(canon$effect, 0.9)
  expect_lt(iso$effect, -0.9)
  expect_equal(iso$effect,
               stats::cor(p$experiments$prot$expr_mat["P04637-2", ],
                          p$experiments$rna$expr_mat["ENSG1", ]))
  expect_identical(iso$feature_id, "TP53 (P04637-2)")
  expect_identical(iso$feature_symbol, "TP53")
})

test_that("ActivePathways scores a gene by its best feature, corrected for how many it has", {
  skip_if_not_installed("ActivePathways")
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("msigdbr")
  p <- fl_project()
  d <- fl_diffs(p)
  captured <- new.env()
  local_mocked_bindings(
    get_gene_set_list = function(...) {
      list(SET1 = c("TP53", "BRCA1", "MDM2", "EGFR", "G1"), SET2 = paste0("G", 1:7))
    })
  local_mocked_bindings(
    ActivePathways = function(scores, gmt, background, geneset_filter, significant,
                              merge_method, cytoscape_file_tag,
                              scores_direction = NULL, constraints_vector = NULL) {
      captured$scores <- scores
      captured$direction <- scores_direction
      NULL
    },
    .package = "ActivePathways")
  b <- run_integration(p, "active_pathways", c("prot", "rna"), diff_bundles = d,
                       geneset_filter = c(1L, 100L))
  ra <- d$prot$results$diff_result_df
  pv <- stats::setNames(ra$p_value, ra$feature_id)
  ef <- stats::setNames(ra$effect, ra$feature_id)
  best <- names(which.min(pv[c("P04637", "P04637-2")]))
  expect_equal(unname(captured$scores["TP53", "prot"]), 1 - (1 - pv[[best]])^2)
  expect_equal(unname(captured$direction["TP53", "prot"]), sign(ef[[best]]))
  # One feature: its own p-value, as before.
  expect_equal(unname(captured$scores["BRCA1", "prot"]), pv[["P38398"]])
  expect_equal(sum(rownames(captured$scores) == "TP53"), 1L)
  expect_true(any(grepl("1 gene(s) are measured by more than one feature in 'prot'",
                        b$warnings, fixed = TRUE)))
  expect_identical(b$params$method_info$n_genes_multi, c(1L, 0L))
})

# ---- reading a link -------------------------------------------------------

test_that("read_feature_link reads the chosen columns and remembers the file", {
  f <- tempfile(fileext = ".csv")
  on.exit(unlink(f), add = TRUE)
  utils::write.csv(data.frame(Note = c("a", "b", "c", "d"),
                              UniProt = c("P04637", " P38398 ", "", "P04637"),
                              Gene = c("TP53", "BRCA1", "MDM2", "TP53")),
                   f, row.names = FALSE)
  link <- read_feature_link(f, c(prot = "UniProt", rna = "Gene"))
  expect_identical(names(link), c("prot", "rna"))
  expect_identical(link$prot, c("P04637", "P38398"))
  expect_identical(attr(link, "source")$path, f)
  expect_identical(attr(link, "source")$columns, c(prot = "UniProt", rna = "Gene"))
  by_pos <- read_feature_link(f, c(prot = 2, rna = 3))
  expect_identical(by_pos$rna, link$rna)
  expect_error(read_feature_link(f, c(prot = "Accession", rna = "Gene")),
               "no column 'Accession'")
  expect_error(read_feature_link(f, c("UniProt", "Gene")), "name two layers")
  # Without columns: the whole file as text, to choose the columns from.
  whole <- read_feature_link(f, columns = NULL)
  expect_identical(names(whole), c("Note", "UniProt", "Gene"))
  expect_identical(whole$UniProt[[2L]], "P38398")

  skip_if_not_installed("openxlsx")
  x <- tempfile(fileext = ".xlsx")
  on.exit(unlink(x), add = TRUE)
  openxlsx::write.xlsx(data.frame(UniProt = c("P04637", "P38398"),
                                  Gene = c("TP53", "BRCA1")), x)
  expect_identical(read_feature_link(x, c(prot = "UniProt", rna = "Gene"))$rna,
                   c("TP53", "BRCA1"))
})

# ---- the exported script --------------------------------------------------

# The integration section of a script, run in-process against the inputs
# it names, from the folder the archive would be unpacked in.
run_integration_section <- function(lines, work, inputs) {
  start <- grep("^# ---- Multi-omics integration", lines)
  end <- grep("^integration <- run_integration\\(", lines)
  end <- end + which(lines[end:length(lines)] == ")")[[1L]] - 1L
  env <- new.env(parent = asNamespace("omicsCore"))
  for (nm in names(inputs)) assign(nm, inputs[[nm]], envir = env)
  old <- setwd(work)
  on.exit(setwd(old), add = TRUE)
  eval(parse(text = lines[start:end]), envir = env)
  env$integration
}

test_that("the exported script reads the archived link and reproduces the pairs", {
  work <- tempfile("fl-script")
  dir.create(file.path(work, "raw"), recursive = TRUE)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  map <- file.path(work, "raw", "protein_map__abc123.tsv")
  utils::write.table(data.frame(Accession = c("P04637", "P38398", "Q00987"),
                                Symbol = c("TP53", "BRCA1", "MDM2")),
                     map, sep = "\t", row.names = FALSE, quote = FALSE)

  p <- fl_project(prot_symbols = rep("WRONG", 12))
  p$feature_link <- read_feature_link(map, c(prot = "Accession", rna = "Symbol"))
  d <- fl_diffs(p)
  b <- run_integration(p, "concordance", c("prot", "rna"), diff_bundles = d)
  expect_identical(b$params$feature_link_source$path, map)
  p$bundles <- list(diff = d$prot, integration = b)

  lines <- export_script(p, include_plots = FALSE)
  txt <- paste(lines, collapse = "\n")
  expect_match(txt, 'feature_link <- read_feature_link(', fixed = TRUE)
  expect_match(txt, '"raw/protein_map__abc123.tsv"', fixed = TRUE)
  expect_match(txt, 'columns = c(prot = "Accession", rna = "Symbol")', fixed = TRUE)
  expect_match(txt, "feature_link\\s+= feature_link")
  expect_false(grepl("feature_link_source", txt, fixed = TRUE))
  expect_true(!inherits(tryCatch(parse(text = lines), error = function(e) e), "error"))

  again <- run_integration_section(lines, work, list(input_prot = p$experiments$prot,
                                                     input_rna = p$experiments$rna))
  expect_equal(again$results$integration_df[c("feature_id", "feature_id_a", "feature_id_b",
                                              "p_value")],
               b$results$integration_df[c("feature_id", "feature_id_a", "feature_id_b",
                                          "p_value")])
})

test_that("a small link from no file is written out; a large one is asked for", {
  p <- fl_project(prot_symbols = rep("WRONG", 12))
  d <- fl_diffs(p)
  link <- data.frame(prot = c("P04637", "P38398"), rna = c("TP53", "BRCA1"))
  b <- run_integration(p, "concordance", c("prot", "rna"), diff_bundles = d,
                       feature_link = link)
  p$bundles <- list(diff = d$prot, integration = b)
  lines <- export_script(p, include_plots = FALSE)
  expect_true(any(grepl('^feature_link <- data.frame\\(prot = c\\("P04637", "P38398"\\)', lines)))
  again <- run_integration_section(lines, tempdir(),
                                   list(input_prot = p$experiments$prot,
                                        input_rna = p$experiments$rna))
  expect_equal(again$results$integration_df$p_value, b$results$integration_df$p_value)

  # Too long to write out and no file behind it: a placeholder that
  # fails where it stands, and a note, rather than symbol matching.
  b$params$feature_link_source$table <- NULL
  b$params$feature_link_source$n_rows <- 5000L
  p$bundles$integration <- b
  lines <- export_script(p, include_plots = FALSE)
  expect_true(any(grepl("<path-to-feature-link-file>", lines, fixed = TRUE)))
  expect_true(any(grepl("NOTE: run_integration() matched features with a feature link (5000 rows)",
                        lines, fixed = TRUE)))

  # A result saved before feature links: nothing about them in the script.
  b0 <- run_integration(fl_project(), "concordance", c("prot", "rna"),
                        diff_bundles = fl_diffs(fl_project()))
  b0$params$feature_link_source <- NULL
  p$bundles$integration <- b0
  expect_false(any(grepl("feature_link", export_script(p, include_plots = FALSE))))
})
