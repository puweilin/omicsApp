# pathway_name was taken from msigdbr's gs_description. For Hallmark
# that is a short sentence, which is why it went unnoticed; for GO it is
# the full definition with its citations --
#   "The chemical reactions and pathways involving
#    10-formyltetrahydrofolate, the formylated derivative of
#    tetrahydrofolate. [GOC:ai]"
# -- and that was arriving as the label on a dotplot axis.

test_that("the collection prefix is stripped", {
  expect_identical(prettify_gene_set_name("HALLMARK_TNFA_SIGNALING_VIA_NFKB"),
                   "TNFA SIGNALING VIA NFKB")
  expect_identical(prettify_gene_set_name("GOBP_KERATINIZATION"),
                   "KERATINIZATION")
  expect_identical(prettify_gene_set_name("REACTOME_SIGNALING_BY_INTERLEUKINS"),
                   "SIGNALING BY INTERLEUKINS")
  expect_identical(prettify_gene_set_name("WP_FERROPTOSIS"), "FERROPTOSIS")
})

test_that("the longer KEGG prefixes win over the shorter one", {
  # KEGG_ is a prefix of KEGG_MEDICUS_, so order matters.
  expect_identical(prettify_gene_set_name("KEGG_MEDICUS_REFERENCE_APOPTOSIS"),
                   "REFERENCE APOPTOSIS")
  expect_identical(prettify_gene_set_name("KEGG_RIBOSOME"), "RIBOSOME")
})

test_that("casing is left as the source wrote it", {
  # Lowercasing reads better on a long GO term, but every rule for
  # deciding which tokens to spare turns some gene symbol into a word,
  # and "Tnfa" is wrong where ALL CAPS is only ugly.
  out <- prettify_gene_set_name("HALLMARK_IL6_JAK_STAT3_SIGNALING")
  expect_identical(out, "IL6 JAK STAT3 SIGNALING")
  expect_true(grepl("STAT3", out, fixed = TRUE))
})

test_that("a name with no known prefix is left alone but still cleaned", {
  expect_identical(prettify_gene_set_name("CUSTOM_SET_ONE"), "CUSTOM SET ONE")
})

test_that("NA and empty survive", {
  expect_true(is.na(prettify_gene_set_name(NA_character_)))
  expect_identical(prettify_gene_set_name(""), "")
  expect_length(prettify_gene_set_name(character(0)), 0L)
})

test_that("term2name carries the name, not the definition", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("msigdbr")
  # This is the one test here that fetches for real, and the fetch is
  # memoised for the session. Left behind, the cached table is what the
  # geneset-cache tests then read instead of their own fixture.
  on.exit(cache_drop("msig::go_bp::Homo sapiens"), add = TRUE)
  t <- build_term_tables("go_bp", "Hs")
  # A GO definition is a sentence and ends in a citation bracket; a name
  # is neither.
  expect_false(any(grepl("^The ", t$term2name$name)))
  expect_false(any(grepl("\\[GOC:", t$term2name$name)))
  # And it is short enough to put on an axis.
  expect_lt(stats::median(nchar(t$term2name$name)), 60)
})
