# The deployment image pins CRAN to a snapshot contemporary with
# Bioconductor 3.20, which carries msigdbr 7.5.1 -- the pre-v10 API. That
# branch of fetch_msigdbr_raw() is never exercised here, because CI and
# development both install a current msigdbr, so the one piece of it that
# can be wrong on its own is tested directly.

test_that("only KEGG is renamed for the pre-v10 msigdbr data", {
  # MSigDB 7.5 calls it CP:KEGG; the rename to CP:KEGG_LEGACY came with
  # KEGG_MEDICUS. DB_MSIGDBR_MAP carries the modern name.
  expect_identical(legacy_subcollection("CP:KEGG_LEGACY"), "CP:KEGG")

  # Everything else is spelled the same in both, and must pass through
  # untouched -- a translation that fired too widely would filter out
  # every gene set rather than error.
  for (sub in c("CP:REACTOME", "CP:WIKIPATHWAYS", "GO:BP", "GO:MF", "GO:CC")) {
    expect_identical(legacy_subcollection(sub), sub, info = sub)
  }
})

test_that("a collection with no subcollection survives the translation", {
  # hallmark carries NA, and the caller tests it with is.na() afterwards,
  # so the NA has to come back out as an NA of the same type.
  expect_true(is.na(legacy_subcollection(NA_character_)))
  expect_type(legacy_subcollection(NA_character_), "character")
})

test_that("every subcollection in the map translates to something usable", {
  # Guards the two from drifting: a new entry in DB_MSIGDBR_MAP whose
  # name differs on old msigdbr would otherwise be found in production.
  for (db in names(DB_MSIGDBR_MAP)) {
    sub <- DB_MSIGDBR_MAP[[db]]$subcollection
    out <- legacy_subcollection(sub)
    expect_length(out, 1L)
    if (!is.na(sub)) expect_true(nzchar(out), info = db)
  }
})
