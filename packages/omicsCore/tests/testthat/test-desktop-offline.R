test_that("offline cache misses and explicit refresh fail without network calls", {
  dir <- withr::local_tempdir()
  withr::local_envvar(OMICSCORE_OFFLINE = "1", OMICSCORE_GENESET_CACHE = dir)
  expect_error(fetch_msigdbr_table("hallmark", "Hs"), "Offline mode.*cache")
  expect_error(refresh_geneset_cache("kegg", force = TRUE), "Offline mode")
  expect_error(kegg_rest_table("/list/hsa"), "Offline mode")
  expect_length(list.files(dir), 0L)
})

test_that("switching directories and replacing tables invalidates session memoisation", {
  first <- withr::local_tempdir()
  second <- withr::local_tempdir()
  table <- function(label) data.frame(gs_name = label, gene_symbol = "TP53")
  file <- "hallmark__Homo_sapiens.qs2"
  qs2::qs_save(table("FIRST"), file.path(first, file))
  qs2::qs_save(table("SECOND"), file.path(second, file))
  withr::local_envvar(OMICSCORE_OFFLINE = "1", OMICSCORE_GENESET_CACHE = first)
  expect_equal(fetch_msigdbr_table("hallmark", "Hs")$gs_name, "FIRST")
  Sys.setenv(OMICSCORE_GENESET_CACHE = second)
  expect_equal(fetch_msigdbr_table("hallmark", "Hs")$gs_name, "SECOND")
  qs2::qs_save(table("REPLACED"), file.path(second, file))
  expect_equal(fetch_msigdbr_table("hallmark", "Hs")$gs_name, "REPLACED")
  unlink(file.path(second, file))
  expect_error(fetch_msigdbr_table("hallmark", "Hs"), "Offline mode")
})

test_that("offline KEGG caches are served without TTL refresh", {
  dir <- withr::local_tempdir()
  withr::local_envvar(OMICSCORE_OFFLINE = "1", OMICSCORE_GENESET_CACHE = dir,
                      OMICSCORE_GENESET_TTL_DAYS = "1")
  path <- file.path(dir, "kegg__Homo_sapiens.qs2")
  table <- data.frame(gs_name = "LOCAL", gene_symbol = "TP53",
                      gs_source = "KEGG REST local use")
  qs2::qs_save(table, path)
  Sys.setFileTime(path, Sys.time() - 90 * 86400)
  testthat::local_mocked_bindings(kegg_rest_table = function(...) stop("network touched"))
  expect_equal(fetch_msigdbr_table("kegg", "Hs"), table)
})

test_that("default HTML reports have no remote MathJax dependency", {
  skip_if_not_installed("rmarkdown")
  skip_if_not(rmarkdown::pandoc_available())
  mat <- matrix(1:12, nrow = 3, dimnames = list(paste0("P", 1:3), paste0("S", 1:4)))
  input <- omics_input(mat, data.frame(group = c("A", "A", "B", "B"),
                                      row.names = colnames(mat)),
                       data.frame(feature_id = rownames(mat)),
                       omics_type = "proteomics", assay_type = "normalized_intensity")
  project <- omics_project("offline report", list(protein = input))
  out <- withr::local_tempfile(fileext = ".html")
  export_report(project, out)
  html <- paste(readLines(out, warn = FALSE), collapse = "\n")
  expect_false(grepl("src=[\"']https?://[^\"']*[Mm]ath[Jj]ax", html))
  expect_gt(file.size(out), 1000)
})
