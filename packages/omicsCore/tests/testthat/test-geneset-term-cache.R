# The built TERM2GENE / TERM2NAME tables are remembered between
# enrichment runs. What has to hold: a remembered table is the table a
# fresh build would give, a different database or organism is a
# different entry, a refreshed gene-set file is never served from the
# old entry, and the cache stays a few entries long.
#
# Most tests here read gene sets from a temporary on-disk cache, so they
# run offline and without msigdbr's data; one checks the msigdbr path.

local_geneset_dir <- function(files = list(), env = parent.frame()) {
  dir <- withr::local_tempdir(.local_envir = env)
  withr::local_envvar(OMICSCORE_GENESET_CACHE = dir, .local_envir = env)
  for (nm in names(files)) qs2::qs_save(files[[nm]], file.path(dir, nm))
  clear_term_table_cache()
  withr::defer(clear_term_table_cache(), envir = env)
  dir
}

sets_table <- function(prefix, n_sets = 3L, per_set = 12L) {
  data.frame(
    gs_name = rep(sprintf("%s_SET_%d", prefix, seq_len(n_sets)), each = per_set),
    gene_symbol = sprintf("G%03d", seq_len(n_sets * per_set)),
    gs_description = "d",
    gs_id = rep(sprintf("M%d", seq_len(n_sets)), each = per_set),
    stringsAsFactors = FALSE
  )
}

# Counts the trips to the raw table, which is what the cache saves.
local_fetch_counter <- function(env = parent.frame()) {
  real <- fetch_msigdbr_table
  calls <- new.env()
  calls$n <- 0L
  testthat::local_mocked_bindings(
    fetch_msigdbr_table = function(...) {
      calls$n <- calls$n + 1L
      real(...)
    },
    .env = env
  )
  calls
}

cached_keys <- function() ls(.term_table_cache)

test_that("a second build is the remembered one, identical to a fresh build", {
  skip_if_not_installed("qs2")
  local_geneset_dir(list("go_bp__Homo_sapiens.qs2" = sets_table("GOBP")))
  calls <- local_fetch_counter()

  first <- build_term_tables("go_bp", "Hs")
  second <- build_term_tables("go_bp", "Homo sapiens")
  expect_identical(second, first)
  expect_identical(first, build_term_tables_uncached("go_bp", "Homo sapiens"))
  # Two builds, one through the cache: the uncached call above is the
  # second fetch.
  expect_identical(calls$n, 2L)
  expect_length(cached_keys(), 1L)
})

test_that("organism and database are each part of the key", {
  skip_if_not_installed("qs2")
  local_geneset_dir(list(
    "go_bp__Homo_sapiens.qs2" = sets_table("HS"),
    "go_bp__Mus_musculus.qs2" = sets_table("MM"),
    "hallmark__Homo_sapiens.qs2" = sets_table("HALLMARK")
  ))
  hs <- build_term_tables("go_bp", "Hs")
  mm <- build_term_tables("go_bp", "Mm")
  hm <- build_term_tables("hallmark", "Hs")
  expect_length(cached_keys(), 3L)
  expect_true(all(startsWith(hs$term2gene$term, "HS_")))
  expect_true(all(startsWith(mm$term2gene$term, "MM_")))
  expect_true(all(startsWith(hm$term2gene$term, "HALLMARK_")))
})

test_that("a gene-set file replaced behind the session's back is read again", {
  skip_if_not_installed("qs2")
  # The monthly cron, or refresh_geneset_cache() in another R process:
  # this session's memo is not told, so the file itself has to be.
  dir <- local_geneset_dir(list("go_bp__Homo_sapiens.qs2" = sets_table("OLD")))
  path <- file.path(dir, "go_bp__Homo_sapiens.qs2")
  old <- build_term_tables("go_bp", "Hs")
  expect_true(all(startsWith(old$term2gene$term, "OLD_")))

  qs2::qs_save(sets_table("NEW", n_sets = 4L), path)
  Sys.setFileTime(path, Sys.time() + 60)
  new <- build_term_tables("go_bp", "Hs")
  expect_true(all(startsWith(new$term2gene$term, "NEW_")))
  expect_identical(new, build_term_tables_uncached("go_bp", "Homo sapiens"))
  # The entry built from the old file went with the old memo.
  expect_length(cached_keys(), 1L)
})

test_that("refresh_geneset_cache() invalidates the built tables", {
  skip_if_not_installed("qs2")
  local_geneset_dir()
  withr::defer(cache_drop("msig::kegg::Homo sapiens"))
  rest <- function(extra) {
    genes <- c("hsa:1", "hsa:2", if (extra) "hsa:5")
    list(
      "/link/pathway/hsa" = data.frame(V1 = genes, V2 = "path:hsa00010"),
      "/list/pathway/hsa" = data.frame(
        V1 = "hsa00010", V2 = "Glycolysis / Gluconeogenesis - Homo sapiens (human)"),
      "/list/hsa" = data.frame(V1 = genes, V2 = "CDS", V3 = "1:1..2",
                               V4 = paste0(c("GENEA", "GENEB", if (extra) "GENEE"), "; g"))
    )
  }
  serve <- function(tables) function(path) tables[[path]]

  testthat::local_mocked_bindings(kegg_rest_table = serve(rest(FALSE)))
  refresh_geneset_cache("kegg", "Hs", quiet = TRUE)
  before <- build_term_tables("kegg", "Hs")
  expect_false("GENEE" %in% before$term2gene$gene)

  testthat::local_mocked_bindings(kegg_rest_table = serve(rest(TRUE)))
  refresh_geneset_cache("kegg", "Hs", force = TRUE, quiet = TRUE)
  after <- build_term_tables("kegg", "Hs")
  expect_true("GENEE" %in% after$term2gene$gene)
})

test_that("the cache holds four entries and lets the least recently used go", {
  skip_if_not_installed("qs2")
  dbs <- c("hallmark", "reactome", "wikipathways", "go_bp", "go_mf")
  files <- stats::setNames(lapply(toupper(dbs), sets_table),
                           sprintf("%s__Homo_sapiens.qs2", dbs))
  local_geneset_dir(files)
  for (db in dbs[1:4]) build_term_tables(db, "Hs")
  expect_length(cached_keys(), 4L)
  # Touching hallmark makes reactome the oldest.
  build_term_tables("hallmark", "Hs")
  build_term_tables("go_mf", "Hs")
  keys <- cached_keys()
  expect_length(keys, TERM_TABLE_CACHE_SIZE)
  expect_false(any(startsWith(keys, "terms::reactome::")))
  expect_true(any(startsWith(keys, "terms::hallmark::")))
})

test_that("from msigdbr, a repeat build does not go back to msigdbr", {
  skip_if_not_installed("msigdbr")
  skip_if_not_installed("qs2")
  # An empty cache directory, so the tables come from msigdbr itself.
  local_geneset_dir()
  probe <- tryCatch(fetch_msigdbr_raw("hallmark", "Homo sapiens"),
                    error = function(e) NULL)
  skip_if(is.null(probe) || nrow(probe) == 0L, "msigdbr data unavailable")
  real_raw <- fetch_msigdbr_raw
  raw_calls <- 0L
  testthat::local_mocked_bindings(fetch_msigdbr_raw = function(...) {
    raw_calls <<- raw_calls + 1L
    real_raw(...)
  })
  cache_drop("msig::hallmark::Homo sapiens")

  first <- build_term_tables("hallmark", "Hs")
  second <- build_term_tables("hallmark", "Hs")
  expect_identical(second, first)
  expect_identical(raw_calls, 1L)
  expect_true(any(grepl("::msigdbr:", cached_keys(), fixed = TRUE)))
  expect_gt(nrow(first$term2gene), 1000L)
})
