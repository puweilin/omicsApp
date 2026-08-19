# The live KEGG refresh exists so `kegg` can mean current pathways instead
# of the 2011 MSigDB snapshot. Every test here runs against a mocked
# network seam (kegg_rest_table); nothing touches rest.kegg.jp. The two
# properties that matter most: a failure never damages an existing cache,
# and no code path fetches unless a live-sourced cache invited it.

local_cache_env <- function(env = parent.frame()) {
  dir <- file.path(tempdir(), paste0("gsr-", as.integer(runif(1, 1, 1e9))))
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  withr::local_envvar(OMICSCORE_GENESET_CACHE = dir, .local_envir = env)
  withr::defer(unlink(dir, recursive = TRUE), envir = env)
  dir
}

fake_rest <- function(extra_gene = FALSE) {
  link <- data.frame(
    V1 = c("hsa:1", "hsa:2", "hsa:2", "hsa:3", "hsa:4"),
    V2 = c("path:hsa00010", "path:hsa00010", "path:hsa00020",
           "path:hsa00020", "path:hsa00020"),
    stringsAsFactors = FALSE
  )
  genes <- data.frame(
    V1 = c("hsa:1", "hsa:2", "hsa:3", "hsa:4"),
    V2 = "CDS",
    V3 = "1:100..200",
    V4 = c("GENEA, ALIASA; gene a", "GENEB; gene b",
           "uncharacterized protein", "GENED; gene d"),
    stringsAsFactors = FALSE
  )
  if (extra_gene) {
    link <- rbind(link, data.frame(V1 = "hsa:5", V2 = "path:hsa00010"))
    genes <- rbind(genes, data.frame(V1 = "hsa:5", V2 = "CDS",
                                     V3 = "1:300..400", V4 = "GENEE; gene e"))
  }
  list(
    "/link/pathway/hsa" = link,
    "/list/pathway/hsa" = data.frame(
      V1 = c("hsa00010", "hsa00020"),
      V2 = c("Glycolysis / Gluconeogenesis - Homo sapiens (human)",
             "Citrate cycle (TCA cycle) - Homo sapiens (human)"),
      stringsAsFactors = FALSE
    ),
    "/list/hsa" = genes
  )
}

mock_rest <- function(tables) {
  function(path) {
    if (!path %in% names(tables)) stop("unexpected KEGG path: ", path)
    tables[[path]]
  }
}

refuse_rest <- function(path) stop("network must not be touched here")

live_table <- function() {
  data.frame(
    gs_name = "OLD_SET", gene_symbol = "G1", gs_description = "old",
    gs_id = "hsa99999", gs_source = "KEGG REST 2026-01-01",
    stringsAsFactors = FALSE
  )
}

kegg_memo_key <- "msig::kegg::Homo sapiens"

# ---- symbol parsing ----------------------------------------------------

test_that("parse_kegg_symbol takes the first symbol before the semicolon", {
  out <- parse_kegg_symbol(c(
    "GENEA, ALIASA; gene a",
    "GENEB; gene b",
    "uncharacterized protein",
    "; odd but symbol-free"
  ))
  expect_equal(out, c("GENEA", "GENEB", NA, NA))
})

test_that("kegg_org_code maps supported organisms and rejects others", {
  expect_equal(kegg_org_code("Homo sapiens"), "hsa")
  expect_equal(kegg_org_code("Mus musculus"), "mmu")
  expect_error(kegg_org_code("Danio rerio"), "supports")
})

# ---- live fetch assembly ----------------------------------------------

test_that("fetch_kegg_live assembles a cache-contract table", {
  testthat::local_mocked_bindings(kegg_rest_table = mock_rest(fake_rest()))
  df <- fetch_kegg_live("Homo sapiens")

  expect_true(all(c("gs_name", "gene_symbol", "gs_description",
                    "gs_id", "gs_source") %in% colnames(df)))
  # The organism suffix goes; parentheses inside the name survive.
  expect_setequal(unique(df$gs_name),
                  c("Glycolysis / Gluconeogenesis", "Citrate cycle (TCA cycle)"))
  # hsa:3 has no symbol and must be dropped, not kept as NA.
  expect_setequal(df$gene_symbol[df$gs_id == "hsa00020"], c("GENEB", "GENED"))
  expect_setequal(df$gene_symbol[df$gs_id == "hsa00010"], c("GENEA", "GENEB"))
  expect_true(all(startsWith(df$gs_source, "KEGG REST ")))
})

test_that("fetch_kegg_live fails loudly on an unexpected response shape", {
  tables <- fake_rest()
  tables[["/list/hsa"]] <- tables[["/list/hsa"]][, 1:2]
  testthat::local_mocked_bindings(kegg_rest_table = mock_rest(tables))
  expect_error(fetch_kegg_live("Homo sapiens"), "response shape")
})

# ---- refresh_geneset_cache --------------------------------------------

test_that("refresh writes a cache the normal read path then serves", {
  skip_if_not_installed("qs2")
  local_cache_env()
  testthat::local_mocked_bindings(kegg_rest_table = mock_rest(fake_rest()))

  # A stale session memo must not survive the refresh.
  cache_set(kegg_memo_key, live_table())
  withr::defer(cache_drop(kegg_memo_key))

  res <- refresh_geneset_cache("kegg", "Hs", quiet = TRUE)
  expect_equal(res$action, "refreshed")
  expect_equal(res$n_sets, 2L)

  served <- fetch_msigdbr_table("kegg", "Hs")
  expect_true(all(c("GENEA", "GENEB", "GENED") %in% served$gene_symbol))
  expect_false("OLD_SET" %in% served$gs_name)
})

test_that("a fresh cache is kept unless force = TRUE", {
  skip_if_not_installed("qs2")
  local_cache_env()
  withr::defer(cache_drop(kegg_memo_key))
  testthat::local_mocked_bindings(kegg_rest_table = mock_rest(fake_rest()))
  refresh_geneset_cache("kegg", "Hs", quiet = TRUE)

  testthat::local_mocked_bindings(
    kegg_rest_table = mock_rest(fake_rest(extra_gene = TRUE))
  )
  expect_equal(refresh_geneset_cache("kegg", "Hs", quiet = TRUE)$action, "fresh")

  res <- refresh_geneset_cache("kegg", "Hs", force = TRUE, quiet = TRUE)
  expect_equal(res$action, "refreshed")
  expect_true("GENEE" %in% fetch_msigdbr_table("kegg", "Hs")$gene_symbol)
})

test_that("a failed fetch keeps the previous cache intact", {
  skip_if_not_installed("qs2")
  dir <- local_cache_env()
  withr::defer(cache_drop(kegg_memo_key))
  qs2::qs_save(live_table(), file.path(dir, "kegg__Homo_sapiens.qs2"))

  testthat::local_mocked_bindings(kegg_rest_table = refuse_rest)
  expect_warning(
    res <- refresh_geneset_cache("kegg", "Hs", force = TRUE, quiet = TRUE),
    "keeping the previous cache"
  )
  expect_equal(res$action, "failed")
  expect_equal(read_geneset_cache("kegg", "Homo sapiens")$gs_name, "OLD_SET")
})

test_that("the msigdbr snapshot path trims columns and stamps provenance", {
  skip_if_not_installed("msigdbr")
  fake_msig <- data.frame(
    gs_name = "HALLMARK_X", gene_symbol = "G1", gs_description = "d",
    gs_id = "M1", gs_pmid = "irrelevant", stringsAsFactors = FALSE
  )
  testthat::local_mocked_bindings(fetch_msigdbr_raw = function(...) fake_msig)
  df <- fetch_msigdbr_snapshot("hallmark", "Homo sapiens")
  expect_false("gs_pmid" %in% colnames(df))
  expect_match(df$gs_source[[1L]], "^MSigDB msigdbr ")
})

# ---- TTL auto-refresh --------------------------------------------------

test_that("a stale live-sourced cache refreshes itself on next use", {
  skip_if_not_installed("qs2")
  dir <- local_cache_env()
  cache_drop(kegg_memo_key)
  withr::defer(cache_drop(kegg_memo_key))
  path <- file.path(dir, "kegg__Homo_sapiens.qs2")
  qs2::qs_save(live_table(), path)
  Sys.setFileTime(path, Sys.time() - 40 * 86400)

  testthat::local_mocked_bindings(kegg_rest_table = mock_rest(fake_rest()))
  expect_message(df <- fetch_msigdbr_table("kegg", "Hs"), "refreshing")
  expect_true("GENEA" %in% df$gene_symbol)
  expect_false("OLD_SET" %in% df$gs_name)
  expect_lt(geneset_cache_age_days(path), 1)
})

test_that("OMICSCORE_GENESET_TTL_DAYS=0 disables the auto-refresh", {
  skip_if_not_installed("qs2")
  dir <- local_cache_env()
  cache_drop(kegg_memo_key)
  withr::defer(cache_drop(kegg_memo_key))
  withr::local_envvar(OMICSCORE_GENESET_TTL_DAYS = "0")
  path <- file.path(dir, "kegg__Homo_sapiens.qs2")
  qs2::qs_save(live_table(), path)
  Sys.setFileTime(path, Sys.time() - 400 * 86400)

  testthat::local_mocked_bindings(kegg_rest_table = refuse_rest)
  df <- fetch_msigdbr_table("kegg", "Hs")
  expect_equal(unique(df$gs_name), "OLD_SET")
})

test_that("a prewarmed msigdbr cache never triggers a network call", {
  skip_if_not_installed("qs2")
  dir <- local_cache_env()
  cache_drop(kegg_memo_key)
  withr::defer(cache_drop(kegg_memo_key))
  prewarmed <- live_table()
  prewarmed$gs_source <- NULL
  path <- file.path(dir, "kegg__Homo_sapiens.qs2")
  qs2::qs_save(prewarmed, path)
  Sys.setFileTime(path, Sys.time() - 400 * 86400)

  testthat::local_mocked_bindings(kegg_rest_table = refuse_rest)
  df <- fetch_msigdbr_table("kegg", "Hs")
  expect_equal(unique(df$gs_name), "OLD_SET")
})

# ---- cache directory resolution ---------------------------------------

test_that("the env var wins, the per-user dir is the fallback", {
  tmp <- withr::local_tempdir()
  withr::local_envvar(OMICSCORE_GENESET_CACHE = NA, R_USER_CACHE_DIR = tmp)

  expect_null(geneset_cache_dir())          # nothing created yet
  wdir <- geneset_cache_dir(write = TRUE)   # creates the per-user dir
  expect_true(dir.exists(wdir))
  expect_true(startsWith(wdir, tmp))
  expect_equal(geneset_cache_dir(), wdir)   # and reads now find it

  explicit <- file.path(tmp, "explicit")
  dir.create(explicit)
  withr::local_envvar(OMICSCORE_GENESET_CACHE = explicit)
  expect_equal(geneset_cache_dir(), explicit)
})

# ---- status and provenance --------------------------------------------

test_that("geneset_cache_status reports source and age", {
  skip_if_not_installed("qs2")
  dir <- local_cache_env()
  qs2::qs_save(live_table(), file.path(dir, "kegg__Homo_sapiens.qs2"))

  status <- geneset_cache_status(c("kegg", "hallmark"), "Hs")
  kegg <- status[status$database == "kegg", ]
  expect_true(kegg$cached)
  expect_equal(kegg$n_sets, 1L)
  expect_equal(kegg$source, "KEGG REST 2026-01-01")
  expect_false(status[status$database == "hallmark", ]$cached)
})

test_that("geneset_table_source reads provenance off the session memo", {
  key <- "msig::hallmark::Homo sapiens"
  cache_drop(key)
  expect_equal(geneset_table_source("hallmark", "Homo sapiens"),
               "MSigDB (msigdbr)")
  cache_set(key, live_table())
  withr::defer(cache_drop(key))
  expect_equal(geneset_table_source("hallmark", "Homo sapiens"),
               "KEGG REST 2026-01-01")
})
