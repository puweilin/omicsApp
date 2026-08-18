# The on-disk gene-set cache exists to remove the ~10s msigdbr lazy-load
# from the first enrichment in a fresh process. It is read on a path
# where a wrong answer would be far worse than a slow one, so every
# failure mode here has to degrade to "no cache" rather than to bad data.

with_cache_dir <- function(code, files = list(), set_env = TRUE) {
  dir <- file.path(tempdir(), paste0("gs-", as.integer(runif(1, 1, 1e9))))
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  for (nm in names(files)) {
    value <- files[[nm]]
    path <- file.path(dir, nm)
    if (is.character(value)) writeLines(value, path) else qs2::qs_save(value, path)
  }
  old <- Sys.getenv("OMICSCORE_GENESET_CACHE", unset = NA)
  if (set_env) Sys.setenv(OMICSCORE_GENESET_CACHE = dir) else
    Sys.unsetenv("OMICSCORE_GENESET_CACHE")
  on.exit({
    if (is.na(old)) Sys.unsetenv("OMICSCORE_GENESET_CACHE")
    else Sys.setenv(OMICSCORE_GENESET_CACHE = old)
    unlink(dir, recursive = TRUE)
  }, add = TRUE)
  force(code)
}

fake_table <- function(n = 30L) {
  data.frame(
    gs_name = rep(c("SET_A", "SET_B"), each = n / 2L),
    gene_symbol = paste0("GENE", seq_len(n)),
    gs_description = rep(c("Set A", "Set B"), each = n / 2L),
    gs_id = rep(c("M1", "M2"), each = n / 2L),
    stringsAsFactors = FALSE
  )
}

test_that("the cache file name encodes database and organism", {
  path <- geneset_cache_file("/cache", "go_bp", "Homo sapiens")
  expect_equal(basename(path), "go_bp__Homo_sapiens.qs2")
})

test_that("no cache directory means no cache", {
  with_cache_dir(expect_null(read_geneset_cache("go_bp", "Homo sapiens")),
                 set_env = FALSE)
})

test_that("a directory that does not exist means no cache", {
  old <- Sys.getenv("OMICSCORE_GENESET_CACHE", unset = NA)
  Sys.setenv(OMICSCORE_GENESET_CACHE = file.path(tempdir(), "definitely-absent"))
  on.exit(if (is.na(old)) Sys.unsetenv("OMICSCORE_GENESET_CACHE")
          else Sys.setenv(OMICSCORE_GENESET_CACHE = old), add = TRUE)
  expect_null(read_geneset_cache("go_bp", "Homo sapiens"))
})

test_that("a missing file for this database means no cache", {
  skip_if_not_installed("qs2")
  with_cache_dir(
    expect_null(read_geneset_cache("go_mf", "Homo sapiens")),
    files = list("go_bp__Homo_sapiens.qs2" = fake_table())
  )
})

test_that("a cached table is returned when it is there", {
  skip_if_not_installed("qs2")
  with_cache_dir({
    df <- read_geneset_cache("go_bp", "Homo sapiens")
    expect_s3_class(df, "data.frame")
    expect_equal(nrow(df), 30L)
    expect_true(all(c("gs_name", "gene_symbol") %in% colnames(df)))
  }, files = list("go_bp__Homo_sapiens.qs2" = fake_table()))
})

test_that("a corrupt file degrades to no cache rather than an error", {
  skip_if_not_installed("qs2")
  with_cache_dir({
    # Truncated or non-qs2 content must not propagate an error into an
    # enrichment run; falling back costs 10s, erroring costs the run.
    expect_null(read_geneset_cache("go_bp", "Homo sapiens"))
  }, files = list("go_bp__Homo_sapiens.qs2" = "not a qs2 archive"))
})

test_that("a table missing the required columns is ignored", {
  skip_if_not_installed("qs2")
  with_cache_dir({
    # Guards against a cache written by an older prewarm whose column
    # contract has since changed.
    expect_null(read_geneset_cache("go_bp", "Homo sapiens"))
  }, files = list("go_bp__Homo_sapiens.qs2" =
                    data.frame(something_else = 1:3)))
})

test_that("an empty table is ignored", {
  skip_if_not_installed("qs2")
  with_cache_dir(
    expect_null(read_geneset_cache("go_bp", "Homo sapiens")),
    files = list("go_bp__Homo_sapiens.qs2" = fake_table()[0, ])
  )
})

test_that("fetch_msigdbr_table prefers the cache over msigdbr", {
  skip_if_not_installed("qs2")
  with_cache_dir({
    # Distinguishable from anything msigdbr would return, so a pass
    # proves the disk path was taken rather than a live fetch.
    df <- fetch_msigdbr_table("go_bp", "Hs")
    expect_equal(nrow(df), 30L)
    expect_equal(sort(unique(df$gs_name)), c("SET_A", "SET_B"))
  }, files = list("go_bp__Homo_sapiens.qs2" = fake_table()))
})

test_that("the cached table still satisfies every downstream reader", {
  skip_if_not_installed("qs2")
  skip_if_not_installed("clusterProfiler")
  with_cache_dir({
    # The prewarm keeps four columns out of msigdbr's twenty. These are
    # the three functions that read the table; if the contract were too
    # narrow, one of them would fail here rather than in production.
    sets <- get_gene_set_list("go_bp", "Hs", min_size = 1L, max_size = 100L)
    expect_length(sets, 2L)

    terms <- build_term_tables("go_bp", "Hs")
    expect_equal(nrow(terms$term2gene), 30L)
    expect_true(all(c("SET_A", "SET_B") %in% terms$term2name$term))

    listed <- list_gene_sets("go_bp", "Hs")
    expect_equal(nrow(listed), 30L)
    expect_equal(sort(unique(listed$pathway_id)), c("M1", "M2"))
  }, files = list("go_bp__Homo_sapiens.qs2" = fake_table()))
})
