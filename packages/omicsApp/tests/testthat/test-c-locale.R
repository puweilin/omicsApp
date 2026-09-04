# The project store in a process whose locale is C.
#
# test-project-store.R proves project_slug() keeps a non-Latin name in
# this process, whatever its locale. This runs the whole store -- slug,
# save, list, archive the upload, open -- in a child R process started
# with LC_ALL=C, which is what a container without a configured locale
# gives the app. Strings are \u escapes so the source stays ASCII.

# Where each package this process runs came from, so the child runs the
# same code: the source trees under load_all, the libraries otherwise.
package_origin <- function(pkg) {
  root <- system.file(package = pkg)
  if (identical(basename(root), "inst")) root <- dirname(root)
  has_sources <- dir.exists(file.path(root, "R")) &&
    length(list.files(file.path(root, "R"), pattern = "\\.[Rr]$")) > 0L
  if (file.exists(file.path(root, "DESCRIPTION")) && has_sources &&
      requireNamespace("pkgload", quietly = TRUE)) {
    return(list(kind = "source", path = normalizePath(root)))
  }
  list(kind = "installed", path = dirname(normalizePath(root)))
}

store_scenario <- function(core, app, store, name, upload, upload_name) {
  load <- function(origin, pkg) {
    if (identical(origin$kind, "source")) pkgload::load_all(origin$path, quiet = TRUE)
    else library(pkg, character.only = TRUE, lib.loc = c(origin$path, .libPaths()))
  }
  load(core, "omicsCore")
  load(app, "omicsApp")
  Sys.setenv(OMICSAPP_DATA_DIR = store)
  ns <- asNamespace("omicsApp")
  slug <- ns$project_slug(name)

  mat <- matrix(as.numeric(1:12), 3, 4,
                dimnames = list(paste0("g", 1:3), paste0("s", 1:4)))
  meta <- data.frame(group = c("A", "A", "B", "B"), row.names = paste0("s", 1:4))
  inp <- omicsCore::omics_input(mat, meta, data.frame(feature_id = rownames(mat)),
                                omics_type = "proteomics",
                                assay_type = "normalized_intensity")
  project <- omicsCore::omics_project(name, experiments = list(proteomics = inp))
  saved <- ns$store_save_project(project, slug)
  listed <- ns$list_saved_projects()
  archived <- ns$store_raw_upload(upload, upload_name, "abcdef0123456789:1")
  opened <- ns$store_load_project(slug)
  list(
    locale = Sys.getlocale("LC_CTYPE"),
    slug = slug,
    saved_ok = isTRUE(saved$ok),
    saved_file_exists = file.exists(saved$path),
    listed_slugs = listed$slug,
    archived_ok = isTRUE(archived$ok),
    archived_name = basename(archived$path),
    opened_ok = isTRUE(opened$ok),
    opened_message = opened$message,
    saved_message = saved$message,
    saved_path = saved$path,
    opened_name = opened$project$name,
    files = list.files(store, recursive = TRUE)
  )
}

test_that("a project named in Chinese is saved, listed, archived and reopened under LC_ALL=C", {
  skip_on_cran()
  skip_if_not_installed("callr")
  skip_if_not_installed("qs2")
  skip_if_not_installed("withr")
  store <- withr::local_tempdir()
  upload <- withr::local_tempfile(fileext = ".xlsx")
  writeBin(as.raw(1:64), upload)
  name <- "蛋白组 G2 vs G1"
  upload_name <- "蛋白组数据.xlsx"

  got <- callr::r(
    store_scenario,
    args = list(package_origin("omicsCore"), package_origin("omicsApp"),
                store, name, upload, upload_name),
    env = c(callr::rcmd_safe_env(), LANG = "C", LC_ALL = "C", LC_CTYPE = "C"),
    timeout = 300
  )

  # Every comparison below is on bytes. This test process may itself be
  # running in a C locale (a bare shell, a CI runner), where comparing
  # or opening a UTF-8 path by name would fail for the runner's reasons
  # rather than the store's.
  utf8_bytes <- function(x) charToRaw(enc2utf8(x))
  expect_identical(got$locale, "C")
  expect_identical(charToRaw(got$slug), utf8_bytes("\u86cb\u767d\u7ec4_G2_vs_G1"))
  expect_true(got$saved_ok)
  expect_true(got$saved_file_exists)
  expect_length(got$listed_slugs, 1L)
  expect_identical(charToRaw(got$listed_slugs), charToRaw(got$slug))
  expect_true(got$archived_ok)
  expect_identical(utils::head(charToRaw(got$archived_name), 12L),
                   utils::head(utf8_bytes(upload_name), 12L))
  expect_true(grepl("__[0-9a-f]{12}\\.xlsx$", got$archived_name))
  expect_true(got$opened_ok)
  expect_identical(charToRaw(got$opened_name), utf8_bytes(name))

  # And the files are on disk under their real names, readable from here
  on_disk <- list.files(store, recursive = TRUE, full.names = TRUE)
  omp <- on_disk[grepl("\\.omp$", on_disk)]
  expect_length(omp, 1L)
  expect_identical(charToRaw(basename(omp)), utf8_bytes(paste0(got$slug, ".omp")))
  reopened <- omicsCore::load_project(omp)
  expect_identical(charToRaw(reopened$name), utf8_bytes(name))
})
