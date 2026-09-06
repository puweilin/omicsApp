# The archive names an upload after the file the user chose, and until
# R 4.5 that went through tools::file_ext(), which now calls basename()
# and therefore fails on a non-ASCII name in a C locale. The split is
# done here instead; these pin its rules to tools' own.

test_that("split_file_name() agrees with tools on ordinary names", {
  split <- omicsApp:::split_file_name
  for (name in c("data.xlsx", "counts.tsv", "archive.tar.gz", "noext",
                 "trailing.", "dots.in.name.csv", "x.a1")) {
    expect_identical(split(name)$ext, tools::file_ext(name), info = name)
    expect_identical(split(name)$stem, tools::file_path_sans_ext(name), info = name)
  }
  # A name that is nothing but an extension has none: R >= 4.5's rule,
  # which older tools::file_ext() did not apply.
  expect_identical(split(".hidden"), list(stem = ".hidden", ext = ""))
})

test_that("split_file_name() keeps a UTF-8 name whole and marked", {
  name <- "蛋白组数据.xlsx"
  parts <- omicsApp:::split_file_name(name)
  expect_identical(parts$ext, "xlsx")
  expect_identical(charToRaw(parts$stem), charToRaw("蛋白组数据"))
  expect_identical(Encoding(parts$stem), "UTF-8")
})

test_that("archiving a UTF-8 name in a C locale does not go through basename()", {
  skip_on_cran()
  skip_if_not_installed("callr")
  skip_if_not_installed("withr")
  # The failure this guards against needs a C locale on the child and
  # R >= 4.5 to reproduce; on older R the assertion still holds.
  got <- callr::r(
    function(name) {
      src <- Sys.getenv("OMICSAPP_SRC")
      if (nzchar(src)) pkgload::load_all(src, quiet = TRUE) else library(omicsApp)
      p <- asNamespace("omicsApp")$split_file_name(name)
      list(stem = charToRaw(p$stem), ext = p$ext, locale = Sys.getlocale("LC_CTYPE"))
    },
    args = list("蛋白组数据.xlsx"),
    env = c(callr::rcmd_safe_env(), LANG = "C", LC_ALL = "C", LC_CTYPE = "C",
            OMICSAPP_SRC = if (nzchar(system.file("DESCRIPTION", package = "omicsApp")) &&
                                dir.exists(file.path(dirname(system.file("DESCRIPTION", package = "omicsApp")), "R")))
              dirname(system.file("DESCRIPTION", package = "omicsApp")) else ""),
    timeout = 120
  )
  expect_identical(got$locale, "C")
  expect_identical(got$ext, "xlsx")
  expect_identical(got$stem, charToRaw("蛋白组数据"))
})
