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

test_that("the split survives a child process whose locale is C", {
  skip_on_cran()
  skip_if_not_installed("callr")
  # The function travels as source text, so the child needs no library:
  # under R CMD check the package is installed, under devtools::test()
  # it is loaded from source, and neither is what is being tested here.
  got <- callr::r(
    function(src, name) {
      split <- eval(parse(text = src))
      p <- split(name)
      list(stem = charToRaw(p$stem), ext = p$ext, locale = Sys.getlocale("LC_CTYPE"))
    },
    args = list(deparse(omicsApp:::split_file_name), "蛋白组数据.xlsx"),
    env = c(callr::rcmd_safe_env(), LANG = "C", LC_ALL = "C", LC_CTYPE = "C"),
    timeout = 120
  )
  expect_identical(got$locale, "C")
  expect_identical(got$ext, "xlsx")
  expect_identical(got$stem, charToRaw("蛋白组数据"))
})
