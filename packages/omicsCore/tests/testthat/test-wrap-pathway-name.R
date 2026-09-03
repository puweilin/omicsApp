# Pathway labels used to be cut at a character count. A truncated GO BP
# term is often indistinguishable from three others sharing its opening
# words, which is the one thing a pathway label must not be.

test_that("a name that fits is left alone", {
  expect_identical(wrap_pathway_name("Inflammatory response"),
                   "Inflammatory response")
  expect_false(grepl("\n", wrap_pathway_name("Apoptosis"), fixed = TRUE))
})

test_that("a long name is wrapped, not shortened", {
  long <- paste("Positive regulation of cytokine production involved in",
                "immune response")
  out <- wrap_pathway_name(long)
  expect_true(grepl("\n", out, fixed = TRUE))
  # Every word survives: that is the whole point of wrapping over
  # truncating.
  expect_identical(strsplit(gsub("\n", " ", out), " +")[[1]],
                   strsplit(long, " +")[[1]])
})

test_that("no line runs past the requested width", {
  long <- strrep("alpha beta gamma ", 6)
  for (w in c(20L, 34L, 50L)) {
    lines <- strsplit(wrap_pathway_name(long, width = w, max_lines = 99L),
                      "\n", fixed = TRUE)[[1]]
    expect_true(all(nchar(lines) <= w), info = paste("width", w))
  }
})

test_that("a runaway name is capped, and says it was", {
  # Without the cap one 200-character Reactome term claims the height of
  # four others and the panel becomes about the label.
  out <- wrap_pathway_name(strrep("word ", 100), width = 20L, max_lines = 3L)
  lines <- strsplit(out, "\n", fixed = TRUE)[[1]]
  expect_length(lines, 3L)
  expect_match(lines[3], "…$")
})

test_that("a name exactly at the cap is not marked as elided", {
  # "one two" / "three" / "four" -- three lines, which is the cap. The
  # ellipsis has to mean "there was more", so it must not appear here.
  out <- wrap_pathway_name("one two three four", width = 9L, max_lines = 3L)
  expect_length(strsplit(out, "\n", fixed = TRUE)[[1]], 3L)
  expect_false(grepl("…", out, fixed = TRUE))
})

test_that("missing and empty names pass through rather than erroring", {
  # These reach the plot from a standardized table that permits them.
  expect_true(is.na(wrap_pathway_name(NA_character_)))
  expect_identical(wrap_pathway_name(""), "")
  expect_length(wrap_pathway_name(character(0)), 0L)
})

test_that("a vector is handled element-wise, with no names attached", {
  out <- wrap_pathway_name(c(a = "short", b = strrep("long ", 20)))
  expect_null(names(out))
  expect_length(out, 2L)
  expect_false(grepl("\n", out[1], fixed = TRUE))
  expect_true(grepl("\n", out[2], fixed = TRUE))
})

test_that("gene symbols still truncate, since a symbol cannot wrap", {
  # plot-integration.R labels points with feature symbols; splitting one
  # over two lines would be worse than shortening it.
  expect_match(truncate_pathway_name(strrep("A", 60L)), "…$")
  expect_identical(truncate_pathway_name("TP53"), "TP53")
})
