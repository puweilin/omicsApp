# Names a person typed, in the files the package writes.
#
# A project name, a layer tag, a group label, a feature id: each is
# text the user chose, and each ends up in a report and a script.
# Markdown passes raw HTML through, so `<script>` in a project name
# ran when the report was opened; a newline in a name ended a comment
# in the script and handed the rest to the parser; a tag with a space
# became a variable name R cannot read. Each is now carried as text.

skip_if_not_installed("limma")

hostile_project <- function() {
  inp <- realistic_input("proteomics")
  feat <- '<img src=x onerror=alert(1)>'
  samp <- "<i>S01</i>"
  rownames(inp$expr_mat)[1] <- feat
  inp$feature_df$feature_id[1] <- feat
  colnames(inp$expr_mat)[1] <- samp
  rownames(inp$meta_df)[1] <- samp
  inp$meta_df$group <- as.character(inp$meta_df$group)
  inp$meta_df$group[inp$meta_df$group == "G1"] <- "<b>G1</b>"
  proj <- omics_project(name = '<script>alert("project")</script>',
                        experiments = list(`<u>prot</u>` = inp))
  proj$bundles <- list(
    `<script>alert("bundle")</script>` = run_diff(
      inp, method = "ttest", analysis_type = "group", group_col = "group",
      control_group = "<b>G1</b>", case_group = "G2"))
  proj
}

test_that("the report shows a hostile name as text", {
  skip_if_not_installed("rmarkdown")
  skip_if_not(rmarkdown::pandoc_available(), "pandoc is not available")
  proj <- hostile_project()
  out <- tempfile(fileext = ".html")
  withr::defer(unlink(out))
  # pandoc reports the `<img src=x>` it can no longer find as a
  # resource, which is the point: it is text now.
  suppressWarnings(export_report(proj, out))
  html <- paste(readLines(out, warn = FALSE), collapse = "\n")
  for (raw in c('<script>alert("project")', '<script>alert("bundle")',
                "<img src=x", "<i>S01</i>", "<b>G1</b>", "<u>prot</u>")) {
    expect_false(grepl(raw, html, fixed = TRUE), label = raw)
  }
  expect_true(grepl("&lt;script&gt;alert", html, fixed = TRUE))
  expect_true(grepl("&lt;u&gt;prot&lt;/u&gt;", html, fixed = TRUE))
})

test_that("a hostile name survives the .omp round trip unchanged", {
  proj <- hostile_project()
  path <- tempfile(fileext = ".omp")
  withr::defer(unlink(path))
  save_project(proj, path)
  back <- load_project(path)
  expect_identical(back$name, proj$name)
  expect_identical(names(back$experiments), names(proj$experiments))
  expect_identical(names(back$bundles), names(proj$bundles))
})

script_parses <- function(proj) {
  lines <- export_script(proj)
  parsed <- tryCatch(parse(text = lines), error = function(e) e)
  expect_false(inherits(parsed, "error"),
               label = paste("script parses for", encodeString(proj$name)))
  lines
}

test_that("the script parses whatever the names carry", {
  inp <- realistic_input("proteomics")
  diff <- run_diff(inp, method = "ttest", analysis_type = "group", group_col = "group",
                   control_group = "G1", case_group = "G2")
  names_try <- c(quote = 'Pro"ject', backslash = "C:\\path\\proj",
                 newline = "line1\nline2", carriage = "a\r\nb", single = "it's",
                 unicode = "\u9879\u76ee", backtick = "a`b", dollar = "cost $5",
                 percent = "100%", hash = "# not a comment", tab = "a\tb")
  for (nm in names(names_try)) {
    proj <- omics_project(name = names_try[[nm]],
                          experiments = stats::setNames(list(inp), names_try[[nm]]))
    proj$bundles <- list(diff = diff)
    lines <- script_parses(proj)
    # The comment header still carries the name, on one line.
    expect_true(any(grepl("^# Project : ", lines)))
    expect_false(any(grepl("^line2|^b$", lines)))
  }
})

test_that("a project with two layers gets variables R can read back", {
  inp <- realistic_input("proteomics")
  rna <- realistic_input("rnaseq")
  proj <- omics_project("two", experiments = list(`rna seq` = rna, `prot-eomics` = inp,
                                                   `<u>x</u>` = inp))
  proj$bundles <- list(
    integration = new_analysis_bundle("run_integration",
                                      params = list(method = "concordance",
                                                    experiments = c("rna seq", "prot-eomics"))))
  lines <- script_parses(proj)
  vars <- sub(" <- read_omics\\($", "", grep(" <- read_omics\\($", lines, value = TRUE))
  expect_length(vars, 3L)
  expect_identical(make.names(vars), vars)
  expect_false(anyDuplicated(vars) > 0L)
  # The experiments list names each layer as it was called.
  exp_line <- grep("^  experiments = list\\(", lines, value = TRUE)
  expect_match(exp_line, "`rna seq` = input_rna.seq", fixed = TRUE)
  expect_match(exp_line, "`<u>x</u>` = ", fixed = TRUE)
  # And evaluates to a list keyed by those tags.
  env <- new.env()
  for (v in vars) assign(v, NULL, envir = env)
  assign("omics_project", function(name, experiments) names(experiments), envir = env)
  keys <- eval(parse(text = grep("^project <- omics_project\\(", lines) |>
                       (\(i) lines[i:(i + 3L)])() |> paste(collapse = "\n") |>
                       sub(pattern = "^project <- ", replacement = "")),
               envir = env)
  expect_setequal(keys, names(proj$experiments))
})
