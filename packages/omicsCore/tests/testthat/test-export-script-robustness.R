# The promise `export_script()` makes is that its output runs. The
# round-trip test in test-export-script.R proves that for one shaped
# project; these prove the properties that have to hold for *any*
# project, including the shapes a hand-built or partially-restored one
# can take. A script that does not parse is at least loud; one that
# parses and quietly omits an argument is the failure worth hunting.

rs_input <- function(omics_type = "proteomics", source_path = "raw/x.xlsx",
                     assay_type = if (omics_type == "rnaseq") "logcpm"
                                  else "normalized_intensity") {
  mat <- matrix(as.numeric(1:24), nrow = 6,
                dimnames = list(paste0("g", 1:6), paste0("s", 1:4)))
  meta <- data.frame(group = c("A", "A", "B", "B"),
                     row.names = paste0("s", 1:4))
  feat <- data.frame(feature_id = paste0("g", 1:6))
  omics_input(mat, meta, feat, omics_type = omics_type,
              assay_type = assay_type, source_path = source_path)
}

rs_project <- function(experiments = list(proteomics = rs_input()),
                       bundles = list()) {
  proj <- omics_project("robustness", experiments = experiments)
  proj$bundles <- bundles
  proj
}

parses <- function(lines) {
  isTRUE(tryCatch({ parse(text = paste(lines, collapse = "\n")); TRUE },
                  error = function(e) FALSE))
}

# ---- the output is always R -------------------------------------------

test_that("every project shape yields a script that parses", {
  shapes <- list(
    "no experiments, no bundles" = rs_project(list()),
    "experiment only"            = rs_project(),
    "no archived file"           = rs_project(
      list(proteomics = rs_input(source_path = NULL))),
    "bundles without experiments" = rs_project(list(), list(
      diff = new_analysis_bundle("run_diff", list(omics_type = "proteomics"),
                                 params = list(method = "limma")))),
    "two layers" = rs_project(list(
      proteomics = rs_input("proteomics"),
      rnaseq = rs_input("rnaseq", assay_type = "raw_count"))),
    "full pipeline" = rs_project(
      list(proteomics = rs_input("proteomics"), rnaseq = rs_input("rnaseq")),
      list(
        qc = new_analysis_bundle("run_qc", list(omics_type = "proteomics"),
               params = list(impute_method = "half_min")),
        diff = new_analysis_bundle("run_diff", list(omics_type = "proteomics"),
                 params = list(method = "limma", analysis_type = "group")),
        enrich = new_analysis_bundle("run_enrichment",
                   list(omics_type = "proteomics"),
                   params = list(type = "ora", database = "hallmark")),
        integration = new_analysis_bundle("run_integration",
                        list(omics_type = c("rnaseq", "proteomics")),
                        params = list(method = "concordance",
                                      experiments = c("rnaseq", "proteomics")))))
  )
  for (nm in names(shapes)) {
    expect_true(parses(export_script(shapes[[nm]])), info = nm)
    expect_true(parses(export_script(shapes[[nm]], include_plots = FALSE)),
                info = paste(nm, "(no plots)"))
  }
})

test_that("a project with analyses but no data defines what it references", {
  # `omics_project(experiments = list())` is legal, and a restored or
  # hand-built project can reach here. Emitting the run_* calls against
  # an `input` no line defines produces a script that fails deep inside
  # run_diff() with nothing to say about why.
  proj <- rs_project(list(), list(
    diff = new_analysis_bundle("run_diff", list(omics_type = "proteomics"),
                               params = list(method = "limma"))))
  lines <- export_script(proj)
  expect_true(any(grepl("^input <- ", lines)))
  expect_true(any(grepl("no imported data", lines)))
  expect_true(parses(lines))
})

# ---- argument fidelity -------------------------------------------------

test_that("a literal containing angle brackets survives", {
  # Cohort studies label groups things like "<30". An earlier version
  # pattern-matched the deparsed text for "<" and dropped the argument,
  # producing an incomplete call with only a NOTE to show for it.
  for (value in list("<30", ">=65", "a<b", "G2 vs G1", "50%", "x'y")) {
    expect_equal(render_value(value), deparse(value), info = value)
  }
})

test_that("literal vectors and numerics round-trip exactly", {
  for (value in list(c("age", "sex"), 0.05, 1L, TRUE, c(1.5, 2.5),
                     NA_character_, c(a = 1, b = 2))) {
    rendered <- render_value(value)
    expect_false(is.na(rendered))
    expect_identical(eval(parse(text = rendered)), value)
  }
})

test_that("a value that is not source becomes a note, not broken code", {
  proj <- rs_project(bundles = list(
    diff = new_analysis_bundle("run_diff", list(omics_type = "proteomics"),
             params = list(method = "limma", covariates = quote(f(x))))))
  lines <- export_script(proj)
  expect_true(any(grepl("not a literal", lines, fixed = TRUE)))
  expect_false(any(grepl("^\\s+covariates\\s+=", lines)))
  expect_true(parses(lines))
})

test_that("an argument named `name` is not rewritten by the integration block", {
  # The integration section used to append its comma by running a sub()
  # over every line, which would have edited any other line beginning
  # with a `name` argument.
  proj <- rs_project(
    list(proteomics = rs_input("proteomics"), rnaseq = rs_input("rnaseq")),
    list(
      diff = new_analysis_bundle("run_diff", list(omics_type = "proteomics"),
               params = list(method = "limma", group_col = "name")),
      integration = new_analysis_bundle("run_integration",
                      list(omics_type = c("rnaseq", "proteomics")),
                      params = list(method = "concordance"))))
  lines <- export_script(proj)
  expect_true(parses(lines))
  # Exactly one comma-terminated `name =` line: the one omics_project()
  # needs.
  expect_equal(sum(grepl("^  name\\s+= .*,$", lines)), 1L)
})

# ---- determinism -------------------------------------------------------

test_that("two exports of one project differ only in the timestamp", {
  proj <- rs_project(bundles = list(
    diff = new_analysis_bundle("run_diff", list(omics_type = "proteomics"),
                               params = list(method = "limma"))))
  drop_stamp <- function(l) l[!grepl("^# Written", l)]
  expect_identical(drop_stamp(export_script(proj)),
                   drop_stamp(export_script(proj)))
})

test_that("writing to a path returns the same lines it wrote", {
  path <- tempfile(fileext = ".R"); on.exit(unlink(path), add = TRUE)
  lines <- export_script(rs_project(), path)
  expect_identical(readLines(path), lines)
})
