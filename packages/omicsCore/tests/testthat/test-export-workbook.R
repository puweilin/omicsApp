# export_workbook(): the tables, one sheet each.
#
# What replaced the rendered report. The questions worth asking of it are
# what a reader would notice: is my analysis in there, are the rows all
# there, and does the sheet say which gene each row is.

wb_input <- function(n_features = 20L, n_samples = 8L, na_frac = 0.1) {
  set.seed(5)
  ids <- paste0("F", seq_len(n_features))
  sids <- paste0("s", seq_len(n_samples))
  m <- matrix(stats::rnorm(n_features * n_samples, 20, 2),
              nrow = n_features, dimnames = list(ids, sids))
  if (na_frac > 0) m[sample(length(m), floor(length(m) * na_frac))] <- NA
  omics_input(m,
              data.frame(group = rep(c("a", "b"), length.out = n_samples),
                         row.names = sids),
              data.frame(feature_id = ids, feature_symbol = ids,
                         row.names = ids),
              omics_type = "proteomics",
              assay_type = "normalized_intensity")
}

wb_project <- function(impute = "MinProb", with_diff = TRUE) {
  x <- wb_input()
  p <- omics_project("wb", experiments = list(prot = x))
  p$bundles <- list(qc = run_qc(x, impute_method = impute))
  if (with_diff) {
    p$bundles$diff <- run_diff(x, method = "ttest", analysis_type = "group",
                               group_col = "group", control_group = "a",
                               case_group = "b")
  }
  p
}

test_that("each analysis that ran contributes a sheet", {
  skip_if_not_installed("openxlsx")
  p <- wb_project()
  s <- workbook_sheets(p)
  expect_true(all(c("Analysis parameters", "QC samples", "QC features",
                    "Differential") %in% names(s)))
  # Nothing was integrated or enriched, so there is no sheet claiming so.
  expect_false(any(c("Enrichment", "Integration") %in% names(s)))
})

# A sheet named "Imputed data" that holds the imported values is worse
# than no sheet: it is an assertion that imputation ran.
test_that("the imputed matrix appears only when imputation ran", {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("imputeLCMD")
  expect_true("Imputed data" %in% names(workbook_sheets(wb_project("MinProb"))))
  expect_false("Imputed data" %in% names(workbook_sheets(wb_project("none"))))
})

test_that("the imputed sheet carries feature ids, not just numbers", {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("imputeLCMD")
  s <- workbook_sheets(wb_project("MinProb"))[["Imputed data"]]
  expect_identical(names(s)[[1L]], "feature_id")
  expect_true(all(grepl("^F[0-9]+$", s$feature_id)))
  # a matrix written straight out would have lost these
  expect_equal(nrow(s), 20L)
})

test_that("a project with no analyses has nothing to export", {
  skip_if_not_installed("openxlsx")
  p <- omics_project("empty", experiments = list(prot = wb_input()))
  # Not "just the parameters": a workbook describing a project nobody has
  # analysed reads as a successful export of nothing.
  expect_length(workbook_sheets(p), 0L)
  expect_error(export_workbook(p, tempfile(fileext = ".xlsx")),
               "nothing to export")
})

test_that("export_workbook writes a readable file at the path given", {
  skip_if_not_installed("openxlsx")
  p <- wb_project()
  # No extension: what Shiny's downloadHandler passes its content
  # function, and the shape that broke the previous export twice.
  out <- withr::local_tempfile()
  export_workbook(p, out)
  expect_true(file.exists(out))
  expect_true("Differential" %in% openxlsx::getSheetNames(out))

  # read.xlsx() insists on the extension where getSheetNames() does not,
  # so the round trip goes through a renamed copy. The handler serves the
  # file under a name of its own choosing, so this is the reader's rule,
  # not something the export has to satisfy.
  named <- withr::local_tempfile(fileext = ".xlsx")
  file.copy(out, named, overwrite = TRUE)
  back <- openxlsx::read.xlsx(named, sheet = "Differential")
  expect_equal(nrow(back), nrow(p$bundles$diff$results$diff_result_df))
})

test_that("export_workbook refuses to overwrite by default", {
  skip_if_not_installed("openxlsx")
  p <- wb_project()
  out <- withr::local_tempfile(fileext = ".xlsx")
  export_workbook(p, out)
  expect_error(export_workbook(p, out), "already exists")
  expect_no_error(export_workbook(p, out, overwrite = TRUE))
})

test_that("export_workbook rejects things that are not projects", {
  expect_error(export_workbook(list(), tempfile()), "omics_project")
  expect_error(workbook_sheets(list()), "omics_project")
})

# Excel stops on a sheet name over 31 characters or holding []:*?/\, and
# openxlsx passes that on as an error -- which would fail an entire export
# over a name this code chose.
test_that("sheet names are ones Excel will accept", {
  long <- paste(rep("a", 60), collapse = "")
  expect_lte(nchar(safe_sheet_name(long)), 31L)
  expect_false(grepl("[\\[\\]:*?/\\\\]", safe_sheet_name("a[b]:c*d?e/f\\g"),
                     perl = TRUE))
  # and two sheets cannot collide, case-insensitively as Excel compares
  first <- safe_sheet_name("Results")
  second <- safe_sheet_name("results", taken = first)
  expect_false(identical(tolower(first), tolower(second)))
  expect_lte(nchar(second), 31L)
})

test_that("the parameters sheet says what each analysis ran with", {
  skip_if_not_installed("openxlsx")
  p <- wb_project()
  s <- workbook_sheets(p)[["Analysis parameters"]]
  expect_true(all(c("analysis", "parameter", "value") %in% names(s)))
  # The contrast is the thing someone asks about three months later.
  expect_true(any(s$analysis == "diff" & s$parameter == "method"))
  expect_true(any(s$value == "ttest"))
  # and which layer produced it
  expect_true(any(grepl("^layer: ", s$analysis)))
})
