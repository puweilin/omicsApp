# omicsCore in a process whose locale is C.
#
# Container base images do not reliably set a UTF-8 locale, and the
# store code already goes out of its way to survive that. Nothing
# checked that the rest of the package does. Every string in this file
# is written with \u escapes, so the test source itself is ASCII and
# what the child sees is exactly the UTF-8 the user typed.
#
# The scenario runs in a child R process started with LC_ALL=C, using
# the helpers in helper-child-process.R. The files it reads are written
# here, in the parent, which is where a user's files come from too.

cjk_ids <- c("基因A", "基因B", "TP53", "MYC", "EGFR", "CDK1",
             paste0("G", 7:30))
cjk_samples <- c("样本1", "样本2", "S3", "S4", "S5", "S6")
cjk_project <- "项目 A"

write_cjk_files <- function() {
  set.seed(1)
  mat <- matrix(round(stats::rnorm(30 * 6, 20, 1), 3), 30, 6,
                dimnames = list(cjk_ids, cjk_samples))
  df <- data.frame(id = cjk_ids, mat, check.names = FALSE,
                   stringsAsFactors = FALSE)
  # Written as UTF-8 bytes, whatever locale this test process has:
  # write.csv() would translate through the native encoding, and under
  # a C locale that turns the ids into "invalid char string" before the
  # child ever sees them.
  csv <- tempfile(fileext = ".csv")
  lines <- c(paste(colnames(df), collapse = ","),
             apply(df, 1L, function(r) paste(r, collapse = ",")))
  con <- file(csv, open = "wb")
  writeBin(charToRaw(paste0(paste(enc2utf8(lines), collapse = "\n"), "\n")), con)
  close(con)
  xlsx <- tempfile(fileext = ".xlsx")
  # openxlsx writes UTF-8 correctly in any locale but grumbles about the
  # conversion when this process itself runs in C; that is the runner's
  # locale, not the package's.
  suppressWarnings(openxlsx::write.xlsx(df, xlsx))
  list(csv = csv, xlsx = xlsx)
}

c_locale_scenario <- function(origin, csv, xlsx, ids, samples, project_name) {
  if (identical(origin$kind, "source")) {
    pkgload::load_all(origin$path, quiet = TRUE)
  } else {
    library(omicsCore, lib.loc = c(origin$path, .libPaths()))
  }
  out <- list(locale = Sys.getlocale("LC_CTYPE"))
  look <- function(x) list(
    identical = identical(x[1:2], ids[1:2]),
    encoding = Encoding(x[1L]),
    # What Shiny would send to the browser
    json = as.character(jsonlite::toJSON(x[1L]))
  )
  from_csv <- omicsCore::read_omics(csv, omics_type = "proteomics",
                                    assay_type = "normalized_intensity")$input
  from_xlsx <- omicsCore::read_omics(xlsx, omics_type = "proteomics",
                                     assay_type = "normalized_intensity")$input
  out$csv <- look(rownames(from_csv$expr_mat))
  out$xlsx <- look(rownames(from_xlsx$expr_mat))
  out$csv_samples <- identical(colnames(from_csv$expr_mat)[1:2], samples[1:2])
  out$same_matrix <- identical(from_csv$expr_mat, from_xlsx$expr_mat)

  p <- omicsCore::omics_project(project_name, experiments = list(proteomics = from_csv))
  f <- tempfile(fileext = ".omp")
  omicsCore::save_project(p, f)
  back <- omicsCore::load_project(f)
  out$omp_name <- identical(back$name, project_name)
  out$omp_ids <- identical(rownames(back$experiments$proteomics$expr_mat), ids)

  from_csv$meta_df$group <- rep(c("A", "B"), each = 3L)
  b <- omicsCore::run_diff(from_csv, method = "ttest", analysis_type = "group",
                           group_col = "group", control_group = "A",
                           case_group = "B")
  out$diff_ran <- is.data.frame(b$results$diff_result_df) &&
    identical(b$results$diff_result_df$feature_id[1:2], ids[1:2])
  script <- omicsCore::export_script(p)
  out$script_parses <- !inherits(try(parse(text = script), silent = TRUE), "try-error")
  out
}

test_that("CJK identifiers survive import, save, load and analysis under LC_ALL=C", {
  skip_on_cran()
  skip_if_not_installed("callr")
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("readxl")
  skip_if_not_installed("jsonlite")
  files <- write_cjk_files()
  on.exit(unlink(unlist(files)), add = TRUE)

  got <- in_child_process(c_locale_scenario, files$csv, files$xlsx,
                          cjk_ids, cjk_samples, cjk_project,
                          env = c(LANG = "C", LC_ALL = "C", LC_CTYPE = "C"))

  expect_identical(got$locale, "C")
  # Both readers hand back the same, correctly marked names ...
  expect_true(got$csv$identical)
  expect_true(got$xlsx$identical)
  expect_identical(got$csv$encoding, "UTF-8")
  expect_true(got$csv_samples)
  expect_true(got$same_matrix)
  # ... and what reaches the browser is the name, not its bytes
  expected_json <- sprintf('["%s"]', cjk_ids[1L])
  expect_identical(charToRaw(got$csv$json), charToRaw(expected_json))
  expect_identical(charToRaw(got$xlsx$json), charToRaw(expected_json))
  expect_false(grepl("<e5>", got$csv$json, fixed = TRUE))

  expect_true(got$omp_name)
  expect_true(got$omp_ids)
  expect_true(got$diff_ran)
  expect_true(got$script_parses)
})
