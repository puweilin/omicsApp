# Run, from the source tree, every test file that reads the repository
# rather than the installed package, and fail if any of them skipped for
# want of it.
#
#   Rscript .github/scripts/run-source-tests.R
#
# Some tests check the package *sources* (no set.seed() outside one
# helper, every UI input wired to the server, the stylesheet's accents)
# or the deploy/ directory beside them. Under R CMD check the tests run
# from a copy in <pkg>.Rcheck/, where neither is next to them, so they
# skip -- correctly for a CRAN-like install, but on CI that meant they
# never ran at all. The deploy contract and the two hygiene files also
# find the checkout through OMICSAPP_REPO_ROOT under check; the rest are
# run here, with the package loaded from source the way devtools::test()
# does.
#
# The files are found rather than listed: any test file whose skip
# message says the sources or deploy/ are missing is one of them, so a
# new test of that kind is picked up without editing this script. Such a
# skip here is a failure -- it would mean the job is not doing its job.

if (!requireNamespace("testthat", quietly = TRUE)) stop("testthat is required")

repo <- normalizePath(Sys.getenv("OMICSAPP_REPO_ROOT", "."), winslash = "/")
missing_source <- paste(c(
  "package source", "package sources", "sources not beside", "installed package has no sources",
  "deploy/ directory not found", "styles\\.scss not found", "workflows not present"
), collapse = "|")

failed <- FALSE
for (pkg in c("omicsCore", "omicsApp")) {
  test_dir <- file.path(repo, "packages", pkg, "tests", "testthat")
  files <- list.files(test_dir, pattern = "^test-.*\\.R$")
  reads_repo <- vapply(files, function(f) {
    any(grepl(missing_source, readLines(file.path(test_dir, f), warn = FALSE)))
  }, logical(1))
  files <- files[reads_repo]
  if (!length(files)) next
  cat(sprintf("\n== %s: %s\n", pkg, paste(files, collapse = ", ")))

  filter <- paste0("^(", paste(sub("^test-(.*)\\.R$", "\\1", files), collapse = "|"), ")$")
  res <- testthat::test_local(file.path(repo, "packages", pkg), filter = filter,
                              reporter = "summary", stop_on_failure = FALSE)
  df <- as.data.frame(res)
  if (any(df$failed > 0) || any(df$error)) failed <- TRUE

  # Skips for the reason this job exists are failures; other skips (a
  # tool missing on the runner) are reported and allowed.
  for (r in res) {
    for (e in r$results) {
      if (inherits(e, "expectation_skip")) {
        msg <- conditionMessage(e)
        repo_skip <- grepl(missing_source, msg)
        cat(sprintf("%s %s / %s: %s\n", if (repo_skip) "FAIL skip" else "skip", r$file, r$test, msg))
        if (repo_skip) failed <- TRUE
      }
    }
  }
}

if (failed) quit(status = 1)
