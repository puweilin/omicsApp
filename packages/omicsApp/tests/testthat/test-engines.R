# The differential view offers whichever engines load, and run_diff()'s
# auto mode falls back to a t-test, with a message, when the preferred
# one does not. That fallback is right for a user and wrong for a test
# suite: on CI it turned "deseq2" into "ttest" and the diff-controls
# tests failed on the symptom without saying why DESeq2 would not load.
# This asks the question directly, with the answer in the failure.

test_that("every engine omicsCore can dispatch to loads in this process", {
  for (pkg in c("limma", "DESeq2", "edgeR")) {
    skip_if(!nzchar(system.file(package = pkg)), paste(pkg, "is not installed"))
    outcome <- tryCatch({
      loadNamespace(pkg)
      "loaded"
    }, error = function(e) conditionMessage(e))
    expect_identical(outcome, "loaded",
                     info = paste0(pkg, " from ", system.file(package = pkg),
                                   "\n.libPaths(): ", paste(.libPaths(), collapse = ", "),
                                   "\n", outcome))
    expect_true(omicsCore:::is_installed(pkg), info = pkg)
  }
})
