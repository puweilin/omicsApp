# Steps of a long computation, said as they start.
#
# A DESeq2 fit on 60,000 genes takes minutes, and the app could only show
# a bar creeping towards an end it did not know; the engines report no
# fraction. What they can say is which step they are on, and "Estimating
# dispersions (2 of 4)" is the difference between a slow run and one
# that looks hung.
#
# The engine does not know who is listening. A caller that wants the
# steps sets `options(omicsCore.progress = function(step) ...)`: the
# Shiny app sets one that updates its notification (in the main process)
# or writes to a file it polls (in a worker). Unset, this is a no-op, so
# scripts and tests see nothing.

report_progress <- function(step) {
  hook <- getOption("omicsCore.progress")
  if (is.function(hook)) tryCatch(hook(step), error = function(e) NULL)
  invisible(NULL)
}

# DESeq() fits in one call and says what it is doing only as messages.
# Those are mapped onto the steps above and then muffled, which is what
# `quiet = TRUE` did: the fit and its results are the same either way.
DESEQ2_STEPS <- c(
  "estimating size factors"  = "Estimating size factors (1 of 4)",
  "using pre-existing"       = "Estimating size factors (1 of 4)",
  "estimating dispersions"   = "Estimating dispersions (2 of 4)",
  "fitting model and testing" = "Fitting the model and testing (3 of 4)",
  "replacing outliers"       = "Refitting genes with outlier counts (4 of 4)"
)

with_deseq2_progress <- function(expr) {
  withCallingHandlers(expr, message = function(m) {
    msg <- tolower(conditionMessage(m))
    hit <- names(DESEQ2_STEPS)[vapply(names(DESEQ2_STEPS), grepl, logical(1),
                                      x = msg, fixed = TRUE)]
    if (length(hit)) report_progress(DESEQ2_STEPS[[hit[1L]]])
    invokeRestart("muffleMessage")
  })
}

# BiocParallel back-end for DESeq2, or NULL for the serial fit.
#
# Off unless asked for, with `options(omicsCore.deseq2_workers = n)`:
# on the production host the app already runs two background workers on
# two cores, and a third set of forks there took a 60k x 300 fit from
# 19.8 s to 15.0 s while slowing every other session. A machine with
# cores to spare can turn it on. Forked workers on Unix, sockets on
# Windows (which cannot fork).
deseq2_bpparam <- function() {
  n <- suppressWarnings(as.integer(getOption("omicsCore.deseq2_workers", 1L)))
  if (is.na(n) || n < 2L || !is_installed("BiocParallel")) {
    return(NULL)
  }
  if (.Platform$OS.type == "windows") {
    BiocParallel::SnowParam(workers = n)
  } else {
    BiocParallel::MulticoreParam(workers = n)
  }
}
