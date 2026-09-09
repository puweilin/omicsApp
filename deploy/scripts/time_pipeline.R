# Where a differential's time actually goes, on this machine and this data.
#
#   Rscript time_pipeline.R /data/_autosave.omp [layer]
#
# Three segments, because they have three different fixes and guessing
# which one is large has been wrong twice:
#
#   1. analysis     run_diff() end to end. Fixed by testing fewer genes
#                   (the QC low-expression filter) or by more CPU.
#   2. return       the bundle crossing back from the worker process.
#                   Fixed by putting less in the bundle.
#   3. autosave     qs_save of the whole project. Same fix as 2, and it
#                   happens on every change rather than once per run.
#
# Read as a whole: if 1 dominates, the analysis is the cost and that is
# arithmetic nobody can avoid -- profile it if you need the split. If 2
# and 3 dominate, something is being carried that does not need to be.

suppressPackageStartupMessages({
  library(omicsCore)
})

# Defined here rather than taken from omicsCore, which does not export it,
# or from base, which only has it since R 4.4. This script is meant to be
# runnable against an older image than the one it shipped with.
`%||%` <- function(a, b) if (is.null(a)) b else a

args <- commandArgs(trailingOnly = TRUE)
# Not `args[[1L]] %||% ...`: an empty vector is not NULL, and [[1]] on one
# is a subscript error rather than a missing default.
path <- if (length(args) >= 1L) args[[1L]] else "/data/_autosave.omp"
if (!file.exists(path)) {
  stop("No project at ", path,
       ". Pass one: Rscript time_pipeline.R /data/<file>.omp")
}

secs <- function(expr) {
  t0 <- Sys.time()
  value <- force(expr)
  list(value = value,
       secs = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}
say <- function(label, s, note = "") {
  cat(sprintf("  %-10s %7.1f s   %s\n", label, s, note))
}

project <- load_project(path)
tags <- experiment_tags(project)
tag <- if (length(args) >= 2L) args[[2L]] else {
  types <- vapply(project$experiments, function(e) e$omics_type %||% "",
                  character(1))
  rna <- names(types)[types == "rnaseq"]
  if (length(rna)) rna[[1L]] else tags[[1L]]
}
if (!tag %in% tags) {
  stop("No layer '", tag, "'. This project has: ", paste(tags, collapse = ", "))
}

# What the app would analyse, not what was imported: a QC apply replaces
# it, and timing the wrong one reports the work of a different run.
eff <- effective_input(project, tag)
input <- eff$input
cat(sprintf("\n%s  layer=%s  %s features x %s samples  assay=%s\n\n",
            project$name %||% "(unnamed)", tag,
            format(nrow(input$expr_mat), big.mark = ","),
            format(ncol(input$expr_mat), big.mark = ","),
            input$assay_type %||% "?"))

# The contrast the project already used, so this times the run they did.
saved <- project$bundles$diff$params
group_col <- saved$group_col %||% NULL
if (is.null(group_col)) {
  cand <- names(input$meta_df)[vapply(input$meta_df, function(x) {
    n <- length(unique(stats::na.omit(x)))
    n >= 2L && n <= 6L
  }, logical(1))]
  if (length(cand) == 0L) stop("No grouping column found; pass one by hand.")
  group_col <- cand[[1L]]
}
lv <- unique(as.character(stats::na.omit(input$meta_df[[group_col]])))
control <- saved$control_group %||% lv[[1L]]
case <- saved$case_group %||% lv[[2L]]
method <- saved$method %||% "deseq2"
cat(sprintf("contrast: %s  %s vs %s   method=%s\n\n",
            group_col, case, control, method))

# ---- 1: the analysis, in one real call --------------------------------
#
# One call, timed as a whole. An earlier version split it into fit and
# extract by running DESeq2 a second time and subtracting -- which
# measures neither: two fits do not cost the same (BLAS caches, memory
# pressure, the dispersion fitType falling back on one run and not the
# other), and the difference of two noisy numbers is noisier than either.
#
# Instrumenting the boundary from outside was tried and did not take:
# neither trace() nor a namespace rebind reliably reached the call from
# inside a sourced script, and the failure printed `1 fit 0.0 s` -- a
# failed measurement wearing the costume of a fast one. A diagnostic must
# not do that, so the split is simply not claimed.
#
# What the four numbers below are for is the coarser question, which is
# the one that has been wrong twice: is the analysis the cost, or is the
# plumbing? For that, "run_diff took N seconds" against segments 2 and 3
# is enough. If run_diff dominates and you then need to know *where* in
# it, profile it: Rprof() on this call reports the split without anything
# having to be patched.
run <- secs(run_diff(input, method = method, analysis_type = "group",
                     group_col = group_col, control_group = control,
                     case_group = case,
                     # The model the project actually used. Omitting these
                     # timed a simpler analysis than the one that was slow.
                     covariates = saved$covariates,
                     paired_col = saved$paired_col))
bundle <- run$value
say("1 analysis", run$secs,
    sprintf("run_diff(%s), the whole call", method))

# ---- 3: back from the worker -----------------------------------------
#
# Serialisation, which is the dominant cost of crossing the boundary --
# the fork itself is copy-on-write and nearly free. It is not the same
# number as a real worker round trip: that also pays for scheduling and
# for the parent reading the result back, and it depends on load. Read it
# as "how much has to be moved", not as "how long the worker took to
# hand it over".
ser <- secs(serialize(bundle, NULL))
say("2 return", ser$secs,
    sprintf("%.0f MB serialised", length(ser$value) / 1024^2))

# ---- 4: autosave ------------------------------------------------------
#
# The whole project, not just this result, and on every change rather
# than once per run.
project$bundles$diff <- bundle
tmp <- tempfile(fileext = ".omp")
on.exit(unlink(tmp), add = TRUE)
sav <- secs(save_project(project, tmp))
say("3 autosave", sav$secs,
    sprintf("%.0f MB on disk", file.size(tmp) / 1024^2))

cat("\n")
