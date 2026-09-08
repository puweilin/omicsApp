# Where a differential's time actually goes, on this machine and this data.
#
#   Rscript time_pipeline.R /data/_autosave.omp [layer]
#
# Four segments, because they have four different fixes and guessing
# which one is large has been wrong twice:
#
#   1. fit          DESeq2's per-gene GLMs. Fixed by testing fewer genes
#                   (the QC low-expression filter) or by more CPU.
#   2. extract      results() plus omicsCore's standardisation, symbol
#                   mapping and schema. Fixed in omicsCore.
#   3. return       the bundle crossing back from the worker process.
#                   Fixed by putting less in the bundle.
#   4. autosave     qs_save of the whole project. Same fix as 3, and it
#                   happens on every change rather than once per run.
#
# Read as a whole: if 1 dominates, the analysis is the cost and that is
# arithmetic nobody can avoid. If 3 and 4 dominate, something is being
# carried that does not need to be.

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

# ---- 1 + 2: fit, and everything after it -----------------------------
#
# run_diff() is fit + extract together. The bare DESeq2 call underneath
# is the fit alone, so the difference is what omicsCore adds.
whole <- secs(run_diff(input, method = method, analysis_type = "group",
                       group_col = group_col, control_group = control,
                       case_group = case))
bundle <- whole$value

fit_secs <- NA_real_
if (identical(method, "deseq2") && requireNamespace("DESeq2", quietly = TRUE)) {
  keep <- rownames(input$meta_df)[
    input$meta_df[[group_col]] %in% c(control, case)]
  cd <- data.frame(g = factor(input$meta_df[keep, group_col],
                              levels = c(control, case)),
                   row.names = keep)
  dds <- DESeq2::DESeqDataSetFromMatrix(
    round(input$expr_mat[, keep, drop = FALSE]), cd, design = ~ g)
  fit <- secs(suppressMessages(DESeq2::DESeq(dds, quiet = TRUE)))
  fit_secs <- fit$secs
}

say("1 fit", fit_secs, "DESeq2's per-gene GLMs")
say("2 extract", whole$secs - (if (is.na(fit_secs)) 0 else fit_secs),
    "results() + standardisation")

# ---- 3: back from the worker -----------------------------------------
#
# Serialisation is what crossing the boundary costs; the fork itself is
# copy-on-write and nearly free. Measured on the bundle as the app sends
# it, which is not the same object as the one above once run_diff has
# stopped keeping the model.
ser <- secs(serialize(bundle, NULL))
say("3 return", ser$secs,
    sprintf("%.0f MB serialised", length(ser$value) / 1024^2))

# ---- 4: autosave ------------------------------------------------------
#
# The whole project, not just this result, and on every change rather
# than once per run.
project$bundles$diff <- bundle
tmp <- tempfile(fileext = ".omp")
on.exit(unlink(tmp), add = TRUE)
sav <- secs(save_project(project, tmp))
say("4 autosave", sav$secs,
    sprintf("%.0f MB on disk", file.size(tmp) / 1024^2))

cat("\n")
