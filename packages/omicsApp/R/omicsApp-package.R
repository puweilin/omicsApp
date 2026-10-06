#' omicsApp: Shiny Web Interface for Multi-Omics Analysis
#'
#' Interactive Shiny application for proteomics and transcriptomics analysis,
#' powered by the [omicsCore](https://github.com/puweilin/omicsApp) engine.
#'
#' The app walks a project through the same steps as the engine, one view
#' each: **Import** (Excel/CSV/TSV workbooks, separate sample sheets,
#' orientation and assay checks), **Quality control** (missing values,
#' depth, outlier samples, imputation, excluding a sample), **Differential**
#' (limma, edgeR, DESeq2, t-test or linear models; several comparisons,
#' paired designs, continuous variables, a global test), **Enrichment**
#' (ORA and GSEA against MSigDB collections, across comparisons),
#' **Integration** (two layers by concordance, correlation or
#' ActivePathways) and **Report** (an HTML/PDF report and an R script that
#' reproduces the analysis, with its data files). Projects are saved to a
#' per-user store and autosaved per browser session.
#'
#' Start it with [launch()]; [shiny_app()] returns the application object
#' for hosting it yourself.
#'
#' @keywords internal
"_PACKAGE"
