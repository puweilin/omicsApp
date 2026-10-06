# omicsApp

Shiny web interface for multi-omics analysis, built on [`omicsCore`](../omicsCore).

A project goes through one view per step: **Import** (Excel/CSV/TSV, a
separate sample sheet, orientation and assay checks) → **Quality control**
(missing values, depth, outlier samples, imputation) → **Differential**
(limma, edgeR, DESeq2, t-test or linear models; several comparisons, paired
designs, continuous variables, a global test) → **Enrichment** (ORA/GSEA on
MSigDB, across comparisons) → **Integration** (two layers) → **Report**
(HTML/PDF report and a reproducible R script with its data files). Projects
are saved per user and autosaved per browser session.

## Install

```r
devtools::install_local("path/to/omicsCore")
devtools::install_local("path/to/omicsApp")
```

## Launch

```r
omicsApp::launch()                      # opens in the browser
omicsApp::launch(host = "0.0.0.0",      # behind a proxy / in a container
                 launch.browser = FALSE)
shiny::runApp(omicsApp::shiny_app())    # the app object, to host yourself
```

The built-in tutorial (Project view → *Load tutorial*) runs the whole
workflow on a small simulated two-treatment dataset.

## Status

0.2.0 — see [NEWS.md](./NEWS.md). Multi-user deployment: [`deploy/`](../../deploy).
