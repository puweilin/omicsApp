# omicsCore

Headless multi-omics analysis engine for proteomics and transcriptomics.

This is the analysis layer of the [omicsApp](https://github.com/puweilin/omicsApp) project.
For the interactive Shiny interface, see the [`omicsApp`](../omicsApp) package.

## Install

```r
# From source
devtools::install_local("path/to/omicsCore")

# Optional: install heavy backends on demand
omicsCore::install_optional("rnaseq")     # DESeq2, edgeR
omicsCore::install_optional("enrichment") # clusterProfiler, msigdbr, fgsea, GSVA
omicsCore::install_optional("all")
```

## Quickstart

```r
library(omicsCore)
inp  <- read_omics("proteomics.xlsx", omics_type = "proteomics")
qc   <- run_qc(inp)
diff <- run_diff(qc$results$cleaned_input, method = "limma",
                 group_col = "group", control_group = "Control",
                 case_group = c("TreatA", "TreatB"))
plot_volcano(select_comparison(diff, "TreatA_vs_Control"))
enr  <- run_enrichment(select_comparison(diff, "TreatA_vs_Control"),
                       type = "ora", database = "hallmark")

proj <- omics_project("my study", list(proteomics = inp))
export_report(proj, "report.html")
export_script(proj, "analysis.R")
```

Every analysis returns an `analysis_bundle` holding its parameters,
results and notes, so a figure, a report and an exported script are all
read from the same object.

## Status

0.2.0 — see [NEWS.md](./NEWS.md).
