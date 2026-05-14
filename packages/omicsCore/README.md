# omicsCore

Headless multi-omics analysis engine for proteomics and transcriptomics.

This is the analysis layer of the [omicsApp](https://github.com/puweilin/omicsApp) project.
For the interactive Shiny interface, see the [`omicsApp`](../omicsApp) package.

## Install

```r
# From source (development)
devtools::install_local("path/to/omicsCore")

# Optional: install heavy backends on demand
omicsCore::install_optional("rnaseq")     # DESeq2, edgeR
omicsCore::install_optional("enrichment") # clusterProfiler, msigdbr, fgsea, GSVA
omicsCore::install_optional("all")
```

## Quickstart

See [docs/export-manifest.md](../../docs/export-manifest.md) for the planned API.
Functions will be added incrementally during Phase 1.

## Status

Phase 0 — skeleton. No functions implemented yet.
