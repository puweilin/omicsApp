# omicsApp

[![R-CMD-check](https://github.com/puweilin/omicsApp/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/puweilin/omicsApp/actions/workflows/R-CMD-check.yaml)

A multi-omics analysis platform for proteomics and transcriptomics data,
providing a unified analysis engine plus a modern Shiny web interface.

## Packages

This monorepo contains two R packages:

| Package | Purpose | Status |
|---|---|---|
| [`omicsCore`](./packages/omicsCore) | Headless analysis engine: QC, differential expression, enrichment (ORA/GSEA/GSVA), and multi-omics integration | Phase 0 — skeleton |
| [`omicsApp`](./packages/omicsApp) | Shiny web interface built on `omicsCore`: smart input parsing, interactive plots, project sessions | Phase 0 — skeleton |

## Design goals

1. **Lightweight install** — base `omicsCore` has only ~15 CRAN dependencies. Heavy Bioconductor packages (DESeq2, limma, clusterProfiler, GSVA, ComplexHeatmap) are `Suggests` and prompted on first use.
2. **Fast interactive UI** — analysis runs once, all views read from a cached `analysis_bundle`. `bindCache` + `future` + `debounce` throughout.
3. **Dual-omics native** — a project can hold both proteomics and RNA-seq experiments side by side, with integration views (RNA-Protein correlation, joint pathway enrichment).
4. **Smart input parsing** — auto-detects expression matrix, sample metadata, and feature annotation across Excel sheets, CSV/TSV, MaxQuant, and Salmon outputs. User confirms inferred schema before running.
5. **Runs in restricted environments** — no Docker, no system installs. `devtools::install_local()` to user library is the supported path.

## Development workflow

```r
# Load without installing (for development)
devtools::load_all("packages/omicsCore")
devtools::load_all("packages/omicsApp")

# Install to user library
devtools::install_local("packages/omicsCore")
devtools::install_local("packages/omicsApp")

# Launch the Shiny app
omicsApp::launch()
```

## Status

Pre-alpha. See [docs/export-manifest.md](./docs/export-manifest.md) for the planned public API and [docs/roadmap.md](./docs/roadmap.md) for the delivery plan.

The legacy `omics_core` framework at `CHISSS/scripts/frameworks/omics_core/` remains **frozen** as the production code path for existing CHISSS analyses. `omicsCore` is a parallel, properly-packaged rewrite that does not affect those scripts.

## License

MIT (c) 2026 Weilin Pu. See [LICENSE](./LICENSE).
