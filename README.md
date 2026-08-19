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
2. **Fast interactive UI** — an analysis runs once and every view reads the resulting `analysis_bundle`, so moving a significance threshold re-masks the existing result rather than re-running it (the threshold sliders are debounced, so a drag costs one re-mask instead of one per pixel). Long-running steps are offloaded to a background `future` worker, and MSigDB gene sets load from an on-disk cache — ~11s to ~0.2s in a cold process.
3. **Dual-omics native** — a project can hold both proteomics and RNA-seq experiments side by side, with integration views (RNA-Protein correlation, joint pathway enrichment).
4. **Smart input parsing** — auto-detects expression matrix, sample metadata, and feature annotation across Excel sheets, CSV/TSV, MaxQuant, and Salmon outputs. User confirms inferred schema before running.
5. **Runs without privileges** — the packages need no Docker and no system installs; `devtools::install_local()` into a user library is the supported path on a locked-down machine. Where root *is* available, [`deploy/`](./deploy) additionally serves the app to a small group over the LAN, one container per logged-in user.

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
