# Roadmap

| Phase | Goal | Estimated effort |
|---|---|---|
| **0** | Repo skeleton, API contract, CI placeholder | 1 week |
| **1** | `omicsCore` packaged: port existing framework, roxygen docs, tests pass, `R CMD check` clean | 2-3 weeks |
| **2** | `omicsApp` UI skeleton: bslib theme, modular structure, mock data flows | 2 weeks |
| **3** | Wire UI to `omicsCore`: import wizard, QC, diff, enrichment, integration views | 2-3 weeks |
| **4** | Polish: perf benchmarks, error states, end-to-end tests, pkgdown sites | 1-2 weeks |
| **5** (optional) | Olink / hormone / ISF modules; collaborator deployment | TBD |

**Total to v0.1.0 (usable internal release): 8-11 weeks.**

## Phase 0 — current

- [x] Decide architecture (monorepo, 2 packages, frozen legacy `omics_core`)
- [x] Create directory skeleton
- [x] DESCRIPTION + NAMESPACE for both packages
- [x] LICENSE + .Rbuildignore + .gitignore
- [x] [export-manifest.md](./export-manifest.md) — API contract
- [x] Placeholder Shiny app
- [ ] Push to GitHub
- [ ] Add GitHub Actions R-CMD-check workflow (deferred to Phase 1)
- [ ] Set up `usethis::use_dev_package()` cross-references (deferred to Phase 1)

## Phase 1 — omicsCore packaging

Migration plan from `scripts/frameworks/omics_core/R/`:

| Source (frozen) | Destination | Strategy |
|---|---|---|
| `R/core/*.R` | `packages/omicsCore/R/` (flat) | Copy, add roxygen, `@export` selected functions |
| `R/data_input/*.R` | `packages/omicsCore/R/io_*.R` | Refactor R6 managers into S3 constructors |
| `R/qc/*.R` | `packages/omicsCore/R/qc_*.R` | Copy, document |
| `R/diff/*.R` | `packages/omicsCore/R/diff_*.R` | Copy, document |
| `R/enrich/*.R` | `packages/omicsCore/R/enrich_*.R` | Copy, document |
| `R/plot/*.R` | `packages/omicsCore/R/plot_*.R` | Copy, add `interactive` param via plotly |
| `R/pipelines/*.R` | `packages/omicsCore/R/pipeline_*.R` | Copy, simplify |
| `R/workflows/*.R` | `packages/omicsCore/R/workflow_*.R` | Copy, simplify |
| `R/legacy/*.R` | **do not copy** | Remain in CHISSS legacy frameworks |
| `R/bootstrap/*.R` | **do not copy** (replaced by NAMESPACE) | |

New additions in Phase 1:

- `R/class-omics_project.R` — multi-omics container
- `R/io_smart_read.R` + `R/io_inference.R` — smart input parser
- `R/integrate_*.R` — dual-omics integration (port from `scripts/analyses/integration/`)
- `R/install_optional.R` — Suggests bootstrapper

## Exit criteria per phase

**Phase 1 done when:** `R CMD check --as-cran packages/omicsCore` returns 0 errors / 0 warnings; all current CHISSS analyses can be reproduced via `library(omicsCore)` without sourcing the legacy bootstrap.

> **Status: met.** Validation script at `CHISSS/scripts/validation/cheek_g2_vs_g1_omicscore.R` reproduces the legacy Cheek G2-vs-G1 limma-with-age-adjustment diff (Pearson r = 1.000 on log2FC and -log10 p-values; volcano + heatmap + export_bundle all round-trip).
>
> **Correction.** That validation started from `readRDS(legacy_rds)` → `manager$get_imputed()`, i.e. from data the legacy framework had already vsn-normalized and imputed. It proved the differential layer faithful and said nothing about the preprocessing layer — which had not been ported, because `R/data_input/*` (migration table above) became the import *wizard* (`read_omics()`, sheet classification) rather than the data *preparation* pipeline the R6 managers ran. Normalization lived in `data_input/proteomics/`, not in `qc/`, so porting "the QC layer" faithfully missed it.
>
> The gap is closed by `normalize_omics()` plus the `assay_type` vocabulary, and re-validated by `packages/omicsApp/tests/testthat/test-golden-raw-to-result.R`, which starts at the raw workbook so that nothing which only holds for pre-processed input can pass.

**Phase 2 done when:** `omicsApp::launch()` opens a styled multi-page UI with all 7 views navigable, populated by built-in example data, no business logic yet.

**Phase 3 done when:** A new user can: upload an Excel file → smart-parse → **declare the value scale and normalize** → QC → run diff → see volcano + enrichment + integration → download a report — without writing any R code.

> The normalization step is named here on purpose. The earlier wording went straight from smart-parse to QC, so nothing on the checklist ever asked where linear intensities become analysis-ready — and limma applies no transform of its own, so skipping it produces a full, plausible, wrong result table rather than an error.

**Phase 4 done when:** UI is < 100 ms input latency on a 10k × 100 dataset; all `shinytest2` end-to-end tests pass; pkgdown sites deployed.
