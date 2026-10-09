# omicsApp

[![R-CMD-check](https://github.com/puweilin/omicsApp/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/puweilin/omicsApp/actions/workflows/R-CMD-check.yaml)

A multi-omics analysis platform for proteomics and transcriptomics data,
providing a unified analysis engine plus a modern Shiny web interface.

## Packages

This monorepo contains two R packages:

| Package | Purpose | Status |
|---|---|---|
| [`omicsCore`](./packages/omicsCore) | Headless analysis engine: QC, differential expression, enrichment (ORA/GSEA/GSVA), and multi-omics integration | 0.2.0 |
| [`omicsApp`](./packages/omicsApp) | Shiny web interface built on `omicsCore`: smart input parsing, interactive plots, project sessions | 0.2.0 |

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

### Tests

```sh
NOT_CRAN=true Rscript -e 'devtools::test("packages/omicsCore")'
NOT_CRAN=true Rscript -e 'devtools::test("packages/omicsApp")'
```

`NOT_CRAN=true` matters: without it the golden test in
`packages/omicsApp/tests/testthat/test-golden-raw-to-result.R` skips silently.
That test runs the real SkinProteomics workbook from the raw file through
import and normalization and checks the result against the legacy
pipeline. It looks for `data/SkinProteomics/Proteomics_Data.xlsx` by walking
up from the working directory, or wherever `OMICSAPP_PROTEOMICS_XLSX` points; it
skips when the file is absent, so a clone without the data still runs green.

Which omicsCore the omicsApp suite tests depends on how omicsApp itself
is loaded (`packages/omicsApp/tests/testthat/setup.R`):

- **From source** (`devtools::test()`, `testthat::test_local()`,
  `load_all()`): the omicsCore *source tree* beside it, loaded with
  pkgload, not the installed copy, so a development run does not depend
  on when `install_local()` was last run. The shinytest2 tests boot
  `inst/app` from source with both packages loaded the same way
  (`OMICSAPP_DEV_ROOT`, see `inst/app/app.R`).
- **Installed** (`R CMD check`, `testthat::test_check()` or
  `test_package()`): the installed packages only; no source is ever
  loaded, and the browser tests run the installed app. Install omicsCore
  from the same checkout first, as CI does (`R CMD INSTALL
  packages/omicsCore` in `.github/workflows/R-CMD-check.yaml`).

`OMICSAPP_TEST_CORE=installed` or `OMICSAPP_TEST_CORE=source` forces one or
the other, e.g. to run the development suite against an installed
omicsCore:

```sh
OMICSAPP_TEST_CORE=installed NOT_CRAN=true Rscript -e 'devtools::test("packages/omicsApp")'
```

A few suites are deliberately opt-in:

| Suite | Runs when | What it holds |
|---|---|---|
| `omicsCore/tests/testthat/test-perf-budget.R` | `OMICSCORE_PERF_TESTS=1` | import, PCA and limma at the size of the real follicle dataset (63k x 258) |
| `omicsCore/tests/testthat/test-concurrent-writers.R` | `callr` installed | two R processes autosaving one `.omp` while a third reads it |
| `omicsApp/tests/testthat/test-app-smoke.R` | Chrome available | the seven views in a browser, from source |
| `omicsApp/tests/testthat/test-deploy-contract.R` | `deploy/` present | the ports, paths and package lists that three deploy files must agree on |
| `omicsApp/tests/testthat/test-app-journey.R` | Chrome available | upload, run, enrich, download the script, replace the data, in a browser |
| `omicsCore/tests/testthat/test-golden-rnaseq-follicle.R` | the parent project's `data/` and `results/`, or `OMICSAPP_FOLLICLE_ROOT` | the follicle counts file through filter, winsorize and edgeR, against the legacy report's table |
| `*/tests/testthat/test-c-locale.R` | `callr` installed | the package in a child process running under `LC_ALL=C`, as a container without a locale would |
| `omicsCore/tests/testthat/test-api-contract.R` (heavy sweep) | `OMICSCORE_FUZZ_TESTS=1` | every argument of the enrichment, GSVA, integration and report functions given every wrong value in turn |
| `omicsApp/tests/testthat/test-two-sessions.R` | Chrome available | two browsers on one R process; what one imports the other must not see |
| `omicsApp/tests/testthat/test-integration-correlation-browser.R` | Chrome available | a saved paired-sample project's correlation scatter in a browser: a hovered point's card, and the dated PNG download |

`omicsCore/tests/testthat/fixtures/omp/` is a corpus of `.omp` files written
by past versions; every file in it must keep opening. When the schema
version or the bundle layout changes, add a file rather than replacing one
(`write_omp_corpus_fixture()` in `test-omp-corpus.R` writes the current
version's).

Every exported omicsCore function answers a wrong argument by name
(`test-api-contract.R`): the light functions are swept on every run with
NULL, NA, a negative number, a string, an empty vector, a list, a function
and a hollow object in each argument, and an error message that is R's own
("argument is of length zero", "$ operator is invalid for atomic vectors")
fails the test. The assertions live in `R/validation-utils.R`; use them in
any new entry point rather than letting a bad value reach the engine.

`test-hygiene.R` (omicsCore) and `test-app-hygiene.R` (omicsApp) hold what a
call may not leave behind in a process that serves every user: the random
stream, a planted seed, an open connection or device, an option, a file in
`tempdir()`. Seed only through `with_fixed_seed()` (core) or `local_seed()`
(app); a bare `set.seed()` in package code fails a test.

Optional-dependency branches are tested by mocking `omicsCore:::is_installed()`
and `omicsApp:::has_pkg()` (`test-lightweight-install.R`,
`test-missing-deps.R`); any new `requireNamespace()` call outside those two
functions fails a test, because it would silently leave that branch untested
again.

## Status

Version 0.2.0 (see each package's `NEWS.md`). The reviews and the plan for
the next releases are in [`docs/`](./docs): [review-2026-10.md](./docs/review-2026-10.md)
and [review-2026-10-round4.md](./docs/review-2026-10-round4.md) record the
October 2026 audits and fixes, and [next-release-plan.md](./docs/next-release-plan.md)
lists what is done and what remains. [docs/export-manifest.md](./docs/export-manifest.md)
describes the public API.

The legacy `omics_core` framework in the parent project's `scripts/frameworks/omics_core/` remains **frozen** as the production code path for the existing analyses. `omicsCore` is a parallel, properly-packaged rewrite that does not affect those scripts.

## License

MIT (c) 2026 Weilin Pu. See [LICENSE](./LICENSE).
