# omicsDesk: the offline, single-machine edition

Status: proposal (2026-09-05). Nothing below is implemented yet.

`omicsApp` today is installed with `devtools::install_local()` on a machine
that can reach CRAN and Bioconductor, or served from `deploy/` on a LAN.
This document plans a third package, **omicsDesk**, for the remaining case:
one user, one laptop, no network at all, no admin rights.

The design principle is the same one the deploy image already follows:
*everything the app can ask for is present before the first click, and
nothing it does at runtime reaches for the network.*

## 1. Why a third package, not a fork and not a flag

omicsCore and omicsApp stay exactly what they are. omicsDesk is a thin
distribution layer on top of them:

| Concern | Lives in |
|---|---|
| analysis code, UI code, tests | omicsCore / omicsApp (unchanged roles) |
| "optional" packages become required | omicsDesk `DESCRIPTION` (promotes Suggests to Imports) |
| bundled gene sets, pandoc, launchers | omicsDesk `inst/` |
| offline defaults (cache dir, data dir, no auto-refresh) | `omicsDesk::launch()` and `.onLoad` |
| self-diagnosis, offline install, cache import | omicsDesk functions |
| building the distributable zip | `tools/build_bundle.R` + a CI workflow |

A merged single package would duplicate 100+ files and both test suites.
A `standalone = TRUE` flag on `omicsApp::launch()` cannot carry 7.5 MB of
gene sets or force Bioconductor packages to be installed. The third
package does both and keeps the existing packages honest: every offline
behaviour is reached through a seam that omicsDesk sets and a test can
set too.

## 2. What breaks offline today (audit of the working tree)

| # | Touchpoint | Where | Offline behaviour today | Fix |
|---|---|---|---|---|
| 1 | MSigDB gene sets | `omicsCore/R/enrich-genesets.R` `fetch_msigdbr_raw()` | msigdbr >= 25 downloads a 34 MB RDS from Zenodo on first call (`curl_download`, 300 s timeout); 10-24 needs `msigdbdf` from r-universe; only < 10 ships its data. With no qs2 cache, every enrichment, GSVA and `list_gene_sets()` call hangs, then fails. | Ship the prewarmed qs2 cache (14 files, 7.5 MB) in omicsDesk; teach omicsCore to read a bundled cache dir; fail fast with an actionable message when offline and uncached. |
| 2 | KEGG live refresh and 30-day TTL | `omicsCore/R/geneset-refresh.R` | `refresh_geneset_cache("kegg")` waits on `rest.kegg.jp`; the TTL only arms for live-sourced caches, so a prewarmed cache is safe, but nothing tells the user KEGG refresh is online-only. | Offline mode short-circuits both; UI labels refresh as online-only; add cache import/export so a networked colleague can refresh KEGG and hand the file over. |
| 3 | Optional-package install | `omicsCore/R/install-optional.R`; eight "Install with: install.packages(...)" hints across both packages | `install_optional()` and the hints point at CRAN/Bioconductor. | Promote every runtime-needed Suggest to an omicsDesk Import so nothing is missing; route hints through one `install_hint()` seam; `install_optional()` honours a local `file://` repo. |
| 4 | Reports | `omicsCore/R/export-report.R`, `omicsApp/R/mod_report_view.R` | Needs pandoc, which plain R does not ship (RStudio does); PDF needs LaTeX, and tinytex fetches missing LaTeX packages from the network; `html_document()` default `mathjax = "default"` injects a CDN script tag. The view gates only on rmarkdown, so a machine with rmarkdown but no pandoc shows a live button that errors. | `mathjax = NULL`; gate HTML on `rmarkdown::pandoc_available()` and PDF on a LaTeX engine, with separate notices; bundle a pandoc binary and set `RSTUDIO_PANDOC`. |
| 5 | Project storage | `omicsApp/R/project_store.R` | `OMICSAPP_DATA_DIR` defaults to `tempdir()`, so on a laptop every saved project and the autosave vanish when R exits. | `omicsApp::launch(data_dir = ...)`; omicsDesk defaults it to a persistent per-user directory. |
| 6 | Dependency closure | `DESCRIPTION` of both packages | 153 packages / 636 MB for Imports alone; 253 packages / ~1.0 GB with the Suggests the app actually uses (measured on macOS arm64, R 4.4.3, Bioc 3.20). Bioconductor packages have no CRAN binaries. | Ship a local CRAN-style repository of binaries per platform and R minor version. |

Already offline-clean, and worth keeping that way: the HGNC map is bundled
(`inst/extdata/hgnc_ensembl.rds`); import templates are generated in code;
the bslib theme names fonts without `font_google()`; `styles.scss` has no
`@import url()`; shiny, plotly, DT and bslib serve their JavaScript from
the packages. On the enrichment path, `clusterProfiler::enricher()`,
`GSEA()` and `enrichplot::gseaplot2()` make no network calls (the
download helpers in `yulab.utils` are reached only from
`clusterProfiler:::stringdb_version`, which nothing here calls).
`future::multisession` workers inherit `.libPaths()`, so a private
library works for background analyses.

## 3. What the user receives

One zip per platform and R minor version, for example
`omicsDesk-0.1.0-R4.5-windows.zip`:

```
omicsDesk-0.1.0-R4.5-windows/
  README.md, README.zh.md      install and first run, in both languages
  INSTALL.R                    installs everything from repo/ into a private library
  START.R, omicsDesk.bat       double-click launchers (omicsDesk.command on macOS)
  repo/                        local CRAN-style repository: PACKAGES index,
    src/contrib/               omicsCore, omicsApp, omicsDesk source tarballs
    bin/windows/contrib/4.5/   ~250 binary packages (CRAN + Bioconductor)
  pandoc/                      pandoc binary for HTML reports
  bundle.json                  R version, platform, Bioc release, package versions,
                               gene-set vintage, git sha
  SHA256SUMS
```

The gene sets travel inside the omicsDesk package (`inst/genesets/`), so
`library(omicsDesk)` alone is enough for enrichment.

Measured sizes on this machine:

| Item | Size |
|---|---|
| Imports closure, installed | 153 packages, 636 MB |
| Imports + promoted Suggests, installed | 253 packages, ~1.0 GB |
| gene-set cache, Hs + Mm, 7 databases | 14 files, 7.5 MB |
| msigdbr's own data file (not needed when the cache ships) | 34 MB |

Expect a zip of roughly half a gigabyte per platform; binaries compress
better than the installed footprint suggests.

## 4. Plan

### Phase 1: make omicsCore and omicsApp offline-safe (2-3 days)

Small, seam-shaped changes in the existing packages. Each one is
independently useful on a networked machine too.

1. **Layered gene-set cache.** Replace `geneset_cache_dir()` with a read
   list and a single write dir. Read order: `OMICSCORE_GENESET_CACHE`,
   then the per-user cache, then `getOption("omicsCore.geneset_cache")`
   (omicsDesk points it at its bundled directory). Writes go only to the
   user cache, never into an installed package. `read_geneset_cache()`
   walks the list; `geneset_cache_status()` reports which directory served
   each table.
2. **Offline mode.** `offline_mode()` reads `getOption("omicsCore.offline")`
   or `OMICSCORE_OFFLINE=1`. When on: `fetch_msigdbr_table()` stops with
   "no cache for hallmark / Homo sapiens; import one with
   import_geneset_cache() or refresh on a networked machine" instead of
   calling msigdbr; `maybe_refresh_kegg_cache()` returns the cached table
   untouched; `refresh_geneset_cache()` errors immediately for `"kegg"`.
   Without the guard an offline laptop sits through msigdbr's 300 s
   download timeout before failing.
3. **Cache export and import.** `export_geneset_cache(path)` zips the
   user cache with its provenance; `import_geneset_cache(path)` validates
   the column contract (`GENESET_CACHE_COLUMNS` plus `gs_source`) and
   copies into the user cache. This is how a KEGG refresh reaches an
   offline machine.
4. **Local repository for installs.** `install_optional()` honours
   `getOption("omicsCore.local_repo")` (a `file://` URL) for both the pak
   and the fallback path; in offline mode with no local repo it stops
   with a message that names the bundle.
5. **One install hint.** `install_hint(pkg)` in omicsCore, used by both
   packages, replaces the eight hard-coded sentences. Default wording is
   today's; offline wording says the package is part of the bundle and
   to run `omicsDesk::doctor()`.
6. **Report capabilities.** `report_capabilities()` returns
   `list(rmarkdown, pandoc, latex)` using `rmarkdown::pandoc_available()`
   and `tinytex::is_tinytex() || nzchar(Sys.which("pdflatex"))`. The
   Report view gates the HTML button on pandoc and the PDF button on a
   LaTeX engine, each with its own notice. `export_report()` passes
   `mathjax = NULL`.
7. **Persistent data dir.** `omicsApp::launch(data_dir = NULL)`: when
   given, sets `OMICSAPP_DATA_DIR` for the session and restores it on
   exit. The default stays `tempdir()`, which is what the server and the
   tests rely on.
8. **Offline contract test**, in both suites: `withr::local_envvar` with
   `http_proxy` and `https_proxy` pointing at `http://127.0.0.1:9`,
   `options(repos = c(CRAN = "file:///nonexistent"))`, offline mode on, a
   small prewarmed cache fixture, and `kegg_rest_table` mocked to fail.
   Then run import, QC, limma, ORA, GSEA, GSVA and the HTML report
   (skipped without pandoc). This is the test that proves the path; the
   bundle CI in phase 3 proves the install.

### Phase 2: `packages/omicsDesk` (3-4 days)

```
packages/omicsDesk/
  DESCRIPTION        Imports: omicsApp, omicsCore, and every runtime Suggest promoted
  R/zzz.R            .onLoad: geneset_cache option -> inst/genesets, offline mode on
                     unless OMICSDESK_ONLINE=1, local_repo if a bundle is registered
  R/launch.R         launch(data_dir = NULL, port = NULL, workers = 2, browser = TRUE)
  R/doctor.R         doctor(): one row per check, what is missing, how to fix it
  R/install.R        install_offline(bundle_dir): the same code as the bundle's INSTALL.R
  R/genesets.R       geneset_status(), import_geneset_cache() re-exported with desk paths
  R/paths.R          desk_data_dir(), desk_cache_dir()
  inst/genesets/     14 qs2 files + MANIFEST.json (MSigDB version, msigdbr version,
                     build date, licence notes)
  inst/bundle/       INSTALL.R, START.R, launchers, README templates
  tests/testthat/    doctor on a healthy library; launch sets and restores env;
                     import rejects a wrong-schema file; the app journey under
                     the offline contract
```

Promoted Imports, from what the app can actually call: DESeq2, edgeR,
GSVA, fgsea, enrichplot, ComplexHeatmap, circlize, ggrepel, vsn,
pcaMethods, impute, ActivePathways, rmarkdown, knitr, writexl, plus
whichever imputation backend `OPTIONAL_GROUPS$imputation` names at
release time. Deliberately not promoted, as the Dockerfile already
argues: ggpubr (unused, drags in a compiler-sensitive chain), tximport
and GenomicFeatures (no Salmon import path yet).

`launch()` does, in order: set `OMICSAPP_DATA_DIR` to `desk_data_dir()`
unless the caller gave one; set `OMICSCORE_GENESET_TTL_DAYS=0`; set
`RSTUDIO_PANDOC` when a bundled pandoc directory exists; call
`omicsApp::launch()`; restore every variable on exit.

`doctor()` checks: R version against `bundle.json`, every promoted
package loads, the 14 cache files read, pandoc found, LaTeX found (as
information, not failure), data dir writable, a `multisession` worker
starts and sees the same library. It prints a table and returns it
invisibly, so the bundle CI can assert on it.

The prewarm script moves from `deploy/docker/` to `tools/` and is
shared: the Dockerfile keeps calling it, and omicsDesk's `inst/genesets`
is built by it. The existing deploy-contract test that checks the
collection map against `DB_MSIGDBR_MAP` covers both consumers.

### Phase 3: the bundle builder and its CI (4-5 days)

`tools/build_bundle.R --platform windows --r-version 4.5 --out dist/`:

1. `R CMD build` the three packages into `repo/src/contrib`.
2. `pak::pkg_download("local::packages/omicsDesk", dependencies = TRUE,
   platforms = ..., r_versions = ...)` into `repo/bin/<platform>/contrib/<minor>`
   and `repo/src/contrib`. pak resolves Bioconductor from the release
   matching the R version.
3. `tools::write_PACKAGES()` on each contrib directory.
4. Copy a pandoc release into `pandoc/`.
5. Write `INSTALL.R`, `START.R`, launchers, both READMEs, `bundle.json`,
   `SHA256SUMS`; zip.

One bundle is built for exactly one R minor version, because binary
packages are tied to it. Bioconductor publishes binaries only for the
current and previous release, so bundles are built while those exist and
the artefacts are kept; a bundle for an old R cannot be rebuilt later.

`.github/workflows/bundle.yaml`, on tag `desk-v*`: a matrix of
windows-latest, macos-latest (arm64), macos-13 (x86_64) and ubuntu-latest
builds and uploads the zips as release assets. A second job installs the
Linux bundle inside `docker run --network none`, runs `omicsDesk::doctor()`
and a headless ORA + GSEA + report, and fails on any error. That job is
the only place the whole promise is tested end to end, so it is not
optional.

For Linux users without compilers, p3m distro binaries can be added to
the matrix later; source tarballs cover the rest.

### Phase 4: users who do not have R (optional, 3-5 days)

The bundle can carry the official R installer for its platform (one
file, ~100 MB). The README then reads: install R, double-click INSTALL,
double-click omicsDesk. No RTools is needed because every package is a
binary. A truly portable Windows build (R-Portable plus a pre-installed
library, no install step) is a later refinement, not a prerequisite.

### Phase 5: updates and data refresh

* A new omicsDesk version is a tag and a rebuilt bundle; `.omp` files
  keep opening because the corpus test already guards the schema.
* Gene sets: a networked colleague runs `refresh_geneset_cache()` and
  `export_geneset_cache()`; the offline user runs
  `import_geneset_cache()`. `geneset_status()` shows the vintage in use,
  and `run_enrichment()` already records `gs_source` in every bundle.

## 5. Testing summary

| Level | Where | Proves |
|---|---|---|
| unit | omicsCore, omicsApp offline-contract tests | no code path reaches the network when offline mode is on |
| package | omicsDesk tests | launch, doctor, cache import behave |
| install | bundle CI, `--network none` | the zip installs and runs on a machine that has never seen the internet |
| browser | existing shinytest2 journey, run under the offline contract | the seven views work from the installed bundle |

## 6. Licensing notes

* MSigDB collections are distributed by msigdbr under MSigDB's terms
  (CC BY 4.0 for most collections; KEGG_LEGACY and BioCarta carry their
  own). Confirm the current terms before the first external release and
  record them in `inst/genesets/MANIFEST.json`.
* KEGG REST-derived caches are for local use only and must never be in
  a bundle; `export_geneset_cache()` should refuse tables whose
  `gs_source` starts with "KEGG REST" unless explicitly asked.
* pandoc is GPL-2; shipping the binary alongside is fine, and
  `bundle.json` names its version.

## 7. Name

**omicsDesk.** It completes the family (omicsCore, omicsApp, omicsDesk),
says "desktop, one machine" without saying "offline" as a limitation, is
short enough to type, and was free on CRAN and Bioconductor on
2026-09-05.

Considered and set aside: omicsBox (BioBam's commercial OmicsBox),
omicsStudio (collides with many "Omics Studio" web tools), omicsHub
(implies a network), omicsOffline (names the constraint, not the
product). Runners-up if omicsDesk does not sit right: omicsSolo,
omicsCrate, omicsBench, all free today.

## 8. Decisions needed before phase 2

1. Which platforms and R minor version ship first. The answer is
   whatever the intended users' laptops run; each extra platform is a CI
   matrix row, not new code.
2. Whether PDF reports matter offline. LaTeX is ~1.5 GB; HTML plus the
   browser's print-to-PDF is the cheap answer.
3. Whether gene sets stay inside omicsDesk or move to a separate data
   package once they need to be updated independently of the code.
