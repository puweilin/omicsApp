# omicsCore Export Manifest

> **Status:** Phase 0 design contract. None of these functions are implemented yet.
> This document defines the planned public API surface so that Phase 1 implementation
> work can proceed against a stable target.

The API is organized into 6 sections:

1. [Project & input](#1-project--input)
2. [Data inspection & QC](#2-data-inspection--qc)
3. [Differential expression](#3-differential-expression)
4. [Pathway enrichment](#4-pathway-enrichment)
5. [Integration (dual-omics)](#5-integration-dual-omics)
6. [Visualization](#6-visualization)
7. [Persistence & I/O](#7-persistence--io)
8. [Installation helpers](#8-installation-helpers)

Conventions:
- **Inputs.** All analysis functions take an `omics_input` or `omics_project` as first argument.
- **Outputs.** Analysis functions return an `analysis_bundle` (S3). Plot functions return `ggplot`, `plotly`, or `ComplexHeatmap` objects.
- **Side effects.** No function writes to disk unless its name starts with `export_` or `save_`.
- **Errors.** Validation errors raise structured conditions via `rlang::abort(class = "omicsCore_error_*")`.

---

## 1. Project & input

### `omics_project(name, experiments = list(), sample_link = NULL, feature_link = NULL)`

Construct a multi-omics project container.

| Argument | Type | Description |
|---|---|---|
| `name` | character(1) | Human-readable project label |
| `experiments` | named list of `omics_input` | Each entry is one omics layer; names are arbitrary tags (e.g. `"proteomics"`, `"rnaseq"`) |
| `sample_link` | data.frame or NULL | Cross-omics sample ID mapping (columns: tag, sample_id, donor_id) |
| `feature_link` | data.frame or NULL | Feature ID mapping (e.g. uniprot ↔ gene_symbol) |

Returns: `omics_project` (S3).

### `omics_input(expr_mat, meta_df, feature_df, omics_type, assay_type = NULL)`

Construct a single-omics input object. **Inherited from existing framework** — same contract as `scripts/frameworks/omics_core/R/core/omics_input.R`.

### `read_omics(path, type = c("auto", "excel", "csv", "salmon", "maxquant", "rds"), ...)`

Smart input reader. Returns a list with two elements:

- `input`: the constructed `omics_input` (or `NULL` if inference failed)
- `report`: an `ImportReport` describing inferred sheet roles, column types, confidence scores, warnings

If `type = "auto"`, the parser:
1. Scans all sheets/files
2. Classifies each as `matrix | metadata | feature_annot | unknown`
3. Detects the orientation (features × samples vs samples × features) and transposes if needed
4. Identifies ID columns by pattern (UniProt, ENSG, gene symbol)
5. Returns the `ImportReport` for the UI to display before committing

### `add_experiment(project, name, input)`

Add an `omics_input` to an `omics_project`. Returns the modified project.

### `subset_omics(input, samples = NULL, features = NULL)`

Subset by sample or feature IDs. Wraps the existing `subset_omics_samples` / `subset_omics_features`.

---

## 2. Data inspection & QC

### `summarize_omics(input)`

Returns a `tibble` with one row summarizing: n_samples, n_features, n_missing, missing_pct, omics_type, assay_type.

### `run_qc(input, missing_threshold = 0.5, impute_method = c("none", "missforest", "knn", "min"), outlier_method = c("none", "pca", "iqr"), ...)`

Runs quality control pipeline. Returns `analysis_bundle` with:
- `bundle$results$qc_summary`: per-sample / per-feature missingness, outlier flags
- `bundle$results$cleaned_input`: post-QC `omics_input`
- `bundle$results$plots`: list of QC ggplots (missingness pattern, PCA, sample connectivity)

---

## 3. Differential expression

### `run_diff(input, design, contrast, method = c("auto", "deseq2", "limma", "edger", "ttest", "lm"), covariates = NULL, ...)`

Dispatcher for differential analysis. Auto-selects method based on `input$omics_type`:
- RNA-seq counts → DESeq2 (default) or edgeR
- Proteomics intensities → limma
- Continuous outcome → lm with splines
- Fallback → Welch t-test

| Argument | Type | Description |
|---|---|---|
| `design` | formula or character | e.g. `~ Group + Age` or `"Group"` |
| `contrast` | list | e.g. `list(group_col = "Group", control = "Young", case = "Old")` |
| `method` | character | Force a specific backend |
| `covariates` | character(n) or NULL | Adjustment variables |

Returns `analysis_bundle` with:
- `bundle$results$diff_result_df`: standardized schema (see existing `core/analysis_result_schema.R`)
- `bundle$results$diff_object`: backend-native object (DESeqDataSet, MArrayLM, etc.)
- `bundle$params`: full parameter record for provenance

### `run_diff_continuous(input, continuous_col, covariates = NULL, smooth = FALSE, ...)`

Convenience wrapper for continuous-outcome differential analysis (e.g. age trajectories).

---

## 4. Pathway enrichment

### `run_enrichment(diff_bundle, type = c("ora", "gsea", "gsva"), database = c("hallmark", "kegg", "reactome", "go_bp", "go_mf", "go_cc", "wikipathways"), organism = "Hs", direction = c("both", "up", "down"), p_cutoff = 0.05, ...)`

Returns `analysis_bundle` with:
- `bundle$results$enrich_result_df`: standardized schema
- `bundle$results$enrich_object`: clusterProfiler `enrichResult` / `gseaResult`

### `run_gsva(input, gene_sets = NULL, database = "hallmark", method = c("gsva", "ssgsea"), ...)`

Returns `analysis_bundle` with:
- `bundle$results$gsva_matrix`: pathways × samples
- `bundle$results$gsva_diff`: per-pathway diff result if a contrast was provided

### `list_gene_sets(database, organism = "Hs")`

Returns a `tibble` of available gene sets in a database. Uses cached MSigDB internally.

---

## 5. Integration (dual-omics)

### `run_integration(project, experiments = NULL, contrast, method = c("correlation", "concordance", "active_pathways"), ...)`

Cross-omics integration. `experiments = NULL` uses all experiments in the project; otherwise a character vector of tags.

| `method` | Output |
|---|---|
| `"correlation"` | Per-feature Pearson/Spearman correlation of expression across paired samples |
| `"concordance"` | Diff-direction concordance: per-feature log2FC in each omics, with quadrant classification |
| `"active_pathways"` | Combined-p-value pathway enrichment across omics (ActivePathways) |

Returns `analysis_bundle` with `bundle$results$integration_df` and (depending on method) `integration_matrix`, `quadrant_table`, etc.

---

## 6. Visualization

All plot functions accept `interactive = FALSE` (returns ggplot) or `interactive = TRUE` (returns plotly htmlwidget).

| Function | Returns | Notes |
|---|---|---|
| `plot_volcano(bundle, top_n = 20, label_features = NULL, ...)` | ggplot/plotly | log2FC vs -log10(padj) |
| `plot_ma(bundle, ...)` | ggplot/plotly | M-A plot |
| `plot_heatmap(bundle, n_top = 50, scale = "row", annotation = NULL, ...)` | ComplexHeatmap | top features by significance |
| `plot_pca(input, color_by, shape_by = NULL, ...)` | ggplot/plotly | PCA scores plot |
| `plot_feature_expression(input, features, group_by, ...)` | ggplot | violin/box for a feature set |
| `plot_enrichment(bundle, top = 20, view = c("dot", "bar", "ridge", "network"), ...)` | ggplot | enrichment summary |
| `plot_gsea(bundle, pathway_id, ...)` | ggplot | GSEA running-score curve |
| `plot_gsva_heatmap(bundle, ...)` | ComplexHeatmap | pathways × samples |
| `plot_integration(int_bundle, view = c("dual_volcano", "scatter", "dotplot", "quadrant"), ...)` | ggplot/plotly | dual-omics views |
| `plot_qc(qc_bundle, view = c("missing", "pca", "connectivity", "imputation"), ...)` | ggplot | QC panels |

### `theme_omicsCore(base_size = 11, base_family = "Inter")`

Project-standard ggplot2 theme. Matches the omicsApp visual identity.

---

## 7. Persistence & I/O

### `save_project(project, path)`

Writes an `omics_project` (with all bundles and artifacts) to a `.omp` file using `qs2::qs_save()`. Atomic write (writes to `path.tmp` then renames).

### `load_project(path)`

Loads a `.omp` file. Validates schema version.

### `export_bundle(bundle, dir, formats = c("xlsx", "tsv", "pdf"), prefix = NULL)`

Writes bundle artifacts to disk:
- Tables → `xlsx` and/or `tsv`
- Plots → `pdf` and/or `png`
- Provenance → `bundle_params.json`

Returns a `tibble` of written file paths (artifact registry).

### `export_report(project, path, format = c("html", "pdf"), template = "default")`

Renders a multi-panel report Rmd using project bundles. HTML uses self-contained mode (offline).

---

## 8. Installation helpers

### `install_optional(group = c("rnaseq", "proteomics", "enrichment", "imputation", "all"))`

Wraps `pak::pkg_install()` (or falls back to `BiocManager::install()` + `install.packages()`) to install groups of Suggests packages on demand.

| Group | Packages |
|---|---|
| `rnaseq` | DESeq2, edgeR, tximport, GenomicFeatures |
| `proteomics` | limma, missForest, pcaMethods |
| `enrichment` | clusterProfiler, msigdbr, fgsea, GSVA |
| `imputation` | missForest, pcaMethods |
| `viz` | ComplexHeatmap, circlize, ggrepel, patchwork, ggpubr |
| `all` | union of the above + `qs2` |

### `check_install(features = c("rnaseq", "proteomics", "enrichment"))`

Reports which optional groups are installed. Returns a `tibble`.

---

## Out of scope for v1

Not in this manifest, to be considered in v2:

- Olink, hormone, ISF modalities
- Single-cell (`SingleCellExperiment`) inputs
- Methylation arrays
- WGCNA module discovery (legacy code retained in `omics_core/legacy/` but not re-exported)
- Aging clock model training (separate package)
- Multi-user / authentication
- Distributed compute backends
