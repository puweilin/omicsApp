# omicsCore 0.2.0

First versioned release.

* Differential analysis with limma, edgeR, DESeq2, t-test and linear
  models; several comparisons from one fit, paired designs (fixed block in
  limma), continuous variables, global tests, tximport length offsets.
  Vectorised t-test/linear-model paths and cached normalisation factors.
* Quality control: missing values, sequencing depth, PCA/connectivity/IQR
  outlier tests, imputation on the log scale.
* Enrichment: ORA, GSEA and GSVA against MSigDB collections with an
  on-disk and in-memory gene-set cache; enrichment across comparisons.
* Integration: concordance, correlation over paired samples,
  ActivePathways.
* Import: Excel/CSV/TSV with orientation detection, sample-sheet matching,
  aligner suffix stripping, `fread` fast path with a `read.table`
  fallback.
* Export: HTML/PDF report for biologists and a reproducible R script.
* Plot helpers: `wrap_label()`, `wrap_comparison()`, `effect_label()`;
  long computations report their steps through
  `options(omicsCore.progress = )`; optional parallel DESeq2 via
  `options(omicsCore.deseq2_workers = )`.

## Added during the 0.2.0 release work

* Import: non-UTF-8 text (Windows-1252, GBK/GB18030), gzip-compressed
  text, `.xlsm`, Salmon/RSEM/kallisto files (`read_quant_files()`, with
  tximport-style length offsets and tx2gene summarising),
  SummarizedExperiment/DESeqDataSet (`read_summarized_experiment()`),
  summary columns dropped, RNA-seq value scale inferred from the values,
  mouse Ensembl ids mapped to MGI symbols, GENCODE `_PAR_Y` rows handled;
  uploaded `.rds` files checked to hold data only.
* QC: a leave-one-out outlier test that works from 4 samples; group-aware
  missing-value filter (`missing_filter`); `winsorize_counts()` leaves
  on/off genes alone and keeps counts whole (`legacy = TRUE` for the old
  rule); MinProb falls back to MinDet when it cannot be estimated; QC
  results store a rebuild record instead of the cleaned matrix
  (`qc_cleaned_input()`), about 94% smaller.
* Normalisation: `normalize_omics(method = "log2", center = "median")`.
* Differential: a signed test statistic for every engine (`signed_stat`);
  edgeR global tests use tximport offsets; paired t-tests pair each
  comparison on its own; optional DESeq2 pre-filter (`prefilter = TRUE`);
  DESeq2 continuous analysis tested and several design errors explained;
  shared standardisation code.
* Enrichment: ORA tests up- and down-regulated genes separately by
  default (`direction = "both"` keeps the pooled test); optional
  correction across databases (`p_adjust_scope = "all"`); more species
  (`enrichment_species()`) with a case-insensitive fallback; GSEA `eps`
  and `n_perm_simple`.
* Integration: directional ActivePathways (DPM); a protein-to-gene
  feature link (`read_feature_link()`, `feature_pairing_preview()`) that
  keeps every isoform.
* Projects: results record their layer (`bundle_layer()`); `.omp` files
  carry a format version with migrations, and can be signed
  (`save_project(signing_key =)`, `load_project(untrusted = TRUE)`,
  `sign_project_file()`).
* `install_optional()` groups cover the packages moved to Suggests
  (limma, clusterProfiler, msigdbr) and a new `"io"` group.
* Integration figures: `plot_integration(view = "effect_pair")` draws
  each layer's effect on equal axes with the hits in both layers counted
  in the legend and the top ones named; the new `view = "top_hits"` lists
  the top hits with a dot for each layer's effect. `view =
  "dual_volcano"` is deprecated (it warns, and new exported scripts no
  longer call it): its x axis, the difference of the two effects, is each
  point's distance from the diagonal in `"effect_pair"`.
* Figures: the volcano and MA plots colour hits by direction (up red,
  down blue) with counts in the legend, on an x axis symmetric about 0;
  the report and exported script draw the volcano at the project's saved
  thresholds. The QC missing-value panel gives samples their own axis
  and shows features as a histogram with the filter cutoff and how many
  were removed. The ORA dot plot puts the fraction of pathway genes on
  x (no longer the overlap twice) and caps the significance colour
  scale so one extreme pathway cannot wash out the rest. The
  correlation scatter names at most eight significant features.
* The QC depth view names at most ten samples and ranks the rest, with
  shallow libraries in amber in both panels; QC cutoff lines stay
  visible over the bars. The enrichment bar view is coloured by
  direction and cuts a far-outlying bar short, marked and labelled with
  its true value.
* QC depth on a proteomics layer is total intensity, summed on the
  linear scale when the layer is log2 (summing logs hid an
  under-loaded sample), and is labelled as such ("Total intensity",
  "Features quantified", "low" rather than "shallow").
* `plot_enrichment_comparison()`: legends no longer collide with the
  subtitle, the filled/hollow key appears only when needed, and point
  size uses the chosen (adjusted) p with a capped scale. `plot_gsea()`
  is a plain ggplot running-score curve that can rebuild its ranking
  from the differential result. `plot_feature_expression()` scales
  counts to log2(CPM + 1) and colours groups as the PCA does;
  `plot_heatmap()` gains group blocks and a highlighted feature.
* Figures: one colour-blind-safe group palette (`group_palette()`),
  distinct from the up/down colours, used for groups everywhere; plain
  titles; `plot_diff_overlap(compact = TRUE)` for narrow screens;
  report figure heights follow their content; enrichment dot plot keys
  under the panel; readable ActivePathways and correlation plots.
  `enrichplot` is no longer a suggested package.
