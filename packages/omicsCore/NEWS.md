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
