# omicsApp 0.2.0

First versioned release. Highlights since the development snapshots:

* Import: separate sample-sheet upload, orientation and matrix-sheet
  confirmation, upload size limits and background parsing, per-session
  autosave, restoring a project brings every view back with its results.
* Quality control: three outlier tests by default, excluding flagged
  samples from a layer, missing values reported before and after
  imputation, step-by-step progress on large layers.
* Differential: several comparisons (all pairs or custom), paired designs,
  continuous variables, a global test, full result downloads (CSV/Excel),
  plots that stay readable with long group names and on phones.
* Enrichment: species selection, enrichment across comparisons, thresholds
  shared with the differential view.
* Report: a report written for biologists (methods paragraph, figures,
  readable p-values) and an R script that reproduces the analysis,
  downloadable with its data files.
* One vocabulary throughout (layer, log2FC, adjusted p); run buttons are
  disabled while their run is in flight; long runs show which step they
  are on and can be cancelled.
* `shiny_app()` returns the application object; `inst/app/app.R` no longer
  reaches into the namespace.

## Added during the 0.2.0 release work

* Several layers of one omics type, named by the user; replacing a layer
  can keep both; results stay with their own layer.
* The Differential and Enrichment views say when the controls no longer
  match the result on screen.
* Import accepts several Salmon/RSEM/kallisto files at once, `.gz` text,
  `.xlsm` and SummarizedExperiment/DESeqDataSet `.rds` files, shows the
  file encoding and the guessed value scale.
* QC: four outlier tests (leave-one-out added), group-aware missing-value
  filter, saved QC results restored with their settings.
* Enrichment: up/down lists separately by default, more species.
* Integration: directional ActivePathways, an optional protein-to-gene
  mapping table.
* Project files are signed and checked before they are opened.
* Internals: module servers split by card; tests use installed packages
  unless omicsApp itself runs from source (`OMICSAPP_TEST_CORE`).
* Integration: the "Mirrored volcano" card is gone; the concordance
  scatter names the top hits, and a new "Top hits in both layers" card
  shows each one's effect in both layers side by side.
* Differential: the volcano follows the p and |log2FC| controls, like
  the table beside it, and its card says what "significant" means.
* Differential: the volcano's top-hit labels are laid out in columns
  beside the points and never overlap ("Label top hits (up to 20)").
* QC: on a proteomics layer the depth panel is called "Intensity" and
  its caption gives the median total intensity, written short (114G).
* Every figure has its own download (PNG 300 dpi, PDF, SVG).
* Hovering a point in the PCA and the integration plots names it.
* Differential: select a gene (table or volcano) to see it by group; a
  heatmap of the hits; the volcano thins dense non-significant points
  above 5,000 features (60,000 features: 5 MB to 0.2 MB).
* Enrichment: select a pathway to see its GSEA curve, or for ORA its
  overlapping genes with their log2FC.
* Tables show p-values to 3 significant digits and sort numerically;
  the enrichment table fits its card. Hover on the enrichment dot plot.
  Card titles replace in-plot titles; phone layouts for the PCA key and
  the overlap plot; the heatmap is drawn when scrolled into view.
  Figure files carry the date; SVGs use svglite when installed.
