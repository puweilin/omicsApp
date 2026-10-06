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
