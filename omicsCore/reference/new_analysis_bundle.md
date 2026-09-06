# Create an analysis bundle

An `analysis_bundle` is the standard return container for every analysis
function in `omicsCore` (QC, differential, enrichment, integration). It
bundles results, parameter provenance, artifact registry, and
informational messages so downstream code (and the Shiny UI) can consume
them uniformly.

## Usage

``` r
new_analysis_bundle(
  analysis_name,
  input_info = list(),
  params = list(),
  results = list(),
  artifacts = NULL,
  messages = character(0),
  warnings = character(0)
)
```

## Arguments

- analysis_name:

  Bundle name (e.g. `"run_diff"`).

- input_info:

  Named list describing inputs (omics_type, sample counts, etc.) for
  provenance.

- params:

  Named list of analysis parameters.

- results:

  Named list of result objects (data frames, matrices, plots).

- artifacts:

  Artifact registry; defaults to an empty registry.

- messages:

  Character vector of informational messages.

- warnings:

  Character vector of warning messages.

## Value

An object of class `analysis_bundle`.
