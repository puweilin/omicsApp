# Limma two-group differential test

Fits `~ 0 + group_col [+ covariates]` and contrasts
`case_group - control_group`. Optionally accommodates a `paired_col` via
[`limma::duplicateCorrelation`](https://rdrr.io/pkg/limma/man/dupcor.html).
Internal — call via
[`run_diff()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_diff.md).

## Usage

``` r
run_limma_group(
  input,
  group_col,
  control_group,
  case_group,
  covariates = NULL,
  paired_col = NULL
)
```

## Arguments

- input:

  A validated proteomics-style `omics_input`.

- group_col:

  Group column in sample metadata.

- control_group:

  Control-group label.

- case_group:

  Case-group label.

- covariates:

  Optional covariate column names.

- paired_col:

  Optional pairing/block column.

## Value

List with `results_raw`, `results_std`, `model_object`, and
`analysis_info`.
