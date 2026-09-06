# Per-feature linear regression for a two-group comparison

Fits `expr ~ group_col [+ covariates]` per feature using
[`stats::lm`](https://rdrr.io/r/stats/lm.html), extracting the
case-group coefficient as the effect estimate. RNA-seq raw counts are
automatically log2(x+1) transformed. Internal — call via
[`run_diff()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_diff.md).

## Usage

``` r
run_lm_group(input, group_col, control_group, case_group, covariates = NULL)
```

## Arguments

- input:

  A validated `omics_input`.

- group_col:

  Group column in sample metadata.

- control_group:

  Control-group label.

- case_group:

  Case-group label.

- covariates:

  Optional character vector of covariate column names.

## Value

List with `results_raw`, `results_std`, `model_object` (`NULL`), and
`analysis_info`.
