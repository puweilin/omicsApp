# Per-feature linear regression against a continuous variable

Fits `expr ~ continuous_col [+ covariates]` per feature, returning the
coefficient on `continuous_col` as the effect, plus a (partial) Spearman
rank correlation. RNA-seq raw counts are automatically log2(x+1)
transformed. Internal — call via
[`run_diff()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_diff.md).

## Usage

``` r
run_lm_continuous(input, continuous_col, covariates = NULL)
```

## Arguments

- input:

  A validated `omics_input`.

- continuous_col:

  Continuous metadata column name.

- covariates:

  Optional character vector of covariate column names.

## Value

List with `results_raw`, `results_std`, `model_object` (`NULL`), and
`analysis_info`.
