# Limma continuous-variable differential test

Either linear or spline (natural cubic via
[`splines::ns`](https://rdrr.io/r/splines/ns.html)) modeling of
`continuous_col`, returning per-feature limma statistics plus a
(partial) Spearman rho and adjusted R^2. Internal — call via
[`run_diff()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_diff.md).

## Usage

``` r
run_limma_continuous(
  input,
  continuous_col,
  method = c("linear", "spline"),
  df = 3,
  covariates = NULL,
  paired_col = NULL
)
```

## Arguments

- input:

  A validated `omics_input`.

- continuous_col:

  Continuous metadata column.

- method:

  Either `"linear"` or `"spline"`.

- df:

  Degrees of freedom for spline fits.

- covariates:

  Optional covariate column names.

- paired_col:

  Optional pairing/block column.

## Value

List with `results_raw`, `results_std`, `model_object`, and
`analysis_info`.
