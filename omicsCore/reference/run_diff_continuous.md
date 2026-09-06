# Continuous-variable differential analysis

Convenience wrapper around
[`run_diff()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_diff.md)
that pins `analysis_type = "continuous"`.

## Usage

``` r
run_diff_continuous(
  input,
  method = "auto",
  continuous_col,
  covariates = NULL,
  paired_col = NULL,
  ...
)
```

## Arguments

- input:

  A validated `omics_input`.

- method:

  Backend name. `"auto"` (default) lets `omicsCore` pick one based on
  `omics_type` and installed Suggests.

- continuous_col:

  Continuous metadata column (continuous only).

- covariates:

  Optional character vector of covariate column names.

- paired_col:

  Optional pairing/block column.

- ...:

  Extra arguments forwarded to the backend, e.g. `var_equal` for t-test
  or `df = 3` for limma spline.

## Value

An
[`analysis_bundle`](https://puweilin.github.io/omicsApp/omicsCore/reference/is_analysis_bundle.md).

## See also

Other diff:
[`SUPPORTED_DIFF_METHODS`](https://puweilin.github.io/omicsApp/omicsCore/reference/SUPPORTED_DIFF_METHODS.md),
[`applicable_diff_methods()`](https://puweilin.github.io/omicsApp/omicsCore/reference/applicable_diff_methods.md),
[`filter_diff_results()`](https://puweilin.github.io/omicsApp/omicsCore/reference/filter_diff_results.md),
[`make_ranked_features()`](https://puweilin.github.io/omicsApp/omicsCore/reference/make_ranked_features.md),
[`plot_feature_expression()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_feature_expression.md),
[`plot_heatmap()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_heatmap.md),
[`plot_ma()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_ma.md),
[`plot_pca()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_pca.md),
[`plot_volcano()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_volcano.md),
[`run_diff()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_diff.md)
