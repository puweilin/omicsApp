# Volcano plot for a diff bundle

x-axis is `effect` (log2FC, beta, correlation, depending on the
backend); y-axis is `-log10(p_value)` by default. Significant features
(`is_significant`) are highlighted; the top `top_n` rows by (adjusted)
p-value are labelled. Optionally supply `label_features` to force-label
a specific set of feature symbols.

## Usage

``` r
plot_volcano(
  bundle,
  top_n = 20,
  label_features = NULL,
  p_basis = c("adjusted", "raw"),
  effect_threshold = NULL,
  p_threshold = 0.05
)
```

## Arguments

- bundle:

  An `analysis_bundle` produced by
  [`run_diff()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_diff.md).

- top_n:

  Number of top features to label, ranked by `p_basis`.

- label_features:

  Optional character vector of `feature_symbol` values to always label.

- p_basis:

  Which p-value column to use for the y-axis, `"adjusted"` or `"raw"`.

- effect_threshold:

  Vertical reference line for absolute effect.

- p_threshold:

  Horizontal reference p-value line (significance cutoff).

## Value

A `ggplot` object.

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
[`run_diff()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_diff.md),
[`run_diff_continuous()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_diff_continuous.md)
