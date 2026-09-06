# MA plot for a diff bundle

Plots `effect` against `base_mean` so the user can spot effect-size
biases concentrated at low- or high-expressed features.

## Usage

``` r
plot_ma(
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

  Number of top features to label.

- label_features:

  Optional character vector of `feature_symbol` values to always label.

- p_basis:

  Whether to threshold on the adjusted or raw p-value.

- effect_threshold:

  Optional absolute-effect cutoff.

- p_threshold:

  P-value cutoff, or `NULL` to make no distinction.

## Value

A `ggplot` object.

## Details

Which features count as significant is decided here, from the thresholds
this call is given — the same contract as
[`plot_volcano()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_volcano.md),
and for the same reason:
[`run_diff()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_diff.md)
applies no cutoff, so its `is_significant` column is `NA` and has
nothing to colour by.

## See also

Other diff:
[`SUPPORTED_DIFF_METHODS`](https://puweilin.github.io/omicsApp/omicsCore/reference/SUPPORTED_DIFF_METHODS.md),
[`applicable_diff_methods()`](https://puweilin.github.io/omicsApp/omicsCore/reference/applicable_diff_methods.md),
[`filter_diff_results()`](https://puweilin.github.io/omicsApp/omicsCore/reference/filter_diff_results.md),
[`make_ranked_features()`](https://puweilin.github.io/omicsApp/omicsCore/reference/make_ranked_features.md),
[`plot_feature_expression()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_feature_expression.md),
[`plot_heatmap()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_heatmap.md),
[`plot_pca()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_pca.md),
[`plot_volcano()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_volcano.md),
[`run_diff()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_diff.md),
[`run_diff_continuous()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_diff_continuous.md)
