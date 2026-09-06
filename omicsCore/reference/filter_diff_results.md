# Filter standardized differential results

Subsets a standardized diff result table to features that pass a
significance cutoff (either adjusted or raw p-value), plus optional
absolute-effect and model-fit thresholds. Returns the filtered table
with `is_significant = TRUE` set on every retained row.

## Usage

``` r
filter_diff_results(
  result_df,
  p_cutoff = 0.05,
  p_preference = c("adjusted", "raw"),
  effect_cutoff = NULL,
  model_fit_cutoff = NULL
)
```

## Arguments

- result_df:

  A standardized diff result `data.frame`, as produced by
  [`run_diff()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_diff.md)
  or any of the backend functions.

- p_cutoff:

  P-value threshold; defaults to 0.05.

- p_preference:

  Whether to threshold on the adjusted or raw p-value.

- effect_cutoff:

  Optional absolute effect cutoff (e.g., 1 for \|log2FC\| \>= 1).

- model_fit_cutoff:

  Optional model-fit cutoff (e.g., adjusted R^2 threshold for continuous
  lm/limma fits).

## Value

Filtered standardized diff result `data.frame`.

## See also

Other diff:
[`SUPPORTED_DIFF_METHODS`](https://puweilin.github.io/omicsApp/omicsCore/reference/SUPPORTED_DIFF_METHODS.md),
[`applicable_diff_methods()`](https://puweilin.github.io/omicsApp/omicsCore/reference/applicable_diff_methods.md),
[`make_ranked_features()`](https://puweilin.github.io/omicsApp/omicsCore/reference/make_ranked_features.md),
[`plot_feature_expression()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_feature_expression.md),
[`plot_heatmap()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_heatmap.md),
[`plot_ma()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_ma.md),
[`plot_pca()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_pca.md),
[`plot_volcano()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_volcano.md),
[`run_diff()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_diff.md),
[`run_diff_continuous()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_diff_continuous.md)
