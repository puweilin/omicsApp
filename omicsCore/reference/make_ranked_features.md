# Build a named ranked feature vector for preranked GSEA

Produces a named numeric vector sorted in decreasing order, suitable as
input to preranked GSEA (e.g., `fgsea::fgsea(stats = ...)`). Duplicate
feature labels are dropped, keeping the first occurrence after the sort.

## Usage

``` r
make_ranked_features(
  result_df,
  feature_col = "feature_symbol",
  rank_col = "effect"
)
```

## Arguments

- result_df:

  A standardized diff result `data.frame`.

- feature_col:

  Column used to name the vector (defaults to `"feature_symbol"`).

- rank_col:

  Numeric column used to rank features (defaults to `"effect"`).

## Value

Named numeric vector sorted decreasing.

## See also

Other diff:
[`SUPPORTED_DIFF_METHODS`](https://puweilin.github.io/omicsApp/omicsCore/reference/SUPPORTED_DIFF_METHODS.md),
[`applicable_diff_methods()`](https://puweilin.github.io/omicsApp/omicsCore/reference/applicable_diff_methods.md),
[`filter_diff_results()`](https://puweilin.github.io/omicsApp/omicsCore/reference/filter_diff_results.md),
[`plot_feature_expression()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_feature_expression.md),
[`plot_heatmap()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_heatmap.md),
[`plot_ma()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_ma.md),
[`plot_pca()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_pca.md),
[`plot_volcano()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_volcano.md),
[`run_diff()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_diff.md),
[`run_diff_continuous()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_diff_continuous.md)
