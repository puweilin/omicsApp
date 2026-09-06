# Per-feature expression box / violin

Plots one or more features as boxplots (with jittered points) split by a
grouping column. RNA-seq raw counts are log2-transformed automatically.

## Usage

``` r
plot_feature_expression(input, features, group_by, color_by = NULL)
```

## Arguments

- input:

  A validated `omics_input`.

- features:

  Character vector of feature IDs or symbols to plot.

- group_by:

  Column in `meta_df` used for the x-axis grouping.

- color_by:

  Optional column for additional color aesthetic.

## Value

A `ggplot` object.

## See also

Other diff:
[`SUPPORTED_DIFF_METHODS`](https://puweilin.github.io/omicsApp/omicsCore/reference/SUPPORTED_DIFF_METHODS.md),
[`applicable_diff_methods()`](https://puweilin.github.io/omicsApp/omicsCore/reference/applicable_diff_methods.md),
[`filter_diff_results()`](https://puweilin.github.io/omicsApp/omicsCore/reference/filter_diff_results.md),
[`make_ranked_features()`](https://puweilin.github.io/omicsApp/omicsCore/reference/make_ranked_features.md),
[`plot_heatmap()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_heatmap.md),
[`plot_ma()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_ma.md),
[`plot_pca()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_pca.md),
[`plot_volcano()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_volcano.md),
[`run_diff()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_diff.md),
[`run_diff_continuous()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_diff_continuous.md)
