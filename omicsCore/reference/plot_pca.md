# PCA scores plot for an omics_input

Mean-imputes NA cells, optionally log2-transforms raw counts, then runs
`prcomp(scale. = TRUE)` on samples-by-features and returns a scatter on
the first two principal components. Sample metadata supplies optional
`color_by` and `shape_by` aesthetics.

## Usage

``` r
plot_pca(input, color_by = NULL, shape_by = NULL, log2 = NULL)
```

## Arguments

- input:

  A validated `omics_input`.

- color_by:

  Optional column in `meta_df` used to color samples.

- shape_by:

  Optional column in `meta_df` used as point shape.

- log2:

  If `TRUE`, apply `log2(x + 1)` before PCA. Defaults to `TRUE` when
  `assay_type == "raw_count"`.

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
[`plot_volcano()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_volcano.md),
[`run_diff()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_diff.md),
[`run_diff_continuous()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_diff_continuous.md)
