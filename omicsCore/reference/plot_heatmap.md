# Top-features expression heatmap

Two dispatch modes:

## Usage

``` r
plot_heatmap(
  x,
  input = NULL,
  n_top = 50L,
  features = NULL,
  scale = c("row", "none", "column"),
  annotation_cols = NULL,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = NULL,
  title = NULL
)
```

## Arguments

- x:

  Either an
  [`omics_input`](https://puweilin.github.io/omicsApp/omicsCore/reference/omics_input.md)
  or an
  [`analysis_bundle`](https://puweilin.github.io/omicsApp/omicsCore/reference/is_analysis_bundle.md)
  from
  [`run_diff()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_diff.md).

- input:

  Required when `x` is an analysis_bundle: the `omics_input` whose
  `expr_mat` should be used for the tiles.

- n_top:

  Number of features to display.

- features:

  Optional character vector of `feature_id`s to force. Overrides
  `n_top`.

- scale:

  One of `"row"` (default), `"none"`, or `"column"`.

- annotation_cols:

  Optional character vector of columns from `input$meta_df` to surface
  as a column annotation.

- cluster_rows, cluster_cols:

  Whether to cluster rows / columns.

- show_rownames:

  If `NULL`, auto-decide based on row count.

- title:

  Plot title.

## Value

A
[`ComplexHeatmap::Heatmap`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
object when `ComplexHeatmap` is available, otherwise a `ggplot` tile
plot.

## Details

- **`omics_input` mode** – pass an `omics_input` as the first argument.
  The function selects the top `n_top` features by row variance (after
  `coerce_to_continuous()`) and draws a samples-by-features heatmap.

- **`analysis_bundle` mode** – pass a
  [`run_diff()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_diff.md)
  bundle and the corresponding `omics_input` via the `input` argument.
  Features are selected by ascending adjusted p-value.

When `ComplexHeatmap` is installed the function returns a
[`ComplexHeatmap::Heatmap`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html);
otherwise it falls back to a ggplot2 tile plot so the function is always
usable from a clean install.

## See also

Other diff:
[`SUPPORTED_DIFF_METHODS`](https://puweilin.github.io/omicsApp/omicsCore/reference/SUPPORTED_DIFF_METHODS.md),
[`applicable_diff_methods()`](https://puweilin.github.io/omicsApp/omicsCore/reference/applicable_diff_methods.md),
[`filter_diff_results()`](https://puweilin.github.io/omicsApp/omicsCore/reference/filter_diff_results.md),
[`make_ranked_features()`](https://puweilin.github.io/omicsApp/omicsCore/reference/make_ranked_features.md),
[`plot_feature_expression()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_feature_expression.md),
[`plot_ma()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_ma.md),
[`plot_pca()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_pca.md),
[`plot_volcano()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_volcano.md),
[`run_diff()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_diff.md),
[`run_diff_continuous()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_diff_continuous.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  # From an omics_input
  plot_heatmap(input, n_top = 30)

  # From a diff bundle
  b <- run_diff(input, ...)
  plot_heatmap(b, input = input, n_top = 50)
} # }
```
