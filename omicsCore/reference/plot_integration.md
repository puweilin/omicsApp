# Integration summary plot

Visualises an `analysis_bundle` produced by
[`run_integration()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_integration.md).
The available views depend on the method that produced the bundle:

## Usage

``` r
plot_integration(
  bundle,
  view = c("scatter", "dual_volcano", "effect_pair", "quadrant", "dotplot"),
  top_n = 20L,
  label_features = NULL,
  p_cutoff = 0.05
)
```

## Arguments

- bundle:

  An
  [`analysis_bundle`](https://puweilin.github.io/omicsApp/omicsCore/reference/is_analysis_bundle.md)
  produced by
  [`run_integration()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_integration.md).

- view:

  One of `"scatter"`, `"dual_volcano"`, `"quadrant"`, `"dotplot"`.

- top_n:

  Number of features / pathways to label or display.

- label_features:

  Optional character vector of `feature_symbol` values to force-label
  (scatter / dual_volcano views).

- p_cutoff:

  Significance cutoff used for highlight color in scatter and
  dual_volcano views.

## Value

A `ggplot` object.

## Details

- `"scatter"` – always available. For `correlation` plots the
  correlation coefficient against `-log10(adj_p_value)`. For
  `concordance` plots `effect_a` vs `effect_b` reconstructed from the
  `effect` column (which carries `effect_a - effect_b`). For
  `active_pathways` plots `-log10(adj_p_value)` against pathway rank.

- `"dual_volcano"` – concordance-only. Plots `effect` (the difference of
  effects) on the x-axis against `-log10(p)` on the y-axis and colors by
  quadrant.

- `"quadrant"` – concordance-only. Bar count of the four
  `(direction_a, direction_b)` sign quadrants.

- `"dotplot"` – active_pathways-only. Dotplot of top pathways, with
  color = adjusted p-value and shape = shared / unique evidence.

## See also

Other integration:
[`run_integration()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_integration.md)
