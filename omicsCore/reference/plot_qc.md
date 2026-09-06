# QC visualizations

Builds standard ggplot panels from a
[`run_qc()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_qc.md)
bundle.

## Usage

``` r
plot_qc(
  bundle,
  view = c("missing", "depth", "pca", "connectivity", "imputation"),
  color_by = NULL,
  ...
)
```

## Arguments

- bundle:

  An `analysis_bundle` produced by
  [`run_qc()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_qc.md).

- view:

  One of `"missing"`, `"pca"`, `"connectivity"`, `"imputation"`.

- color_by:

  Optional name of a column in the cleaned input's `meta_df` used to
  color samples in the `"pca"` view.

- ...:

  Reserved for future arguments.

## Value

A `ggplot` object.

## See also

Other qc:
[`IMPUTE_METHOD_ASSUMPTION`](https://puweilin.github.io/omicsApp/omicsCore/reference/IMPUTE_METHOD_ASSUMPTION.md),
[`impute_matrix()`](https://puweilin.github.io/omicsApp/omicsCore/reference/impute_matrix.md),
[`qc_depth()`](https://puweilin.github.io/omicsApp/omicsCore/reference/qc_depth.md),
[`qc_depth_outliers()`](https://puweilin.github.io/omicsApp/omicsCore/reference/qc_depth_outliers.md),
[`qc_missingness()`](https://puweilin.github.io/omicsApp/omicsCore/reference/qc_missingness.md),
[`qc_outliers()`](https://puweilin.github.io/omicsApp/omicsCore/reference/qc_outliers.md),
[`resolve_impute_method()`](https://puweilin.github.io/omicsApp/omicsCore/reference/resolve_impute_method.md),
[`run_qc()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_qc.md),
[`winsorize_counts()`](https://puweilin.github.io/omicsApp/omicsCore/reference/winsorize_counts.md)
