# Samples whose depth is far from the rest

Flagged on the ratio to the median rather than an absolute count: what
counts as a shallow library depends entirely on the experiment, and a
fixed threshold would be wrong for every study but one.

## Usage

``` r
qc_depth_outliers(depth_df, min_ratio = 0.3)
```

## Arguments

- depth_df:

  Output of
  [`qc_depth()`](https://puweilin.github.io/omicsApp/omicsCore/reference/qc_depth.md).

- min_ratio:

  Flag samples below this fraction of the median.

## Value

Character vector of sample ids.

## See also

Other qc:
[`IMPUTE_METHOD_ASSUMPTION`](https://puweilin.github.io/omicsApp/omicsCore/reference/IMPUTE_METHOD_ASSUMPTION.md),
[`impute_matrix()`](https://puweilin.github.io/omicsApp/omicsCore/reference/impute_matrix.md),
[`plot_qc()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_qc.md),
[`qc_depth()`](https://puweilin.github.io/omicsApp/omicsCore/reference/qc_depth.md),
[`qc_missingness()`](https://puweilin.github.io/omicsApp/omicsCore/reference/qc_missingness.md),
[`qc_outliers()`](https://puweilin.github.io/omicsApp/omicsCore/reference/qc_outliers.md),
[`resolve_impute_method()`](https://puweilin.github.io/omicsApp/omicsCore/reference/resolve_impute_method.md),
[`run_qc()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_qc.md),
[`winsorize_counts()`](https://puweilin.github.io/omicsApp/omicsCore/reference/winsorize_counts.md)
