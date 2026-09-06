# Detect sample-level outliers

Outlier detection on an
[`omics_input()`](https://puweilin.github.io/omicsApp/omicsCore/reference/omics_input.md).
Three methods are supported:

## Usage

``` r
qc_outliers(input, method = c("pca", "connectivity", "iqr"), sd_threshold = 3)
```

## Arguments

- input:

  An `omics_input`.

- method:

  One of `"pca"`, `"connectivity"`, `"iqr"`. Multiple methods may be
  supplied to run them in parallel and take the union of the flagged
  sets.

- sd_threshold:

  Z-score or IQR multiplier used to flag outliers (default `3`).

## Value

When `method` is length 1, a list with `method`, `stats`,
`flagged_samples`. When length \> 1, the same shape with `stats`
row-bound across methods and a `by_method` field holding per-method
results.

## Details

- `"pca"` — flags samples whose absolute z-score on PC1 or PC2 exceeds
  `sd_threshold`.

- `"connectivity"` — flags samples whose mean inter-sample correlation
  is more than `sd_threshold` standard deviations below the cohort mean
  (i.e. poorly-connected samples).

- `"iqr"` — flags samples whose mean log-intensity falls outside
  `[Q1 - k*IQR, Q3 + k*IQR]`; `sd_threshold` is reused as `k`.

All three methods are pure base R; no Suggests packages required.

## See also

Other qc:
[`IMPUTE_METHOD_ASSUMPTION`](https://puweilin.github.io/omicsApp/omicsCore/reference/IMPUTE_METHOD_ASSUMPTION.md),
[`impute_matrix()`](https://puweilin.github.io/omicsApp/omicsCore/reference/impute_matrix.md),
[`plot_qc()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_qc.md),
[`qc_depth()`](https://puweilin.github.io/omicsApp/omicsCore/reference/qc_depth.md),
[`qc_depth_outliers()`](https://puweilin.github.io/omicsApp/omicsCore/reference/qc_depth_outliers.md),
[`qc_missingness()`](https://puweilin.github.io/omicsApp/omicsCore/reference/qc_missingness.md),
[`resolve_impute_method()`](https://puweilin.github.io/omicsApp/omicsCore/reference/resolve_impute_method.md),
[`run_qc()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_qc.md),
[`winsorize_counts()`](https://puweilin.github.io/omicsApp/omicsCore/reference/winsorize_counts.md)
