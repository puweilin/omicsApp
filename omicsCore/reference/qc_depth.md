# Per-sample depth and detection

Per-sample depth and detection

## Usage

``` r
qc_depth(input)
```

## Arguments

- input:

  An
  [omics_input](https://puweilin.github.io/omicsApp/omicsCore/reference/omics_input.md).

## Value

A data frame with one row per sample: `sample_id`, `library_size`,
`n_detected`, `detection_rate`, and `library_size_ratio` (each library
over the median, so "half the depth of a typical sample" is readable
without arithmetic).

## See also

Other qc:
[`IMPUTE_METHOD_ASSUMPTION`](https://puweilin.github.io/omicsApp/omicsCore/reference/IMPUTE_METHOD_ASSUMPTION.md),
[`impute_matrix()`](https://puweilin.github.io/omicsApp/omicsCore/reference/impute_matrix.md),
[`plot_qc()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_qc.md),
[`qc_depth_outliers()`](https://puweilin.github.io/omicsApp/omicsCore/reference/qc_depth_outliers.md),
[`qc_missingness()`](https://puweilin.github.io/omicsApp/omicsCore/reference/qc_missingness.md),
[`qc_outliers()`](https://puweilin.github.io/omicsApp/omicsCore/reference/qc_outliers.md),
[`resolve_impute_method()`](https://puweilin.github.io/omicsApp/omicsCore/reference/resolve_impute_method.md),
[`run_qc()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_qc.md),
[`winsorize_counts()`](https://puweilin.github.io/omicsApp/omicsCore/reference/winsorize_counts.md)
