# Which imputation a layer gets when the caller does not say

Proteomics gets `MinProb`; anything counted gets `none`. Mirrors how
[`run_qc()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_qc.md)
already resolves `outlier_method` per omics type, and for the same
reason: the right answer differs by measurement, and a single global
default is wrong for one of them whichever it is.

## Usage

``` r
resolve_impute_method(omics_type)
```

## Arguments

- omics_type:

  Modality string.

## Value

A single method name.

## See also

Other qc:
[`IMPUTE_METHOD_ASSUMPTION`](https://puweilin.github.io/omicsApp/omicsCore/reference/IMPUTE_METHOD_ASSUMPTION.md),
[`impute_matrix()`](https://puweilin.github.io/omicsApp/omicsCore/reference/impute_matrix.md),
[`plot_qc()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_qc.md),
[`qc_depth()`](https://puweilin.github.io/omicsApp/omicsCore/reference/qc_depth.md),
[`qc_depth_outliers()`](https://puweilin.github.io/omicsApp/omicsCore/reference/qc_depth_outliers.md),
[`qc_missingness()`](https://puweilin.github.io/omicsApp/omicsCore/reference/qc_missingness.md),
[`qc_outliers()`](https://puweilin.github.io/omicsApp/omicsCore/reference/qc_outliers.md),
[`run_qc()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_qc.md),
[`winsorize_counts()`](https://puweilin.github.io/omicsApp/omicsCore/reference/winsorize_counts.md)
