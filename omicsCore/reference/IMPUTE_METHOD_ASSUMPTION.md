# Imputation methods, grouped by what they assume

`IMPUTE_METHOD_ASSUMPTION` maps each method to the missingness it
assumes; `IMPUTE_METHODS` is its names, in the order the app lists them,
and is what
[`impute_matrix()`](https://puweilin.github.io/omicsApp/omicsCore/reference/impute_matrix.md)
and
[`run_qc()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_qc.md)
accept.

## Usage

``` r
IMPUTE_METHOD_ASSUMPTION

IMPUTE_METHODS
```

## Format

A named character vector: names are the method, values are `"MNAR"`,
`"MAR"`, `"either"` or `"none"`.

An object of class `character` of length 11.

## See also

Other qc:
[`impute_matrix()`](https://puweilin.github.io/omicsApp/omicsCore/reference/impute_matrix.md),
[`plot_qc()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_qc.md),
[`qc_depth()`](https://puweilin.github.io/omicsApp/omicsCore/reference/qc_depth.md),
[`qc_depth_outliers()`](https://puweilin.github.io/omicsApp/omicsCore/reference/qc_depth_outliers.md),
[`qc_missingness()`](https://puweilin.github.io/omicsApp/omicsCore/reference/qc_missingness.md),
[`qc_outliers()`](https://puweilin.github.io/omicsApp/omicsCore/reference/qc_outliers.md),
[`resolve_impute_method()`](https://puweilin.github.io/omicsApp/omicsCore/reference/resolve_impute_method.md),
[`run_qc()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_qc.md),
[`winsorize_counts()`](https://puweilin.github.io/omicsApp/omicsCore/reference/winsorize_counts.md)
