# Winsorize per-gene outliers in a count matrix

For each gene (row), values exceeding `Q3 + k * IQR` are clipped to that
threshold. This removes extreme leverage points that can drive spurious
differential expression results while preserving the overall
distribution.

## Usage

``` r
winsorize_counts(count_mat, k = 20)
```

## Arguments

- count_mat:

  Integer or numeric count matrix (genes × samples).

- k:

  IQR multiplier for the upper fence. Larger values are more
  conservative (fewer values clipped). Default `20`.

## Value

A list with components:

- `count_mat`:

  Winsorized count matrix (same dimensions / names).

- `stats`:

  `data.frame` with per-gene outlier statistics.

- `n_clipped`:

  Total number of values clipped.

- `n_genes_affected`:

  Number of genes with at least one clipped value.

- `k`:

  The multiplier used.

## Details

Pure base R; no `matrixStats` dependency.

## See also

Other qc:
[`IMPUTE_METHOD_ASSUMPTION`](https://puweilin.github.io/omicsApp/omicsCore/reference/IMPUTE_METHOD_ASSUMPTION.md),
[`impute_matrix()`](https://puweilin.github.io/omicsApp/omicsCore/reference/impute_matrix.md),
[`plot_qc()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_qc.md),
[`qc_depth()`](https://puweilin.github.io/omicsApp/omicsCore/reference/qc_depth.md),
[`qc_depth_outliers()`](https://puweilin.github.io/omicsApp/omicsCore/reference/qc_depth_outliers.md),
[`qc_missingness()`](https://puweilin.github.io/omicsApp/omicsCore/reference/qc_missingness.md),
[`qc_outliers()`](https://puweilin.github.io/omicsApp/omicsCore/reference/qc_outliers.md),
[`resolve_impute_method()`](https://puweilin.github.io/omicsApp/omicsCore/reference/resolve_impute_method.md),
[`run_qc()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_qc.md)
