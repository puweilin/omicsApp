# Impute missing values in an expression matrix

The methods are DEP's, so a choice made here means the same thing it
means in the proteomics literature.

## Usage

``` r
impute_matrix(mat, method = IMPUTE_METHODS, ...)
```

## Arguments

- mat:

  A numeric matrix, features in rows.

- method:

  One of
  [IMPUTE_METHODS](https://puweilin.github.io/omicsApp/omicsCore/reference/IMPUTE_METHOD_ASSUMPTION.md).

- ...:

  Forwarded to the backend (`q` for MinProb/MinDet, `k` for knn,
  `shift`/`scale` for man).

## Value

A numeric matrix with the same dimensions and names as `mat`.

## Details

- `"none"` — leave `NA`. Not neutral: downstream this becomes
  complete-case analysis.

- `"MinProb"` — random draws from a narrow gaussian near the observed
  minimum (MNAR). The default for proteomics.

- `"MinDet"` — a low quantile of the observed values, deterministic
  (MNAR).

- `"QRILC"` — quantile regression for left-censored data (MNAR).

- `"min"` — the observed minimum of that feature (MNAR).

- `"zero"` — zero (MNAR). Offered for parity; it distorts variance.

- `"knn"` — k-nearest neighbours (MAR).

- `"MLE"` — maximum likelihood (MAR).

- `"bpca"` — Bayesian PCA (MAR).

- `"mixed"` — MAR or MNAR per feature, chosen by a test.

- `"man"` — manual shift/scale, as in DEP's `shift`/`scale` arguments.

MNAR methods assume log-scale data, which is the scale DEP's workflow
imputes on. On linear intensities they will place imputed values on a
scale the observed values are not on.

## See also

Other qc:
[`IMPUTE_METHOD_ASSUMPTION`](https://puweilin.github.io/omicsApp/omicsCore/reference/IMPUTE_METHOD_ASSUMPTION.md),
[`plot_qc()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_qc.md),
[`qc_depth()`](https://puweilin.github.io/omicsApp/omicsCore/reference/qc_depth.md),
[`qc_depth_outliers()`](https://puweilin.github.io/omicsApp/omicsCore/reference/qc_depth_outliers.md),
[`qc_missingness()`](https://puweilin.github.io/omicsApp/omicsCore/reference/qc_missingness.md),
[`qc_outliers()`](https://puweilin.github.io/omicsApp/omicsCore/reference/qc_outliers.md),
[`resolve_impute_method()`](https://puweilin.github.io/omicsApp/omicsCore/reference/resolve_impute_method.md),
[`run_qc()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_qc.md),
[`winsorize_counts()`](https://puweilin.github.io/omicsApp/omicsCore/reference/winsorize_counts.md)
