# Per-sample / per-feature missingness summary

Computes per-sample and per-feature missing rates for an
[`omics_input()`](https://puweilin.github.io/omicsApp/omicsCore/reference/omics_input.md)
and flags samples / features whose missing rate exceeds the supplied
thresholds.

## Usage

``` r
qc_missingness(
  input,
  sample_missing_cutoff = NULL,
  feature_missing_cutoff = 0.5
)
```

## Arguments

- input:

  An `omics_input`.

- sample_missing_cutoff:

  Optional sample missing-rate threshold in `[0, 1]`. `NULL` (default)
  leaves all samples unflagged.

- feature_missing_cutoff:

  Feature missing-rate threshold in `[0, 1]`.

## Value

A list with:

- `sample_metrics`:

  `data.frame` with `sample_id`, `missing_rate`.

- `feature_metrics`:

  `data.frame` with `feature_id`, `missing_rate`.

- `flagged_samples`:

  Character vector of sample IDs exceeding `sample_missing_cutoff`.

- `flagged_features`:

  Character vector of feature IDs exceeding `feature_missing_cutoff`.

- `settings`:

  Echo of the thresholds used.

## See also

Other qc:
[`IMPUTE_METHOD_ASSUMPTION`](https://puweilin.github.io/omicsApp/omicsCore/reference/IMPUTE_METHOD_ASSUMPTION.md),
[`impute_matrix()`](https://puweilin.github.io/omicsApp/omicsCore/reference/impute_matrix.md),
[`plot_qc()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_qc.md),
[`qc_depth()`](https://puweilin.github.io/omicsApp/omicsCore/reference/qc_depth.md),
[`qc_depth_outliers()`](https://puweilin.github.io/omicsApp/omicsCore/reference/qc_depth_outliers.md),
[`qc_outliers()`](https://puweilin.github.io/omicsApp/omicsCore/reference/qc_outliers.md),
[`resolve_impute_method()`](https://puweilin.github.io/omicsApp/omicsCore/reference/resolve_impute_method.md),
[`run_qc()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_qc.md),
[`winsorize_counts()`](https://puweilin.github.io/omicsApp/omicsCore/reference/winsorize_counts.md)
