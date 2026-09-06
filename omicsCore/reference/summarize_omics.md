# Summarize an omics_input or omics_project

Returns a one-row-per-input `tibble` reporting sample / feature counts,
missingness, and modality labels. This is the inspection counterpart to
[`run_qc()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_qc.md)
and is intended for quick pre-flight checks in the Shiny app's import
wizard.

## Usage

``` r
summarize_omics(x)
```

## Arguments

- x:

  An
  [`omics_input`](https://puweilin.github.io/omicsApp/omicsCore/reference/omics_input.md)
  or
  [`omics_project`](https://puweilin.github.io/omicsApp/omicsCore/reference/is_omics_project.md).

## Value

A
[`tibble::tibble`](https://tibble.tidyverse.org/reference/tibble.html).

## Details

For an `omics_input`, returns a single-row `tibble`. For an
`omics_project`, returns one row per experiment with an extra `tag`
column.

Columns: `tag` (project only), `omics_type`, `assay_type`, `n_samples`,
`n_features`, `n_missing`, `missing_pct`, `n_meta_cols`,
`feature_columns`.

## See also

Other omics_input:
[`LOG_SCALE_ASSAY_TYPES`](https://puweilin.github.io/omicsApp/omicsCore/reference/LOG_SCALE_ASSAY_TYPES.md),
[`SUPPORTED_ASSAY_TYPES`](https://puweilin.github.io/omicsApp/omicsCore/reference/SUPPORTED_ASSAY_TYPES.md),
[`check_assay_scale()`](https://puweilin.github.io/omicsApp/omicsCore/reference/check_assay_scale.md),
[`drop_meta_na()`](https://puweilin.github.io/omicsApp/omicsCore/reference/drop_meta_na.md),
[`infer_assay_type()`](https://puweilin.github.io/omicsApp/omicsCore/reference/infer_assay_type.md),
[`is_omics_input()`](https://puweilin.github.io/omicsApp/omicsCore/reference/is_omics_input.md),
[`normalize_omics()`](https://puweilin.github.io/omicsApp/omicsCore/reference/normalize_omics.md),
[`omics_input()`](https://puweilin.github.io/omicsApp/omicsCore/reference/omics_input.md),
[`select_complete_cases()`](https://puweilin.github.io/omicsApp/omicsCore/reference/select_complete_cases.md),
[`subset_omics()`](https://puweilin.github.io/omicsApp/omicsCore/reference/subset_omics.md),
[`subset_omics_features()`](https://puweilin.github.io/omicsApp/omicsCore/reference/subset_omics_features.md),
[`subset_omics_samples()`](https://puweilin.github.io/omicsApp/omicsCore/reference/subset_omics_samples.md),
[`validate_omics_input()`](https://puweilin.github.io/omicsApp/omicsCore/reference/validate_omics_input.md)
