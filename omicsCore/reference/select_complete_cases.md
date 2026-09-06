# Select features passing a missingness threshold

Select features passing a missingness threshold

## Usage

``` r
select_complete_cases(omics_input, feature_missing_cutoff = 1)
```

## Arguments

- omics_input:

  A validated `omics_input`.

- feature_missing_cutoff:

  Maximum allowed per-feature missing rate (0–1).

## Value

A new `omics_input` filtered by feature missingness.

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
[`subset_omics()`](https://puweilin.github.io/omicsApp/omicsCore/reference/subset_omics.md),
[`subset_omics_features()`](https://puweilin.github.io/omicsApp/omicsCore/reference/subset_omics_features.md),
[`subset_omics_samples()`](https://puweilin.github.io/omicsApp/omicsCore/reference/subset_omics_samples.md),
[`summarize_omics()`](https://puweilin.github.io/omicsApp/omicsCore/reference/summarize_omics.md),
[`validate_omics_input()`](https://puweilin.github.io/omicsApp/omicsCore/reference/validate_omics_input.md)
