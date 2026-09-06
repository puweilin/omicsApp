# Assay types that already sit on a log-like scale

Used by
[`check_assay_scale()`](https://puweilin.github.io/omicsApp/omicsCore/reference/check_assay_scale.md)
to decide which direction a scale mismatch points in, and by callers
deciding whether a layer still needs
[`normalize_omics()`](https://puweilin.github.io/omicsApp/omicsCore/reference/normalize_omics.md).
Everything else in
[SUPPORTED_ASSAY_TYPES](https://puweilin.github.io/omicsApp/omicsCore/reference/SUPPORTED_ASSAY_TYPES.md)
is linear.

## Usage

``` r
LOG_SCALE_ASSAY_TYPES
```

## Format

An object of class `character` of length 5.

## See also

Other omics_input:
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
[`summarize_omics()`](https://puweilin.github.io/omicsApp/omicsCore/reference/summarize_omics.md),
[`validate_omics_input()`](https://puweilin.github.io/omicsApp/omicsCore/reference/validate_omics_input.md)
