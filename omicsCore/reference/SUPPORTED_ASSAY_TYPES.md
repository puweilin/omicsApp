# Assay types recognised per omics modality

A named list with one entry per
[SUPPORTED_OMICS_TYPES](https://puweilin.github.io/omicsApp/omicsCore/reference/omicsCore-constants.md)
value. The names are deliberately explicit about scale, because nothing
downstream re-derives it: `raw_intensity` is linear instrument output,
`normalized_intensity` has been variance-stabilised (vsn) and is
log-like, `raw_count` is untransformed RNA-seq counts that DESeq2/edgeR
model directly.

## Usage

``` r
SUPPORTED_ASSAY_TYPES
```

## Format

An object of class `list` of length 2.

## Details

Values outside this vocabulary are allowed –
[`validate_omics_input()`](https://puweilin.github.io/omicsApp/omicsCore/reference/validate_omics_input.md)
only warns – so callers can drive omicsCore with modalities it does not
ship dispatchers for (CyTOF `arcsinh_intensity`, for example).

## See also

Other omics_input:
[`LOG_SCALE_ASSAY_TYPES`](https://puweilin.github.io/omicsApp/omicsCore/reference/LOG_SCALE_ASSAY_TYPES.md),
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
