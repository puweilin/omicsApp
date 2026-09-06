# Guess an assay type from the values

A starting point for an import wizard, not a verdict: it separates
linear proteomics intensities from already-transformed ones by
magnitude, the same heuristic
[`check_assay_scale()`](https://puweilin.github.io/omicsApp/omicsCore/reference/check_assay_scale.md)
uses, and assumes counts for RNA-seq because that is what gets uploaded.
The caller is expected to show the guess and let the user correct it –
nothing else recovers the scale if this is wrong.

## Usage

``` r
infer_assay_type(expr_mat, omics_type)
```

## Arguments

- expr_mat:

  Expression matrix.

- omics_type:

  Omics modality.

## Value

A single assay type from
[SUPPORTED_ASSAY_TYPES](https://puweilin.github.io/omicsApp/omicsCore/reference/SUPPORTED_ASSAY_TYPES.md),
or `NA_character_` for a modality with no vocabulary.

## See also

Other omics_input:
[`LOG_SCALE_ASSAY_TYPES`](https://puweilin.github.io/omicsApp/omicsCore/reference/LOG_SCALE_ASSAY_TYPES.md),
[`SUPPORTED_ASSAY_TYPES`](https://puweilin.github.io/omicsApp/omicsCore/reference/SUPPORTED_ASSAY_TYPES.md),
[`check_assay_scale()`](https://puweilin.github.io/omicsApp/omicsCore/reference/check_assay_scale.md),
[`drop_meta_na()`](https://puweilin.github.io/omicsApp/omicsCore/reference/drop_meta_na.md),
[`is_omics_input()`](https://puweilin.github.io/omicsApp/omicsCore/reference/is_omics_input.md),
[`normalize_omics()`](https://puweilin.github.io/omicsApp/omicsCore/reference/normalize_omics.md),
[`omics_input()`](https://puweilin.github.io/omicsApp/omicsCore/reference/omics_input.md),
[`select_complete_cases()`](https://puweilin.github.io/omicsApp/omicsCore/reference/select_complete_cases.md),
[`subset_omics()`](https://puweilin.github.io/omicsApp/omicsCore/reference/subset_omics.md),
[`subset_omics_features()`](https://puweilin.github.io/omicsApp/omicsCore/reference/subset_omics_features.md),
[`subset_omics_samples()`](https://puweilin.github.io/omicsApp/omicsCore/reference/subset_omics_samples.md),
[`summarize_omics()`](https://puweilin.github.io/omicsApp/omicsCore/reference/summarize_omics.md),
[`validate_omics_input()`](https://puweilin.github.io/omicsApp/omicsCore/reference/validate_omics_input.md)

## Examples

``` r
infer_assay_type(matrix(2^rnorm(40, 20, 2), nrow = 10), "proteomics")
#> [1] "raw_intensity"
infer_assay_type(matrix(rpois(40, 200), nrow = 10), "rnaseq")
#> [1] "raw_count"
```
