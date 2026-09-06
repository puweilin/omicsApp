# Normalize a proteomics `omics_input`

Turns linear intensities into analysis-ready values and records that in
`assay_type`. Every differential backend reads `expr_mat` as-is – limma
applies no transform at all – so this step is what makes the numbers
comparable across samples. Without it, limma runs on raw instrument
output: heavily right-skewed, variance scaling with the mean, and fold
changes that are differences of linear intensities rather than log
ratios.

## Usage

``` r
normalize_omics(input, method = c("vsn", "log2"), offset = 1)
```

## Arguments

- input:

  An `omics_input` with `omics_type = "proteomics"` carrying linear
  intensities.

- method:

  `"vsn"` (default) or `"log2"`.

- offset:

  Added before the log for `method = "log2"`, so that zeros survive.
  Ignored by `"vsn"`.

## Value

The `omics_input` with `expr_mat` normalized, `raw_mat` holding the
input matrix, and `assay_type` set to `"normalized_intensity"`.

## Details

`"vsn"` is the method the legacy proteomics pipeline has always used. It
fits a variance-stabilising transform on the linear matrix and returns
glog2-scale values, which behave like log2 intensities at the high end
while staying finite near zero.

The original matrix is kept in `raw_mat`, so the un-normalized values
remain available for QC views.

## See also

Other omics_input:
[`LOG_SCALE_ASSAY_TYPES`](https://puweilin.github.io/omicsApp/omicsCore/reference/LOG_SCALE_ASSAY_TYPES.md),
[`SUPPORTED_ASSAY_TYPES`](https://puweilin.github.io/omicsApp/omicsCore/reference/SUPPORTED_ASSAY_TYPES.md),
[`check_assay_scale()`](https://puweilin.github.io/omicsApp/omicsCore/reference/check_assay_scale.md),
[`drop_meta_na()`](https://puweilin.github.io/omicsApp/omicsCore/reference/drop_meta_na.md),
[`infer_assay_type()`](https://puweilin.github.io/omicsApp/omicsCore/reference/infer_assay_type.md),
[`is_omics_input()`](https://puweilin.github.io/omicsApp/omicsCore/reference/is_omics_input.md),
[`omics_input()`](https://puweilin.github.io/omicsApp/omicsCore/reference/omics_input.md),
[`select_complete_cases()`](https://puweilin.github.io/omicsApp/omicsCore/reference/select_complete_cases.md),
[`subset_omics()`](https://puweilin.github.io/omicsApp/omicsCore/reference/subset_omics.md),
[`subset_omics_features()`](https://puweilin.github.io/omicsApp/omicsCore/reference/subset_omics_features.md),
[`subset_omics_samples()`](https://puweilin.github.io/omicsApp/omicsCore/reference/subset_omics_samples.md),
[`summarize_omics()`](https://puweilin.github.io/omicsApp/omicsCore/reference/summarize_omics.md),
[`validate_omics_input()`](https://puweilin.github.io/omicsApp/omicsCore/reference/validate_omics_input.md)

## Examples

``` r
expr <- matrix(2^rnorm(200, 20, 2), nrow = 50, dimnames = list(
  paste0("P", 1:50), paste0("s", 1:4)
))
meta <- data.frame(group = c("A", "A", "B", "B"), row.names = paste0("s", 1:4))
feat <- data.frame(feature_id = paste0("P", 1:50), row.names = paste0("P", 1:50))
input <- omics_input(expr, meta, feat, omics_type = "proteomics",
                     assay_type = "raw_intensity")
normalized <- normalize_omics(input, method = "log2")
#> Normalized 50 features x 4 samples with log2
normalized$assay_type
#> [1] "normalized_intensity"
```
