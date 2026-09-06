# Check that the assay values look like the declared scale

`assay_type` is the only record of what scale `expr_mat` is on, and the
analysis backends act on it without re-deriving anything: limma runs
straight on `expr_mat`, while `log2(x + 1)` is applied only for
`"raw_count"`. Mislabelling therefore changes results silently – feeding
vsn-normalised values to a `raw_intensity` path, or running limma on
linear intensities, produces plausible numbers that are wrong.

## Usage

``` r
check_assay_scale(x)
```

## Arguments

- x:

  An `omics_input`.

## Value

Invisibly, `TRUE` when the values match the declared scale, `FALSE` when
they do not.

## Details

The check is a magnitude heuristic: log2 proteomics intensities top out
around 35, linear intensities run to millions. It warns rather than
fails, because a small or unusual matrix can legitimately land either
side of
[MAX_PLAUSIBLE_LOG_SCALE_VALUE](https://puweilin.github.io/omicsApp/omicsCore/reference/omicsCore-constants.md).

Not called from
[`validate_omics_input()`](https://puweilin.github.io/omicsApp/omicsCore/reference/validate_omics_input.md)
on purpose: validation runs at every analysis entry point, and this
scans the matrix.
[`omics_input()`](https://puweilin.github.io/omicsApp/omicsCore/reference/omics_input.md)
runs it once at construction, which is where a wrong label enters.

## See also

Other omics_input:
[`LOG_SCALE_ASSAY_TYPES`](https://puweilin.github.io/omicsApp/omicsCore/reference/LOG_SCALE_ASSAY_TYPES.md),
[`SUPPORTED_ASSAY_TYPES`](https://puweilin.github.io/omicsApp/omicsCore/reference/SUPPORTED_ASSAY_TYPES.md),
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
