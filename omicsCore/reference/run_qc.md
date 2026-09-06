# Run the standard quality-control pipeline

End-to-end QC for an
[`omics_input()`](https://puweilin.github.io/omicsApp/omicsCore/reference/omics_input.md).
Computes missingness, detects sample-level outliers, optionally imputes
the expression matrix, and returns an
[analysis_bundle](https://puweilin.github.io/omicsApp/omicsCore/reference/is_analysis_bundle.md)
containing the cleaned input plus the QC summary so the result can be
plotted by
[`plot_qc()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_qc.md)
or consumed directly by downstream analysis functions.

## Usage

``` r
run_qc(
  input,
  missing_threshold = 0.5,
  sample_missing_threshold = NULL,
  impute_method = NULL,
  outlier_method = NULL,
  outlier_sd_threshold = 3,
  ...
)
```

## Arguments

- input:

  An `omics_input`.

- missing_threshold:

  Feature missing-rate cutoff in `[0, 1]`. Features above this are
  flagged and removed from `cleaned_input`. Default `0.5`.

- sample_missing_threshold:

  Optional sample missing-rate cutoff. Samples above this are flagged
  and removed.

- impute_method:

  One of
  [IMPUTE_METHODS](https://puweilin.github.io/omicsApp/omicsCore/reference/IMPUTE_METHOD_ASSUMPTION.md)
  – DEP's method set. Applied to the expression matrix of
  `cleaned_input` after filtering. `NULL` (the default) resolves per
  modality via
  [`resolve_impute_method()`](https://puweilin.github.io/omicsApp/omicsCore/reference/resolve_impute_method.md):
  `"MinProb"` for proteomics, `"none"` for counts.

- outlier_method:

  One of `"none"`, `"pca"`, `"connectivity"`, `"iqr"`, or a vector of
  those (other than `"none"`) to union their flags.

- outlier_sd_threshold:

  Z-score / IQR multiplier passed to
  [`qc_outliers()`](https://puweilin.github.io/omicsApp/omicsCore/reference/qc_outliers.md).
  Default `3`.

- ...:

  Forwarded to the imputation backend.

## Value

An `analysis_bundle` with the following fields under `results`:

- `qc_summary`:

  List with `missingness`, `outliers`, and `recommended_filters`
  (sample/feature IDs to remove).

- `cleaned_input`:

  `omics_input` with flagged samples/features removed and (optionally)
  imputed expression matrix. `raw_mat` carries the pre-imputation matrix
  when imputation occurred.

## Details

Sensible defaults:

- `omics_type == "proteomics"` → `outlier_method = "pca"`,
  `impute_method = "MinProb"`. Missingness in DIA/DDA is mostly
  left-censored – a protein is absent because it fell below the
  detection limit – and leaving `NA` is not the neutral choice it looks
  like: limma drops what it cannot fit, so "none" is complete-case
  analysis taken silently.

- `omics_type == "rnaseq"` → `outlier_method = "connectivity"`,
  `impute_method = "none"`. A zero count is an observation, and imputing
  it feeds a negative-binomial model numbers it never saw.

Pass explicit arguments to override the defaults.

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
[`winsorize_counts()`](https://puweilin.github.io/omicsApp/omicsCore/reference/winsorize_counts.md)

## Examples

``` r
set.seed(1)
expr <- matrix(rnorm(60), nrow = 6,
               dimnames = list(paste0("g", 1:6), paste0("s", 1:10)))
expr[, 1] <- NA  # entirely missing sample
meta <- data.frame(group = rep(c("A", "B"), each = 5),
                   row.names = colnames(expr))
feat <- data.frame(feature_id = rownames(expr),
                   row.names = rownames(expr))
input <- omics_input(expr, meta, feat, omics_type = "proteomics")
#> Warning: `assay_type` is missing. It is the only record of what scale `expr_mat` is on; analyses read it to decide whether to transform. Set one of: 'raw_intensity', 'normalized_intensity', 'imputed_intensity', 'filtered_intensity', 'raw_count', 'tpm', 'fpkm', 'vst', 'logcpm'
bundle <- run_qc(input, sample_missing_threshold = 0.9,
                 outlier_method = "iqr")
#> Warning: `assay_type` is missing. It is the only record of what scale `expr_mat` is on; analyses read it to decide whether to transform. Set one of: 'raw_intensity', 'normalized_intensity', 'imputed_intensity', 'filtered_intensity', 'raw_count', 'tpm', 'fpkm', 'vst', 'logcpm'
#> Warning: `assay_type` is missing. It is the only record of what scale `expr_mat` is on; analyses read it to decide whether to transform. Set one of: 'raw_intensity', 'normalized_intensity', 'imputed_intensity', 'filtered_intensity', 'raw_count', 'tpm', 'fpkm', 'vst', 'logcpm'
#> Warning: `assay_type` is missing. It is the only record of what scale `expr_mat` is on; analyses read it to decide whether to transform. Set one of: 'raw_intensity', 'normalized_intensity', 'imputed_intensity', 'filtered_intensity', 'raw_count', 'tpm', 'fpkm', 'vst', 'logcpm'
#> Warning: `assay_type` is missing. It is the only record of what scale `expr_mat` is on; analyses read it to decide whether to transform. Set one of: 'raw_intensity', 'normalized_intensity', 'imputed_intensity', 'filtered_intensity', 'raw_count', 'tpm', 'fpkm', 'vst', 'logcpm'
#> Warning: `assay_type` is missing. It is the only record of what scale `expr_mat` is on; analyses read it to decide whether to transform. Set one of: 'raw_intensity', 'normalized_intensity', 'imputed_intensity', 'filtered_intensity', 'raw_count', 'tpm', 'fpkm', 'vst', 'logcpm'
#> Warning: `assay_type` is missing. It is the only record of what scale `expr_mat` is on; analyses read it to decide whether to transform. Set one of: 'raw_intensity', 'normalized_intensity', 'imputed_intensity', 'filtered_intensity', 'raw_count', 'tpm', 'fpkm', 'vst', 'logcpm'
#> Warning: `assay_type` is missing. It is the only record of what scale `expr_mat` is on; analyses read it to decide whether to transform. Set one of: 'raw_intensity', 'normalized_intensity', 'imputed_intensity', 'filtered_intensity', 'raw_count', 'tpm', 'fpkm', 'vst', 'logcpm'
bundle
#> <analysis_bundle>
#>   analysis : run_qc 
#>   results  : qc_summary, cleaned_input 
#>   params   : 5 entries
#>   artifacts: 0 rows
```
