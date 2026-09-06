# Single-omics input container

Constructs the canonical input object consumed by all `omicsCore`
analysis functions. The argument order is `expr_mat`, `meta_df`,
`feature_df`, `omics_type`, `assay_type` to match the contract in
`docs/export-manifest.md`.

## Usage

``` r
omics_input(
  expr_mat,
  meta_df,
  feature_df,
  omics_type,
  assay_type = NULL,
  raw_mat = NULL,
  normalized_mat = NULL,
  raw_object = NULL,
  source_fingerprint = NULL,
  source_path = NULL
)
```

## Arguments

- expr_mat:

  Main expression matrix; rows are features, columns are samples. A
  `data.frame` is coerced to a matrix.

- meta_df:

  Sample metadata; `rownames(meta_df)` must include
  `colnames(expr_mat)`.

- feature_df:

  Feature metadata containing a `feature_id` column.

- omics_type:

  Omics modality, one of
  [SUPPORTED_OMICS_TYPES](https://puweilin.github.io/omicsApp/omicsCore/reference/omicsCore-constants.md).

- assay_type:

  Assay semantic label, ideally one of the values listed for this
  `omics_type` in
  [SUPPORTED_ASSAY_TYPES](https://puweilin.github.io/omicsApp/omicsCore/reference/SUPPORTED_ASSAY_TYPES.md)
  (e.g. `"raw_intensity"`, `"normalized_intensity"`, `"raw_count"`).
  This is the only record of what scale `expr_mat` is on – nothing
  downstream re-derives it – so an inaccurate label silently changes
  what the analysis backends do. Superseded spellings in
  [DEPRECATED_ASSAY_TYPE_ALIASES](https://puweilin.github.io/omicsApp/omicsCore/reference/omicsCore-constants.md)
  are rewritten with a warning.

- raw_mat:

  Optional raw matrix carried alongside `expr_mat`.

- normalized_mat:

  Optional normalized matrix carried alongside `expr_mat`.

- raw_object:

  Optional upstream raw object (e.g. a `DESeqDataSet` or
  `SummarizedExperiment`) preserved for reference.

- source_fingerprint:

  Optional string identifying the file this input was parsed from.
  Callers that re-import into a live project use it to tell "the user
  picked the same file again" from "the data changed"; analyses computed
  on the previous data are only valid for the former. `NULL` for inputs
  built directly from matrices.

- source_path:

  Optional path to the archived copy of that file.
  [`export_script()`](https://puweilin.github.io/omicsApp/omicsCore/reference/export_script.md)
  points its
  [`read_omics()`](https://puweilin.github.io/omicsApp/omicsCore/reference/read_omics.md)
  line at it, which is what makes an exported script runnable rather
  than illustrative.

## Value

An object of class `omics_input`.

## Details

`expr_mat` must be a numeric matrix with sample column names and feature
row names. `meta_df` must have rownames matching the sample column names
of `expr_mat`. `feature_df` must contain a `feature_id` column.

## See also

Other omics_input:
[`LOG_SCALE_ASSAY_TYPES`](https://puweilin.github.io/omicsApp/omicsCore/reference/LOG_SCALE_ASSAY_TYPES.md),
[`SUPPORTED_ASSAY_TYPES`](https://puweilin.github.io/omicsApp/omicsCore/reference/SUPPORTED_ASSAY_TYPES.md),
[`check_assay_scale()`](https://puweilin.github.io/omicsApp/omicsCore/reference/check_assay_scale.md),
[`drop_meta_na()`](https://puweilin.github.io/omicsApp/omicsCore/reference/drop_meta_na.md),
[`infer_assay_type()`](https://puweilin.github.io/omicsApp/omicsCore/reference/infer_assay_type.md),
[`is_omics_input()`](https://puweilin.github.io/omicsApp/omicsCore/reference/is_omics_input.md),
[`normalize_omics()`](https://puweilin.github.io/omicsApp/omicsCore/reference/normalize_omics.md),
[`select_complete_cases()`](https://puweilin.github.io/omicsApp/omicsCore/reference/select_complete_cases.md),
[`subset_omics()`](https://puweilin.github.io/omicsApp/omicsCore/reference/subset_omics.md),
[`subset_omics_features()`](https://puweilin.github.io/omicsApp/omicsCore/reference/subset_omics_features.md),
[`subset_omics_samples()`](https://puweilin.github.io/omicsApp/omicsCore/reference/subset_omics_samples.md),
[`summarize_omics()`](https://puweilin.github.io/omicsApp/omicsCore/reference/summarize_omics.md),
[`validate_omics_input()`](https://puweilin.github.io/omicsApp/omicsCore/reference/validate_omics_input.md)

## Examples

``` r
expr <- matrix(rnorm(20), nrow = 4, dimnames = list(
  paste0("g", 1:4), paste0("s", 1:5)
))
meta <- data.frame(group = c("A", "A", "B", "B", "B"),
                   row.names = paste0("s", 1:5))
feat <- data.frame(feature_id = paste0("g", 1:4),
                   row.names = paste0("g", 1:4))
omics_input(expr, meta, feat, omics_type = "proteomics",
            assay_type = "normalized_intensity")
#> <omics_input>
#>   omics_type : proteomics 
#>   assay_type : normalized_intensity 
#>   features   : 4 
#>   samples    : 5 
#>   missing %  : 0.00 
```
