# Low-level `omics_input` constructor

Builds an `omics_input` without running
[`validate_omics_input()`](https://puweilin.github.io/omicsApp/omicsCore/reference/validate_omics_input.md).
Prefer
[`omics_input()`](https://puweilin.github.io/omicsApp/omicsCore/reference/omics_input.md)
in user-facing code; this constructor is used internally when validation
is deferred (e.g. inside subset operations).

## Usage

``` r
new_omics_input(
  omics_type,
  assay_type,
  expr_mat,
  meta_df,
  feature_df,
  raw_mat = NULL,
  normalized_mat = NULL,
  raw_object = NULL,
  source_fingerprint = NULL,
  source_path = NULL
)
```

## Arguments

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

- expr_mat:

  Main expression matrix; rows are features, columns are samples. A
  `data.frame` is coerced to a matrix.

- meta_df:

  Sample metadata; `rownames(meta_df)` must include
  `colnames(expr_mat)`.

- feature_df:

  Feature metadata containing a `feature_id` column.

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
