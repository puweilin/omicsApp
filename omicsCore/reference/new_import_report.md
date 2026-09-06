# Construct an `ImportReport`

Side-output object emitted by
[`read_omics()`](https://puweilin.github.io/omicsApp/omicsCore/reference/read_omics.md).
Bundles a per-sheet classification table, free-form warnings, and a
`suggested_input` slot that captures the assignment used to (try to)
build an `omics_input`.

## Usage

``` r
new_import_report(
  sheets = NULL,
  warnings = character(0),
  suggested_input = list(),
  source = NA_character_
)
```

## Arguments

- sheets:

  A `data.frame` (or `tibble`) with one row per sheet/tab detected in
  the input file. Required columns: `name`, `role`, `n_rows`, `n_cols`,
  `confidence`, `orientation`, `notes`. The constructor coerces partial
  inputs and fills in missing columns.

- warnings:

  Character vector of warnings to surface to the user.

- suggested_input:

  Named list describing the assignment used to build the returned
  `omics_input`. Conventional fields: `matrix_sheet`, `metadata_sheet`,
  `feature_sheet`, `orientation`, `id_column`, `omics_type`,
  `assay_type`. May be empty.

- source:

  Optional file path or other origin label.

## Value

An object of class `ImportReport`.

## See also

Other io:
[`classify_sheet_role()`](https://puweilin.github.io/omicsApp/omicsCore/reference/classify_sheet_role.md),
[`detect_id_columns()`](https://puweilin.github.io/omicsApp/omicsCore/reference/detect_id_columns.md),
[`detect_orientation()`](https://puweilin.github.io/omicsApp/omicsCore/reference/detect_orientation.md),
[`import_report_sheets()`](https://puweilin.github.io/omicsApp/omicsCore/reference/import_report_sheets.md),
[`import_report_warnings()`](https://puweilin.github.io/omicsApp/omicsCore/reference/import_report_warnings.md),
[`is_import_report()`](https://puweilin.github.io/omicsApp/omicsCore/reference/is_import_report.md),
[`read_omics()`](https://puweilin.github.io/omicsApp/omicsCore/reference/read_omics.md)
