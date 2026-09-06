# Classify a sheet/data.frame's role in an omics import

Walks a small decision tree to label one sheet as `"matrix"`,
`"metadata"`, `"feature_annot"`, or `"unknown"`. Returns the confidence
in `[0, 1]` so callers can present alternatives.

## Usage

``` r
classify_sheet_role(df, name = NA_character_)
```

## Arguments

- df:

  A `data.frame`.

- name:

  Optional sheet name (used as a tiebreaker hint).

## Value

A list with `role`, `confidence`, `orientation`, `notes`.

## See also

Other io:
[`detect_id_columns()`](https://puweilin.github.io/omicsApp/omicsCore/reference/detect_id_columns.md),
[`detect_orientation()`](https://puweilin.github.io/omicsApp/omicsCore/reference/detect_orientation.md),
[`import_report_sheets()`](https://puweilin.github.io/omicsApp/omicsCore/reference/import_report_sheets.md),
[`import_report_warnings()`](https://puweilin.github.io/omicsApp/omicsCore/reference/import_report_warnings.md),
[`is_import_report()`](https://puweilin.github.io/omicsApp/omicsCore/reference/is_import_report.md),
[`new_import_report()`](https://puweilin.github.io/omicsApp/omicsCore/reference/new_import_report.md),
[`read_omics()`](https://puweilin.github.io/omicsApp/omicsCore/reference/read_omics.md)
