# Detect the orientation of a candidate expression matrix

Two-way guess: are samples in columns (rows = features) or in rows (rows
= samples)? The heuristic favours the axis whose labels look more like
biological identifiers (long, alphanumeric, no spaces) and tiebreaks on
the typical "many features, few samples" shape.

## Usage

``` r
detect_orientation(df)
```

## Arguments

- df:

  A `data.frame` whose first column may be an ID column.

## Value

A list with `orientation` (one of `"features_in_rows"`,
`"samples_in_rows"`, `"ambiguous"`), `confidence`, and `notes`.

## See also

Other io:
[`classify_sheet_role()`](https://puweilin.github.io/omicsApp/omicsCore/reference/classify_sheet_role.md),
[`detect_id_columns()`](https://puweilin.github.io/omicsApp/omicsCore/reference/detect_id_columns.md),
[`import_report_sheets()`](https://puweilin.github.io/omicsApp/omicsCore/reference/import_report_sheets.md),
[`import_report_warnings()`](https://puweilin.github.io/omicsApp/omicsCore/reference/import_report_warnings.md),
[`is_import_report()`](https://puweilin.github.io/omicsApp/omicsCore/reference/is_import_report.md),
[`new_import_report()`](https://puweilin.github.io/omicsApp/omicsCore/reference/new_import_report.md),
[`read_omics()`](https://puweilin.github.io/omicsApp/omicsCore/reference/read_omics.md)
