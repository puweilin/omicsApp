# Read an omics workbook and report classifier findings

Front door for the omicsApp import wizard. Reads an Excel workbook, CSV,
or saved R object, classifies each sheet/data frame, and tries to build
a candidate `omics_input`. Always returns both the built object (or
`NULL`) and a structured
[`ImportReport`](https://puweilin.github.io/omicsApp/omicsCore/reference/new_import_report.md)
so the UI can show what was detected before the user commits.

## Usage

``` r
read_omics(
  path,
  type = c("auto", "excel", "csv", "rds"),
  omics_type = NULL,
  assay_type = NULL,
  sheet_roles = NULL,
  ...
)
```

## Arguments

- path:

  Path to a file. The extension drives auto-detection.

- type:

  One of `"auto"`, `"excel"`, `"csv"`, or `"rds"`.

- omics_type:

  Optional omics modality. If `NULL`, callers (the Shiny wizard) are
  expected to set it after inspecting the report.

- assay_type:

  Optional assay semantic label.

- sheet_roles:

  Optional named character vector overriding what the classifier
  decided, as `c("<sheet name>" = "<role>")` with roles drawn from
  [IMPORT_REPORT_ROLES](https://puweilin.github.io/omicsApp/omicsCore/reference/IMPORT_REPORT_ROLES.md).
  Sheets left out keep their inferred role.

  The classifier is a heuristic, and a confident wrong answer is worse
  than an unconfident one: a metadata sheet read as the matrix produces
  an `omics_input` that analyses cleanly and means nothing. This is how
  a caller lets the user correct it without re-implementing the reader.

- ...:

  Forwarded to the underlying reader
  ([`readxl::read_excel`](https://readxl.tidyverse.org/reference/read_excel.html),
  [`utils::read.csv`](https://rdrr.io/r/utils/read.table.html)).

## Value

A list with two elements:

- `input`: an `omics_input` or `NULL`.

- `report`: an
  [`ImportReport`](https://puweilin.github.io/omicsApp/omicsCore/reference/new_import_report.md).

## See also

Other io:
[`classify_sheet_role()`](https://puweilin.github.io/omicsApp/omicsCore/reference/classify_sheet_role.md),
[`detect_id_columns()`](https://puweilin.github.io/omicsApp/omicsCore/reference/detect_id_columns.md),
[`detect_orientation()`](https://puweilin.github.io/omicsApp/omicsCore/reference/detect_orientation.md),
[`import_report_sheets()`](https://puweilin.github.io/omicsApp/omicsCore/reference/import_report_sheets.md),
[`import_report_warnings()`](https://puweilin.github.io/omicsApp/omicsCore/reference/import_report_warnings.md),
[`is_import_report()`](https://puweilin.github.io/omicsApp/omicsCore/reference/is_import_report.md),
[`new_import_report()`](https://puweilin.github.io/omicsApp/omicsCore/reference/new_import_report.md)
