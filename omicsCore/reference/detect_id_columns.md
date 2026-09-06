# Detect biological ID columns in a data.frame

Scans column names and the first non-missing value in each column; tags
any column that matches a UniProt / Ensembl / RefSeq / HGNC regex.

## Usage

``` r
detect_id_columns(df, max_check = 200L)
```

## Arguments

- df:

  A `data.frame`.

- max_check:

  Number of rows to sample per column when probing values. Defaults to
  200 to keep large workbooks responsive.

## Value

A `data.frame` with columns `column`, `pattern`, `match_rate`. Empty
data frame if no matches were found.

## See also

Other io:
[`classify_sheet_role()`](https://puweilin.github.io/omicsApp/omicsCore/reference/classify_sheet_role.md),
[`detect_orientation()`](https://puweilin.github.io/omicsApp/omicsCore/reference/detect_orientation.md),
[`import_report_sheets()`](https://puweilin.github.io/omicsApp/omicsCore/reference/import_report_sheets.md),
[`import_report_warnings()`](https://puweilin.github.io/omicsApp/omicsCore/reference/import_report_warnings.md),
[`is_import_report()`](https://puweilin.github.io/omicsApp/omicsCore/reference/is_import_report.md),
[`new_import_report()`](https://puweilin.github.io/omicsApp/omicsCore/reference/new_import_report.md),
[`read_omics()`](https://puweilin.github.io/omicsApp/omicsCore/reference/read_omics.md)
