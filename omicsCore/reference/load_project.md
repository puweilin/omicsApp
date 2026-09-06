# Load an omics_project from disk

Reads an `.omp` archive written by
[`save_project()`](https://puweilin.github.io/omicsApp/omicsCore/reference/save_project.md).
Validates the schema version and that the payload is an `omics_project`.
Future schema changes will add migrations here.

## Usage

``` r
load_project(path)
```

## Arguments

- path:

  Path to the `.omp` file.

## Value

The deserialized `omics_project`.

## Details

Requires the `qs2` package; install it with
`install_optional("persistence")`.

## See also

Other persistence:
[`export_bundle()`](https://puweilin.github.io/omicsApp/omicsCore/reference/export_bundle.md),
[`export_report()`](https://puweilin.github.io/omicsApp/omicsCore/reference/export_report.md),
[`export_script()`](https://puweilin.github.io/omicsApp/omicsCore/reference/export_script.md),
[`save_project()`](https://puweilin.github.io/omicsApp/omicsCore/reference/save_project.md)
