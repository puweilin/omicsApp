# Export an analysis_bundle to disk

Writes the data-bearing fields of a bundle (`results$*_df` /
`results$*_matrix`), provenance (`params`), and – if `plots` is
non-empty – the supplied ggplot/heatmap objects, in a target directory.

## Usage

``` r
export_bundle(
  bundle,
  dir,
  formats = c("xlsx", "tsv", "pdf"),
  prefix = NULL,
  plots = NULL,
  width = 7,
  height = 5
)
```

## Arguments

- bundle:

  An
  [`analysis_bundle`](https://puweilin.github.io/omicsApp/omicsCore/reference/is_analysis_bundle.md).

- dir:

  Output directory. Created if it does not exist.

- formats:

  Character vector of table/plot formats. Tables support `"xlsx"` and
  `"tsv"`; plots support `"pdf"` and `"png"`. Defaults to
  `c("xlsx", "tsv", "pdf")`.

- prefix:

  Optional file-name prefix. Useful when exporting multiple bundles to
  the same directory.

- plots:

  Optional named list of plot objects to write alongside the tables
  (e.g. `list(volcano = plot_volcano(bundle))`). Both `ggplot` and
  `ComplexHeatmap` objects are supported.

- width, height:

  Plot dimensions in inches. Defaults `7 x 5`.

## Value

A `data.frame` artifact registry (one row per written file).

## Details

Output layout under `dir`:

    <dir>/
      <prefix><bundle_name>_<table>.xlsx
      <prefix><bundle_name>_<table>.tsv
      <prefix><bundle_name>_<matrix>.tsv
      <prefix><bundle_name>_<plot>.pdf      (or .png)
      <prefix><bundle_name>_params.json

`prefix` defaults to `""` (no prefix). The bundle's existing artifact
registry is merged with the rows added by this call so downstream code
can keep tracking every file from one entry point.

## See also

Other persistence:
[`export_report()`](https://puweilin.github.io/omicsApp/omicsCore/reference/export_report.md),
[`export_script()`](https://puweilin.github.io/omicsApp/omicsCore/reference/export_script.md),
[`load_project()`](https://puweilin.github.io/omicsApp/omicsCore/reference/load_project.md),
[`save_project()`](https://puweilin.github.io/omicsApp/omicsCore/reference/save_project.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  d <- run_diff(input, method = "ttest", analysis_type = "group",
                group_col = "g", control_group = "ctrl", case_group = "case")
  export_bundle(d, dir = "out/", plots = list(volcano = plot_volcano(d)))
} # }
```
