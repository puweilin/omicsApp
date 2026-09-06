# Render an HTML report for an omics_project

Renders a self-contained HTML (or PDF) report summarising every
experiment and any `analysis_bundle` objects attached to a project. The
bundled `inst/rmd/default-report.Rmd` template is used unless a custom
path is supplied.

## Usage

``` r
export_report(
  project,
  path,
  format = c("html", "pdf"),
  template = NULL,
  overwrite = FALSE
)
```

## Arguments

- project:

  An
  [`omics_project`](https://puweilin.github.io/omicsApp/omicsCore/reference/is_omics_project.md).

- path:

  Output path (file extension is inferred from `format`).

- format:

  `"html"` (default) or `"pdf"`. PDF requires a working LaTeX engine in
  the user environment.

- template:

  Optional path to a custom Rmd template. Defaults to the package's
  `inst/rmd/default-report.Rmd`.

- overwrite:

  If `FALSE` (default) and `path` exists, raise an error.

## Value

Invisibly returns `path`.

## Details

Bundles can be attached to the project in a `bundles` slot (e.g.
`project$bundles$diff <- run_diff(...)`); the template walks over
`names(project$bundles)` and emits per-bundle parameter and top-rows
tables.

Requires the `rmarkdown` package.

## See also

Other persistence:
[`export_bundle()`](https://puweilin.github.io/omicsApp/omicsCore/reference/export_bundle.md),
[`export_script()`](https://puweilin.github.io/omicsApp/omicsCore/reference/export_script.md),
[`load_project()`](https://puweilin.github.io/omicsApp/omicsCore/reference/load_project.md),
[`save_project()`](https://puweilin.github.io/omicsApp/omicsCore/reference/save_project.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  p <- omics_project("demo", experiments = list(proteo = my_input))
  p$bundles <- list(diff = run_diff(my_input, ...))
  export_report(p, "demo.html")
} # }
```
