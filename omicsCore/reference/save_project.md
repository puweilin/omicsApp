# Save an omics_project to disk

Writes the full project (every experiment, sample_link, and any attached
`analysis_bundle`s) to a single `.omp` file using `qs2`. The write is
atomic: the payload first goes to `path.tmp`, then renames over `path`,
so an interrupted save will not corrupt an existing file.

## Usage

``` r
save_project(project, path, overwrite = FALSE)
```

## Arguments

- project:

  An
  [`omics_project`](https://puweilin.github.io/omicsApp/omicsCore/reference/is_omics_project.md).

- path:

  Output path. The file extension is up to the caller, but the
  convention is `.omp` for projects produced by `omicsCore`.

- overwrite:

  If `FALSE` (default) and `path` already exists, `save_project()`
  raises an error. Set to `TRUE` to replace.

## Value

Invisibly returns `path`.

## Details

Requires the `qs2` package; install it with
`install_optional("persistence")`.

## See also

Other persistence:
[`export_bundle()`](https://puweilin.github.io/omicsApp/omicsCore/reference/export_bundle.md),
[`export_report()`](https://puweilin.github.io/omicsApp/omicsCore/reference/export_report.md),
[`export_script()`](https://puweilin.github.io/omicsApp/omicsCore/reference/export_script.md),
[`load_project()`](https://puweilin.github.io/omicsApp/omicsCore/reference/load_project.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  p <- omics_project("demo", experiments = list(proteo = my_input))
  save_project(p, "demo.omp")
  q <- load_project("demo.omp")
} # }
```
