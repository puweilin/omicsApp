# Export a project's analysis history as a runnable R script

Renders the calls that produced a project's analyses, in the order they
depend on each other, as an R script that reproduces them. Arguments
appear as they were *resolved*: an analysis run with `method = "auto"`
is emitted with the engine that auto selected, so the script and the
report agree.

## Usage

``` r
export_script(project, path = NULL, include_plots = TRUE)
```

## Arguments

- project:

  An
  [`omics_project`](https://puweilin.github.io/omicsApp/omicsCore/reference/is_omics_project.md).

- path:

  Optional file to write to. The lines are returned either way.

- include_plots:

  Whether to append the plotting calls. On by default: the Shiny views
  draw through these same functions, so the figures the script produces
  are the figures the user was looking at. Set `FALSE` for a script that
  only recomputes the numbers.

## Value

Character vector of script lines, invisibly when `path` is given.

## Details

Whatever cannot be reconstructed faithfully is flagged with a `NOTE`
comment in the script rather than approximated — a script that reads
correctly but computes something else is worse than none.

The input line points at the archived upload when the project records
one (see `source_path` on
[`omics_input()`](https://puweilin.github.io/omicsApp/omicsCore/reference/omics_input.md));
otherwise it emits a placeholder path for the reader to fill in.

## See also

Other persistence:
[`export_bundle()`](https://puweilin.github.io/omicsApp/omicsCore/reference/export_bundle.md),
[`export_report()`](https://puweilin.github.io/omicsApp/omicsCore/reference/export_report.md),
[`load_project()`](https://puweilin.github.io/omicsApp/omicsCore/reference/load_project.md),
[`save_project()`](https://puweilin.github.io/omicsApp/omicsCore/reference/save_project.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  export_script(project, "reproduce.R")
} # }
```
