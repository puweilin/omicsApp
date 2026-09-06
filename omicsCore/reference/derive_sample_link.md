# Build a sample link from donor columns

Build a sample link from donor columns

## Usage

``` r
derive_sample_link(project)
```

## Arguments

- project:

  An
  [omics_project](https://puweilin.github.io/omicsApp/omicsCore/reference/omics_project.md).

## Value

A data frame with `tag`, `sample_id`, `donor_id`, or `NULL` when fewer
than two layers name a donor.

## See also

Other integrate:
[`sample_pairing_preview()`](https://puweilin.github.io/omicsApp/omicsCore/reference/sample_pairing_preview.md),
[`suggest_sample_link()`](https://puweilin.github.io/omicsApp/omicsCore/reference/suggest_sample_link.md)
