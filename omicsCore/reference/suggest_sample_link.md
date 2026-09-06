# Suggest a sample link from the ids

Only when the stems pair one-to-one. `A-1` and `A-2` both stem to `A`,
which would pair two samples of one layer to the same donor and pair the
wrong people in the other – a suggestion worse than none, because the
result of acting on it still looks like a result.

## Usage

``` r
suggest_sample_link(project, tag_a, tag_b)
```

## Arguments

- project:

  An
  [omics_project](https://puweilin.github.io/omicsApp/omicsCore/reference/omics_project.md).

- tag_a, tag_b:

  Experiment tags.

## Value

A data frame shaped like
[`derive_sample_link()`](https://puweilin.github.io/omicsApp/omicsCore/reference/derive_sample_link.md),
or `NULL`.

## See also

Other integrate:
[`derive_sample_link()`](https://puweilin.github.io/omicsApp/omicsCore/reference/derive_sample_link.md),
[`sample_pairing_preview()`](https://puweilin.github.io/omicsApp/omicsCore/reference/sample_pairing_preview.md)
