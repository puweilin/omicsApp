# The pairing two layers would be integrated on

What the Integration view shows before it runs anything. One row per
donor, one column per layer, and a `source` saying where the pairing
came from – because "these two samples are the same person" is an
assertion the reader should be able to check rather than infer.

## Usage

``` r
sample_pairing_preview(project, tag_a, tag_b)
```

## Arguments

- project:

  An
  [omics_project](https://puweilin.github.io/omicsApp/omicsCore/reference/omics_project.md).

- tag_a, tag_b:

  Experiment tags.

## Value

A list with `pairs` (data frame, possibly zero rows) and `source`, one
of `"linked"` (a pairing saved on the project), `"donor"` (a donor
column in both layers), `"sample_id"` (ids that match outright),
`"suggested"` (guessed from the ids) or `"none"`.

## See also

Other integrate:
[`derive_sample_link()`](https://puweilin.github.io/omicsApp/omicsCore/reference/derive_sample_link.md),
[`suggest_sample_link()`](https://puweilin.github.io/omicsApp/omicsCore/reference/suggest_sample_link.md)
