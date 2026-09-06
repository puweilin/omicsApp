# Warn when an assay type is unknown for this modality

Mirrors how
[`validate_omics_input()`](https://puweilin.github.io/omicsApp/omicsCore/reference/validate_omics_input.md)
treats `omics_type`: a warning, not an error, so a modality omicsCore
ships no dispatcher for still works.

## Usage

``` r
validate_assay_type(assay_type, omics_type)
```

## Arguments

- assay_type:

  Assay type label.

- omics_type:

  Omics modality.

## Value

Invisibly `TRUE`.
