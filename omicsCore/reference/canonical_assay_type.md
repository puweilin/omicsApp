# Rewrite a superseded assay-type spelling

Maps the names in
[DEPRECATED_ASSAY_TYPE_ALIASES](https://puweilin.github.io/omicsApp/omicsCore/reference/omicsCore-constants.md)
onto their replacements, warning so the caller can update. Anything else
passes through untouched, including labels outside
[SUPPORTED_ASSAY_TYPES](https://puweilin.github.io/omicsApp/omicsCore/reference/SUPPORTED_ASSAY_TYPES.md).

## Usage

``` r
canonical_assay_type(assay_type)
```

## Arguments

- assay_type:

  Assay type label, or `NULL`.

## Value

The canonical label, or `assay_type` unchanged.
