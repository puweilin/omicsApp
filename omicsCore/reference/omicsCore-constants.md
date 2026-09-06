# Package constants

Internal vocabulary used across `omicsCore`. Not exported.

Names are the old values, elements the replacement.
[`omics_input()`](https://puweilin.github.io/omicsApp/omicsCore/reference/omics_input.md)
rewrites them;
[`validate_omics_input()`](https://puweilin.github.io/omicsApp/omicsCore/reference/validate_omics_input.md)
accepts them with a warning so projects saved before the vocabulary
existed still load.

log2 proteomics intensities top out near 35; linear intensities run into
the millions.
[`check_assay_scale()`](https://puweilin.github.io/omicsApp/omicsCore/reference/check_assay_scale.md)
uses this to tell them apart.

[`check_assay_scale()`](https://puweilin.github.io/omicsApp/omicsCore/reference/check_assay_scale.md)
separates linear from log values by how large they get, which only works
for intensity-like assays. Counts, TPM and FPKM are legitimately small –
a gene with five reads is ordinary – so applying the same rule to them
would flag correct data.

## Usage

``` r
SUPPORTED_OMICS_TYPES

DEPRECATED_ASSAY_TYPE_ALIASES

MAX_PLAUSIBLE_LOG_SCALE_VALUE

SCALE_CHECKED_ASSAY_TYPES

SUPPORTED_DIFF_ANALYSIS_TYPES

SUPPORTED_PREFERENCE

SUPPORTED_ENRICH_PREFERENCE

DEFAULT_CACHE_MAX_AGE_HOURS

CANONICAL_DATABASES
```

## Format

An object of class `character` of length 2.

An object of class `character` of length 2.

An object of class `numeric` of length 1.

An object of class `character` of length 4.

An object of class `character` of length 3.

An object of class `character` of length 2.

An object of class `character` of length 3.

An object of class `numeric` of length 1.

An object of class `character` of length 6.

## Details

`"intensity"` was what the import view stamped on every proteomics
upload, and `"raw_counts"` appeared in the documentation while the code
only ever matched `"raw_count"`.
