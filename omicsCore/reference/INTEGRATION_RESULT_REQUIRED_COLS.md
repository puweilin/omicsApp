# Standardized integration result columns

All
[`run_integration()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_integration.md)
methods return a `data.frame` with these columns so downstream
visualisation and reporting code can treat correlation, concordance, and
pathway integration uniformly.

## Usage

``` r
INTEGRATION_RESULT_REQUIRED_COLS
```

## Format

An object of class `character` of length 15.

## Details

`effect` carries the method-specific score (Pearson/Spearman r, the
difference of effect sizes for concordance, or the combined -log10(p)
for ActivePathways) and `effect_type` records which one. `direction` is
`"concordant"` / `"discordant"` for concordance, `"positive"` /
`"negative"` for correlation, and `"shared"` / `"unique"` for pathway
integration. `quadrant` is filled for concordance (e.g. `"up_up"`).
