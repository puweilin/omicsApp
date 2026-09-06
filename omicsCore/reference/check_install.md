# Report which optional dependency groups are installed

Returns a `data.frame` with one row per requested group listing total
package count, how many are installed, and the names of any missing
packages. Useful as a pre-flight check in the Shiny app.

## Usage

``` r
check_install(
  features = c("rnaseq", "proteomics", "enrichment", "imputation", "viz", "persistence")
)
```

## Arguments

- features:

  Character vector of group names to check.

## Value

A `data.frame` with columns `group`, `n_total`, `n_installed`,
`is_ready`, `missing`.

## See also

Other install:
[`install_optional()`](https://puweilin.github.io/omicsApp/omicsCore/reference/install_optional.md)
