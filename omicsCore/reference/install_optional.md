# Install optional dependency groups

Installs the on-demand `Suggests` packages required by a given backend
group. Uses
[`pak::pkg_install()`](https://pak.r-lib.org/reference/pkg_install.html)
if available, otherwise falls back to
[`BiocManager::install()`](https://bioconductor.github.io/BiocManager/reference/install.html)
(for Bioconductor packages) and
[`utils::install.packages()`](https://rdrr.io/r/utils/install.packages.html)
(for CRAN packages).

## Usage

``` r
install_optional(
  group = c("rnaseq", "proteomics", "enrichment", "imputation", "viz", "persistence",
    "all"),
  ask = interactive(),
  upgrade = FALSE
)
```

## Arguments

- group:

  One of `"rnaseq"`, `"proteomics"`, `"enrichment"`, `"imputation"`,
  `"viz"`, `"persistence"`, or `"all"`.

- ask:

  If `TRUE` (default in interactive sessions), prompt before installing.

- upgrade:

  If `TRUE`, allow upgrading already-installed packages. Defaults to
  `FALSE`.

## Value

Invisibly returns a character vector of the packages that were targeted
for installation.

## Details

This function exists so that a fresh `install.packages("omicsCore")` can
stay small on restricted environments where Docker / system installs are
not available. Users opt in to heavy backends only when they need them.

## See also

Other install:
[`check_install()`](https://puweilin.github.io/omicsApp/omicsCore/reference/check_install.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  install_optional("rnaseq")
  install_optional("all", ask = FALSE)
} # }
```
