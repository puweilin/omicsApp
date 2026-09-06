# Standardized differential analysis result columns

All differential analysis backends (DESeq2, edgeR, limma, t-test, lm)
return a `data.frame` with these columns so that downstream
visualization and integration functions can operate uniformly.

## Usage

``` r
DIFF_RESULT_REQUIRED_COLS
```

## Format

An object of class `character` of length 17.
