# Differential methods that are valid for an input

`SUPPORTED_DIFF_METHODS` lists every backend that exists; this lists the
ones whose assumptions the data actually meets.

## Usage

``` r
applicable_diff_methods(input, analysis_type = "group")
```

## Arguments

- input:

  An
  [`omics_input`](https://puweilin.github.io/omicsApp/omicsCore/reference/omics_input.md).

- analysis_type:

  One of
  [SUPPORTED_DIFF_ANALYSIS_TYPES](https://puweilin.github.io/omicsApp/omicsCore/reference/omicsCore-constants.md).

## Value

Character vector of method names, always including `"auto"`.

## Details

DESeq2 and edgeR model raw counts as negative binomial. Given continuous
intensities they do not refuse — DESeq2 rounds to integers ("converting
counts to integer mode") and reports p-values for a model the data never
fitted. Conversely, this package's limma backend does not apply voom, so
it has no business being handed raw counts. Both mistakes produce a
full, plausible result table and no error, which is the failure mode
worth engineering against.

Written next to `auto_select_diff_method()` on purpose: one decides what
to run by default and the other what a caller may choose, and the two
disagreeing would be its own bug.

## See also

Other diff:
[`SUPPORTED_DIFF_METHODS`](https://puweilin.github.io/omicsApp/omicsCore/reference/SUPPORTED_DIFF_METHODS.md),
[`filter_diff_results()`](https://puweilin.github.io/omicsApp/omicsCore/reference/filter_diff_results.md),
[`make_ranked_features()`](https://puweilin.github.io/omicsApp/omicsCore/reference/make_ranked_features.md),
[`plot_feature_expression()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_feature_expression.md),
[`plot_heatmap()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_heatmap.md),
[`plot_ma()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_ma.md),
[`plot_pca()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_pca.md),
[`plot_volcano()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_volcano.md),
[`run_diff()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_diff.md),
[`run_diff_continuous()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_diff_continuous.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  applicable_diff_methods(rnaseq_counts_input)  # deseq2, edger, ...
  applicable_diff_methods(proteomics_input)     # limma, ttest, lm
} # }
```
