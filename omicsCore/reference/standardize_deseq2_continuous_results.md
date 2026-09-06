# Standardize DESeq2 continuous results

Standardize DESeq2 continuous results

## Usage

``` r
standardize_deseq2_continuous_results(
  raw_df,
  feature_df,
  comparison,
  omics_type = "rnaseq"
)
```

## Arguments

- raw_df:

  Raw DESeq2 result with `feature_id`, `log2FoldChange`, `pvalue`,
  `padj`.

- feature_df:

  Feature metadata.

- comparison:

  Continuous variable name.

- omics_type:

  Omics modality label.

## Value

Standardized diff result data frame.
