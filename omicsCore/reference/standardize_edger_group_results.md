# Standardize edgeR two-group results

Standardize edgeR two-group results

## Usage

``` r
standardize_edger_group_results(
  raw_df,
  feature_df,
  comparison,
  omics_type = "rnaseq"
)
```

## Arguments

- raw_df:

  Raw edgeR topTags result with `feature_id`, `logFC`, `PValue`, `FDR`.

- feature_df:

  Feature metadata.

- comparison:

  Comparison label.

- omics_type:

  Omics modality label.

## Value

Standardized diff result data frame.
