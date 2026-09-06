# Standardize limma two-group results

Standardize limma two-group results

## Usage

``` r
standardize_limma_group_results(
  raw_df,
  feature_df,
  comparison,
  omics_type = "proteomics"
)
```

## Arguments

- raw_df:

  Raw limma topTable result with `feature_id`, `logFC`, `P.Value`,
  `adj.P.Val`.

- feature_df:

  Feature metadata.

- comparison:

  Comparison label.

- omics_type:

  Omics modality label.

## Value

Standardized diff result data frame.
