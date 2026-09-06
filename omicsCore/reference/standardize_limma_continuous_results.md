# Standardize limma continuous results

Standardize limma continuous results

## Usage

``` r
standardize_limma_continuous_results(
  raw_df,
  feature_df,
  comparison,
  omics_type = "proteomics"
)
```

## Arguments

- raw_df:

  Raw limma topTable result with `feature_id`, `P.Value`, `adj.P.Val`,
  `spearman_rho`, `adj_r_squared`.

- feature_df:

  Feature metadata.

- comparison:

  Continuous variable name.

- omics_type:

  Omics modality label.

## Value

Standardized diff result data frame.
