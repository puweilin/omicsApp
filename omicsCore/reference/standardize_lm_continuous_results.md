# Standardize per-feature lm continuous results

Standardize per-feature lm continuous results

## Usage

``` r
standardize_lm_continuous_results(
  raw_df,
  feature_df,
  comparison,
  omics_type = "proteomics"
)
```

## Arguments

- raw_df:

  Per-feature lm table with `feature_id`, `spearman_rho`, `t_stat`,
  `p_value`, `adj_p_value`, `adj_r_squared`, `base_mean`.

- feature_df:

  Feature metadata.

- comparison:

  Continuous variable name.

- omics_type:

  Omics modality label.

## Value

Standardized diff result data frame.
