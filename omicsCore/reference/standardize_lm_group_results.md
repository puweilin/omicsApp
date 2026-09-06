# Standardize per-feature lm group results

Standardize per-feature lm group results

## Usage

``` r
standardize_lm_group_results(
  raw_df,
  feature_df,
  comparison,
  omics_type = "proteomics"
)
```

## Arguments

- raw_df:

  Per-feature lm table with `feature_id`, `beta`, `t_stat`, `p_value`,
  `adj_p_value`, `adj_r_squared`, `base_mean`.

- feature_df:

  Feature metadata.

- comparison:

  Comparison label.

- omics_type:

  Omics modality label.

## Value

Standardized diff result data frame.
