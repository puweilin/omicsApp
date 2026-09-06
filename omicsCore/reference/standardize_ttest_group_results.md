# Standardize per-feature t-test group results

Standardize per-feature t-test group results

## Usage

``` r
standardize_ttest_group_results(
  raw_df,
  feature_df,
  comparison,
  omics_type = "proteomics"
)
```

## Arguments

- raw_df:

  Per-feature t-test table with `feature_id`, `mean_diff`, `t_stat`,
  `p_value`, `adj_p_value`, optionally `mean_ctrl`, `mean_case`.

- feature_df:

  Feature metadata.

- comparison:

  Comparison label.

- omics_type:

  Omics modality label.

## Value

Standardized diff result data frame.
