# edgeR QL F-test for a two-group comparison

edgeR QL F-test for a two-group comparison

## Usage

``` r
run_edger_group(
  input,
  group_col,
  control_group,
  case_group,
  covariates = NULL,
  paired_col = NULL
)
```

## Arguments

- input:

  A validated RNA-seq `omics_input` built from raw counts.

- group_col:

  Group column in sample metadata.

- control_group:

  Control-group label.

- case_group:

  Case-group label.

- covariates:

  Optional covariate column names.

- paired_col:

  Optional pairing column.

## Value

List with `results_raw`, `results_std`, `model_object` (`DGEGLM`), and
`analysis_info`.
