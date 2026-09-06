# Limma multi-group ANOVA-style test

Builds the standardized diff schema in place (no shared standardizer) so
that F-statistic and AveExpr fields are passed through cleanly. Internal
— call via
[`run_diff()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_diff.md)
with `analysis_type = "anova"`.

## Usage

``` r
run_limma_anova(
  input,
  group_col,
  covariates = NULL,
  selected_groups = NULL,
  paired_col = NULL
)
```

## Arguments

- input:

  A validated `omics_input`.

- group_col:

  Grouping column in sample metadata.

- covariates:

  Optional covariate column names.

- selected_groups:

  Optional subset of groups to retain.

- paired_col:

  Optional pairing/block column.

## Value

List with `results_raw`, `results_std`, `model_object`, and
`analysis_info`.
