# Per-feature two-group t-test

Welch (default) or paired t-test per feature, comparing `case_group`
against `control_group`. RNA-seq raw counts are automatically log2(x+1)
transformed before testing. Internal — call via
[`run_diff()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_diff.md).

## Usage

``` r
run_ttest_group(
  input,
  group_col,
  control_group,
  case_group,
  var_equal = FALSE,
  paired_col = NULL
)
```

## Arguments

- input:

  A validated `omics_input`.

- group_col:

  Group column in sample metadata.

- control_group:

  Control-group label.

- case_group:

  Case-group label.

- var_equal:

  If `TRUE`, use equal-variance t-test. Default `FALSE` (Welch).

- paired_col:

  Optional pairing column for paired t-test.

## Value

List with `results_raw`, `results_std`, `model_object` (`NULL`), and
`analysis_info`.
