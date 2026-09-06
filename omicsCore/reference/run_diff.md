# Run a differential-expression analysis

Single entry point for proteomics and RNA-seq differential analysis.
Dispatches to the appropriate backend (DESeq2, edgeR, limma, t-test, or
linear model) and returns an
[`analysis_bundle`](https://puweilin.github.io/omicsApp/omicsCore/reference/is_analysis_bundle.md)
wrapping the standardized result, the backend's raw table, and the
fitted model object.

## Usage

``` r
run_diff(
  input,
  method = "auto",
  analysis_type = c("group", "continuous", "anova"),
  group_col = NULL,
  control_group = NULL,
  case_group = NULL,
  continuous_col = NULL,
  covariates = NULL,
  paired_col = NULL,
  selected_groups = NULL,
  ...
)
```

## Arguments

- input:

  A validated `omics_input`.

- method:

  Backend name. `"auto"` (default) lets `omicsCore` pick one based on
  `omics_type` and installed Suggests.

- analysis_type:

  One of `"group"`, `"continuous"`, or `"anova"`.

- group_col:

  Group column in sample metadata (group/anova).

- control_group:

  Control-group label (group only).

- case_group:

  Case-group label (group only).

- continuous_col:

  Continuous metadata column (continuous only).

- covariates:

  Optional character vector of covariate column names.

- paired_col:

  Optional pairing/block column.

- selected_groups:

  Optional subset of groups to retain (anova only).

- ...:

  Extra arguments forwarded to the backend, e.g. `var_equal` for t-test
  or `df = 3` for limma spline.

## Value

An
[`analysis_bundle`](https://puweilin.github.io/omicsApp/omicsCore/reference/is_analysis_bundle.md)
with `results$diff_result_df` (standardized schema),
`results$diff_raw_df` (backend-native), and `results$diff_object`
(fitted model, may be `NULL`).

## Details

`method = "auto"` picks limma for proteomics and DESeq2 for raw-count
RNA-seq. When the preferred Bioconductor backend is not installed, it
silently falls back to t-test / lm (with a message) so analyses still
run in restricted environments without
[`omicsCore::install_optional()`](https://puweilin.github.io/omicsApp/omicsCore/reference/install_optional.md)
being invoked first.

For paired designs supply `paired_col`; for ANOVA-style multi-group
tests set `analysis_type = "anova"` (currently limma-backed only). The
`continuous` analysis type requires `continuous_col` instead of
`group_col` + `control_group` + `case_group`.

## See also

Other diff:
[`SUPPORTED_DIFF_METHODS`](https://puweilin.github.io/omicsApp/omicsCore/reference/SUPPORTED_DIFF_METHODS.md),
[`applicable_diff_methods()`](https://puweilin.github.io/omicsApp/omicsCore/reference/applicable_diff_methods.md),
[`filter_diff_results()`](https://puweilin.github.io/omicsApp/omicsCore/reference/filter_diff_results.md),
[`make_ranked_features()`](https://puweilin.github.io/omicsApp/omicsCore/reference/make_ranked_features.md),
[`plot_feature_expression()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_feature_expression.md),
[`plot_heatmap()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_heatmap.md),
[`plot_ma()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_ma.md),
[`plot_pca()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_pca.md),
[`plot_volcano()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_volcano.md),
[`run_diff_continuous()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_diff_continuous.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  res <- run_diff(input,
                  analysis_type = "group",
                  group_col = "treatment",
                  control_group = "DMSO",
                  case_group = "Drug")
  head(res$results$diff_result_df)
} # }
```
