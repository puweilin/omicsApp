# Run a cross-omics integration analysis

Single entry point for dual-omics integration. Three methods are
supported:

## Usage

``` r
run_integration(
  project,
  method = c("correlation", "concordance", "active_pathways"),
  experiments = NULL,
  diff_bundles = NULL,
  by = "feature_symbol",
  ...
)
```

## Arguments

- project:

  An
  [`omics_project`](https://puweilin.github.io/omicsApp/omicsCore/reference/is_omics_project.md)
  containing the two experiments to integrate.

- method:

  One of `"correlation"`, `"concordance"`, `"active_pathways"`.

- experiments:

  Length-2 character vector of experiment tags. If `NULL` and the
  project has exactly two layers, both are used.

- diff_bundles:

  Named list of
  [`run_diff()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_diff.md)
  bundles keyed by experiment tag. Required for `"concordance"` and
  `"active_pathways"`.

- by:

  Feature key used to join layers (default `"feature_symbol"`).

- ...:

  Method-specific arguments forwarded to the backend. See the per-method
  sections below.

## Value

An
[`analysis_bundle`](https://puweilin.github.io/omicsApp/omicsCore/reference/is_analysis_bundle.md)
with `results$integration_df` (standardized schema) and, for
`"active_pathways"`, `results$integration_raw` carrying the raw
`ActivePathways` table.

## Details

- `"correlation"` – per-feature Pearson or Spearman correlation across
  paired samples from two experiments. Useful for RNA / protein layers
  that share donors.

- `"concordance"` – per-feature agreement between two differential
  analyses, classified into the four sign quadrants and combined via
  Fisher's method.

- `"active_pathways"` – pathway-level combined-p enrichment via the
  `ActivePathways` package, fed by the p-values of two differential
  analyses.

All methods return an `analysis_bundle` whose `results$integration_df`
follows the schema documented in
[`check_integration_result_schema()`](https://puweilin.github.io/omicsApp/omicsCore/reference/check_integration_result_schema.md).

## Correlation arguments

- `cor_method` – `"spearman"` (default) or `"pearson"`.

- `p_adjust_method` – defaults to `"BH"`.

- `min_samples` – minimum paired samples (default `4`).

- `p_cutoff` – significance cutoff (default `0.05`).

## Concordance arguments

- `p_preference` – `"adjusted"` (default) or `"raw"`.

- `p_cutoff` – significance cutoff (default `0.05`).

- `p_adjust_method` – defaults to `"BH"`.

## ActivePathways arguments

- `database` – MSigDB shorthand (default `"hallmark"`).

- `organism` – defaults to `"Hs"`.

- `p_preference` – `"raw"` (default) or `"adjusted"`. ActivePathways
  prefers raw p-values since it applies its own multiple-testing
  correction across pathways.

- `significant` – pathway-level cutoff (default `0.05`).

- `geneset_filter` – length-2 integer vector of min/max gene-set sizes
  (default `c(5L, 1000L)`).

- `merge_method` – `"Brown"` (default) or another method accepted by
  [`ActivePathways::ActivePathways()`](https://rdrr.io/pkg/ActivePathways/man/ActivePathways.html).

## See also

Other integration:
[`plot_integration()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_integration.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  # Correlation across paired RNA / protein samples
  cor_bundle <- run_integration(project, method = "correlation",
                                experiments = c("rna", "prot"))

  # Concordance from two diff bundles
  diff_rna <- run_diff(project$experiments$rna, ...)
  diff_prot <- run_diff(project$experiments$prot, ...)
  con_bundle <- run_integration(
    project, method = "concordance",
    experiments = c("rna", "prot"),
    diff_bundles = list(rna = diff_rna, prot = diff_prot)
  )
} # }
```
