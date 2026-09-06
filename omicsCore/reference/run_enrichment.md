# Run pathway enrichment from a differential bundle

Single entry point for over-representation analysis (`type = "ora"`) and
rank-based gene-set enrichment analysis (`type = "gsea"`). Both backends
run against the requested MSigDB database in symbol space via
[`clusterProfiler::enricher()`](https://rdrr.io/pkg/clusterProfiler/man/enricher.html)
/ `GSEA()`, so no `org.*` annotation package is required.

## Usage

``` r
run_enrichment(
  diff_bundle,
  type = c("ora", "gsea"),
  database = c("hallmark", "kegg", "reactome", "go_bp", "go_mf", "go_cc", "wikipathways"),
  organism = "Hs",
  direction = c("both", "up", "down"),
  p_cutoff = 0.05,
  output_p_cutoff = NULL,
  p_preference = c("adjusted", "raw"),
  effect_cutoff = NULL,
  p_adjust_method = "BH",
  min_size = 10L,
  max_size = 500L,
  ...
)
```

## Arguments

- diff_bundle:

  An
  [`analysis_bundle`](https://puweilin.github.io/omicsApp/omicsCore/reference/is_analysis_bundle.md)
  produced by
  [`run_diff()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_diff.md).

- type:

  One of `"ora"` or `"gsea"`.

- database:

  One of `"hallmark"`, `"kegg"`, `"reactome"`, `"wikipathways"`,
  `"go_bp"`, `"go_mf"`, `"go_cc"`. Pass a character vector to query
  multiple databases.

- organism:

  Organism shorthand (e.g. `"Hs"`).

- direction:

  One of `"both"` (default), `"up"`, or `"down"`.

- p_cutoff:

  Significance cutoff for selecting diff features (ORA only; GSEA ranks
  the whole list).

- output_p_cutoff:

  Bound on the returned table. Defaults to `p_cutoff`. Pass `1` to keep
  every pathway, so a caller can choose raw or adjusted p at display
  time without re-running.

- p_preference:

  For ORA feature selection: `"adjusted"` (default) or `"raw"`.

- effect_cutoff:

  Optional \|effect\| cutoff for ORA feature selection.

- p_adjust_method:

  Multiple-testing correction method.

- min_size, max_size:

  GSEA min/max gene-set sizes.

- ...:

  Reserved for backend-specific extensions.

## Value

An
[`analysis_bundle`](https://puweilin.github.io/omicsApp/omicsCore/reference/is_analysis_bundle.md)
with `results$enrich_result_df` (standardized schema) and
`results$enrich_object` (named list of clusterProfiler objects keyed by
database, or for ORA by `<direction>__<database>`).

## Details

For ORA, features are split into up/down sets unless
`direction = "both"`; for GSEA, the full ranked vector is always passed
to the backend and `direction` filters the standardized output.

For GSVA-style sample-level scoring, see
[`run_gsva()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_gsva.md).

## See also

Other enrich:
[`filter_enrich_results()`](https://puweilin.github.io/omicsApp/omicsCore/reference/filter_enrich_results.md),
[`geneset_cache_status()`](https://puweilin.github.io/omicsApp/omicsCore/reference/geneset_cache_status.md),
[`hgnc_map_provenance()`](https://puweilin.github.io/omicsApp/omicsCore/reference/hgnc_map_provenance.md),
[`list_gene_sets()`](https://puweilin.github.io/omicsApp/omicsCore/reference/list_gene_sets.md),
[`map_ensembl_symbols()`](https://puweilin.github.io/omicsApp/omicsCore/reference/map_ensembl_symbols.md),
[`plot_enrichment()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_enrichment.md),
[`plot_gsea()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_gsea.md),
[`plot_gsva_heatmap()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_gsva_heatmap.md),
[`refresh_geneset_cache()`](https://puweilin.github.io/omicsApp/omicsCore/reference/refresh_geneset_cache.md),
[`run_gsva()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_gsva.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  diff <- run_diff(input, method = "limma", analysis_type = "group",
                   group_col = "treatment",
                   control_group = "DMSO", case_group = "Drug")
  enr <- run_enrichment(diff, type = "ora", database = "hallmark")
  head(enr$results$enrich_result_df)
} # }
```
