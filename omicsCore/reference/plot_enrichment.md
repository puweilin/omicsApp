# Enrichment summary plot

Visualises the standardized enrichment table produced by
[`run_enrichment()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_enrichment.md).
The `"dot"` view draws a dotplot of the top `top_n` pathways per
database (size = overlap, color = adjusted p-value); the `"bar"` view
draws a horizontal bar chart of the same selection. For GSEA bundles a
`"gsea_dot"` view is available that splits pathways by direction (up /
down).

## Usage

``` r
plot_enrichment(
  bundle,
  top_n = 20L,
  view = c("dot", "bar", "gsea_dot"),
  p_preference = c("adjusted", "raw", "qvalue"),
  p_cutoff = NULL,
  database = NULL
)
```

## Arguments

- bundle:

  An
  [`analysis_bundle`](https://puweilin.github.io/omicsApp/omicsCore/reference/is_analysis_bundle.md)
  produced by
  [`run_enrichment()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_enrichment.md).

- top_n:

  Number of pathways to display per database.

- view:

  One of `"dot"` (default), `"bar"`, or `"gsea_dot"`.

- p_preference:

  `"adjusted"` (default), `"raw"`, or `"qvalue"`.

- p_cutoff:

  Optional significance cutoff. If `NULL`, all rows are shown (subject
  to `top_n`).

- database:

  Optional vector of databases to restrict to.

## Value

A `ggplot` object.

## See also

Other enrich:
[`filter_enrich_results()`](https://puweilin.github.io/omicsApp/omicsCore/reference/filter_enrich_results.md),
[`geneset_cache_status()`](https://puweilin.github.io/omicsApp/omicsCore/reference/geneset_cache_status.md),
[`hgnc_map_provenance()`](https://puweilin.github.io/omicsApp/omicsCore/reference/hgnc_map_provenance.md),
[`list_gene_sets()`](https://puweilin.github.io/omicsApp/omicsCore/reference/list_gene_sets.md),
[`map_ensembl_symbols()`](https://puweilin.github.io/omicsApp/omicsCore/reference/map_ensembl_symbols.md),
[`plot_gsea()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_gsea.md),
[`plot_gsva_heatmap()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_gsva_heatmap.md),
[`refresh_geneset_cache()`](https://puweilin.github.io/omicsApp/omicsCore/reference/refresh_geneset_cache.md),
[`run_enrichment()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_enrichment.md),
[`run_gsva()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_gsva.md)
