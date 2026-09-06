# GSVA heatmap

Draws a pathways-by-samples heatmap from a
[`run_gsva()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_gsva.md)
bundle. When `ComplexHeatmap` is installed the function returns a
`Heatmap` object (recommended); otherwise it falls back to a ggplot2
tile plot so the function is always usable from a clean install.

## Usage

``` r
plot_gsva_heatmap(
  bundle,
  top_n = 25L,
  pathways = NULL,
  meta_df = NULL,
  annotation_cols = NULL,
  scale = c("row", "none", "column"),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  title = "GSVA"
)
```

## Arguments

- bundle:

  An
  [`analysis_bundle`](https://puweilin.github.io/omicsApp/omicsCore/reference/is_analysis_bundle.md)
  produced by
  [`run_gsva()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_gsva.md).

- top_n:

  Number of pathways to plot when `pathways` is `NULL`. Selected by
  variance across samples.

- pathways:

  Optional character vector of pathway names to plot.

- meta_df:

  Optional sample metadata `data.frame` (or named list of columns) used
  for column annotations. Rows must be indexable by sample id (column
  names of the score matrix).

- annotation_cols:

  Optional character vector of columns from `meta_df` to surface as a
  `HeatmapAnnotation`.

- scale:

  One of `"row"` (default), `"none"`, or `"column"`.

- cluster_rows, cluster_cols:

  Whether to cluster rows / columns.

- title:

  Heatmap title.

## Value

A
[`ComplexHeatmap::Heatmap`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
object when `ComplexHeatmap` is available, otherwise a `ggplot` tile
plot.

## See also

Other enrich:
[`filter_enrich_results()`](https://puweilin.github.io/omicsApp/omicsCore/reference/filter_enrich_results.md),
[`geneset_cache_status()`](https://puweilin.github.io/omicsApp/omicsCore/reference/geneset_cache_status.md),
[`hgnc_map_provenance()`](https://puweilin.github.io/omicsApp/omicsCore/reference/hgnc_map_provenance.md),
[`list_gene_sets()`](https://puweilin.github.io/omicsApp/omicsCore/reference/list_gene_sets.md),
[`map_ensembl_symbols()`](https://puweilin.github.io/omicsApp/omicsCore/reference/map_ensembl_symbols.md),
[`plot_enrichment()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_enrichment.md),
[`plot_gsea()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_gsea.md),
[`refresh_geneset_cache()`](https://puweilin.github.io/omicsApp/omicsCore/reference/refresh_geneset_cache.md),
[`run_enrichment()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_enrichment.md),
[`run_gsva()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_gsva.md)
