# GSEA running-score plot

Plots the `clusterProfiler::gseaplot2()` running-score curve for a
single pathway in a GSEA enrichment bundle. Requires the `enrichplot`
Bioconductor package.

## Usage

``` r
plot_gsea(bundle, pathway_id, database = NULL)
```

## Arguments

- bundle:

  An
  [`analysis_bundle`](https://puweilin.github.io/omicsApp/omicsCore/reference/is_analysis_bundle.md)
  produced by
  [`run_enrichment()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_enrichment.md)
  with `type = "gsea"`.

- pathway_id:

  Pathway ID (matches `pathway_id` in the standardized table) or pathway
  name.

- database:

  Database key. Required when the bundle holds multiple databases;
  ignored when there's only one.

## Value

A `ggplot` object.

## See also

Other enrich:
[`filter_enrich_results()`](https://puweilin.github.io/omicsApp/omicsCore/reference/filter_enrich_results.md),
[`geneset_cache_status()`](https://puweilin.github.io/omicsApp/omicsCore/reference/geneset_cache_status.md),
[`hgnc_map_provenance()`](https://puweilin.github.io/omicsApp/omicsCore/reference/hgnc_map_provenance.md),
[`list_gene_sets()`](https://puweilin.github.io/omicsApp/omicsCore/reference/list_gene_sets.md),
[`map_ensembl_symbols()`](https://puweilin.github.io/omicsApp/omicsCore/reference/map_ensembl_symbols.md),
[`plot_enrichment()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_enrichment.md),
[`plot_gsva_heatmap()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_gsva_heatmap.md),
[`refresh_geneset_cache()`](https://puweilin.github.io/omicsApp/omicsCore/reference/refresh_geneset_cache.md),
[`run_enrichment()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_enrichment.md),
[`run_gsva()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_gsva.md)
