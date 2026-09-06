# Filter standardized enrichment results

Subsets a standardized enrichment `data.frame` (the schema returned by
[`run_enrichment()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_enrichment.md))
to pathways that pass a significance cutoff on either the adjusted, raw,
or q-value column, with an optional minimum gene-set size and direction
filter.

## Usage

``` r
filter_enrich_results(
  enrich_df,
  p_cutoff = 0.05,
  p_preference = c("adjusted", "raw", "qvalue"),
  min_genes = NULL,
  direction = NULL
)
```

## Arguments

- enrich_df:

  Standardized enrichment `data.frame`.

- p_cutoff:

  Significance cutoff.

- p_preference:

  One of `"adjusted"` (default), `"raw"`, or `"qvalue"`.

- min_genes:

  Optional minimum number of overlapping / leading genes.

- direction:

  Optional direction filter (`"up"` or `"down"`). Useful for GSEA where
  pathways carry a sign.

## Value

Filtered standardized enrichment `data.frame`.

## See also

Other enrich:
[`geneset_cache_status()`](https://puweilin.github.io/omicsApp/omicsCore/reference/geneset_cache_status.md),
[`hgnc_map_provenance()`](https://puweilin.github.io/omicsApp/omicsCore/reference/hgnc_map_provenance.md),
[`list_gene_sets()`](https://puweilin.github.io/omicsApp/omicsCore/reference/list_gene_sets.md),
[`map_ensembl_symbols()`](https://puweilin.github.io/omicsApp/omicsCore/reference/map_ensembl_symbols.md),
[`plot_enrichment()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_enrichment.md),
[`plot_gsea()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_gsea.md),
[`plot_gsva_heatmap()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_gsva_heatmap.md),
[`refresh_geneset_cache()`](https://puweilin.github.io/omicsApp/omicsCore/reference/refresh_geneset_cache.md),
[`run_enrichment()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_enrichment.md),
[`run_gsva()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_gsva.md)
