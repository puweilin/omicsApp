# Map Ensembl gene ids to HGNC symbols

Map Ensembl gene ids to HGNC symbols

## Usage

``` r
map_ensembl_symbols(ids)
```

## Arguments

- ids:

  Character vector of Ensembl gene ids, with or without the version
  suffix.

## Value

Character vector the same length as `ids`; `NA` where the id has no
approved HGNC symbol.

## See also

Other enrich:
[`filter_enrich_results()`](https://puweilin.github.io/omicsApp/omicsCore/reference/filter_enrich_results.md),
[`geneset_cache_status()`](https://puweilin.github.io/omicsApp/omicsCore/reference/geneset_cache_status.md),
[`hgnc_map_provenance()`](https://puweilin.github.io/omicsApp/omicsCore/reference/hgnc_map_provenance.md),
[`list_gene_sets()`](https://puweilin.github.io/omicsApp/omicsCore/reference/list_gene_sets.md),
[`plot_enrichment()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_enrichment.md),
[`plot_gsea()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_gsea.md),
[`plot_gsva_heatmap()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_gsva_heatmap.md),
[`refresh_geneset_cache()`](https://puweilin.github.io/omicsApp/omicsCore/reference/refresh_geneset_cache.md),
[`run_enrichment()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_enrichment.md),
[`run_gsva()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_gsva.md)

## Examples

``` r
map_ensembl_symbols(c("ENSG00000141510", "ENSG00000012048"))
#> [1] "TP53"  "BRCA1"
```
