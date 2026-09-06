# List the gene sets available in a pathway database

Returns a `tibble` with one row per gene-set member, sourced from MSigDB
(`msigdbr`). All gene identifiers are returned as HGNC / official gene
symbols. Results are cached for the duration of the R session so
repeated calls are cheap.

## Usage

``` r
list_gene_sets(database, organism = "Hs")
```

## Arguments

- database:

  One of the supported database keys (see above).

- organism:

  Organism, accepting `"Hs"`, `"human"`, `"Homo sapiens"`, `"Mm"`,
  `"mouse"`, or `"Mus musculus"`.

## Value

A `tibble` with columns `database`, `pathway_id`, `pathway_name`,
`gene_symbol`.

## Details

Supported databases (lowercase): `hallmark`, `kegg`, `reactome`,
`wikipathways`, `go_bp`, `go_mf`, `go_cc`.

## See also

Other enrich:
[`filter_enrich_results()`](https://puweilin.github.io/omicsApp/omicsCore/reference/filter_enrich_results.md),
[`geneset_cache_status()`](https://puweilin.github.io/omicsApp/omicsCore/reference/geneset_cache_status.md),
[`hgnc_map_provenance()`](https://puweilin.github.io/omicsApp/omicsCore/reference/hgnc_map_provenance.md),
[`map_ensembl_symbols()`](https://puweilin.github.io/omicsApp/omicsCore/reference/map_ensembl_symbols.md),
[`plot_enrichment()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_enrichment.md),
[`plot_gsea()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_gsea.md),
[`plot_gsva_heatmap()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_gsva_heatmap.md),
[`refresh_geneset_cache()`](https://puweilin.github.io/omicsApp/omicsCore/reference/refresh_geneset_cache.md),
[`run_enrichment()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_enrichment.md),
[`run_gsva()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_gsva.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  list_gene_sets("hallmark")
  list_gene_sets("go_bp", organism = "Hs")
} # }
```
