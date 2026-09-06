# Report the state of the on-disk gene-set cache

One row per database: whether a cache file is present and readable, how
many gene sets it holds, where it came from (`gs_source`, e.g.
`"KEGG REST 2026-08-19"`), and how old it is. Useful as a pre-flight
check before deciding whether
[`refresh_geneset_cache()`](https://puweilin.github.io/omicsApp/omicsCore/reference/refresh_geneset_cache.md)
is worth a network round-trip.

## Usage

``` r
geneset_cache_status(databases = SUPPORTED_ENRICH_DATABASES, organism = "Hs")
```

## Arguments

- databases:

  Database keys to report on. Defaults to all supported.

- organism:

  Organism shorthand (e.g. `"Hs"`, `"Mm"`).

## Value

A `data.frame` with columns `database`, `cached`, `n_sets`, `source`,
`cached_at`, `age_days`.

## See also

Other enrich:
[`filter_enrich_results()`](https://puweilin.github.io/omicsApp/omicsCore/reference/filter_enrich_results.md),
[`hgnc_map_provenance()`](https://puweilin.github.io/omicsApp/omicsCore/reference/hgnc_map_provenance.md),
[`list_gene_sets()`](https://puweilin.github.io/omicsApp/omicsCore/reference/list_gene_sets.md),
[`map_ensembl_symbols()`](https://puweilin.github.io/omicsApp/omicsCore/reference/map_ensembl_symbols.md),
[`plot_enrichment()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_enrichment.md),
[`plot_gsea()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_gsea.md),
[`plot_gsva_heatmap()`](https://puweilin.github.io/omicsApp/omicsCore/reference/plot_gsva_heatmap.md),
[`refresh_geneset_cache()`](https://puweilin.github.io/omicsApp/omicsCore/reference/refresh_geneset_cache.md),
[`run_enrichment()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_enrichment.md),
[`run_gsva()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_gsva.md)
