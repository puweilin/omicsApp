# Refresh the on-disk gene-set cache

Rebuilds the qs2 gene-set cache that
[`run_enrichment()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_enrichment.md),
[`run_gsva()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_gsva.md),
and
[`list_gene_sets()`](https://puweilin.github.io/omicsApp/omicsCore/reference/list_gene_sets.md)
read. For `"kegg"` the table is fetched live from the KEGG REST API
(current pathway definitions, in symbol space – about twice the pathways
of the bundled 2011 `KEGG_LEGACY` snapshot); all other databases are
re-snapshotted from the installed `msigdbr`. Every table is stamped with
a `gs_source` so downstream results can name the pathway definitions
they used.

## Usage

``` r
refresh_geneset_cache(
  databases = SUPPORTED_ENRICH_DATABASES,
  organism = "Hs",
  max_age_days = 30,
  force = FALSE,
  quiet = FALSE
)
```

## Arguments

- databases:

  Database keys to refresh. Defaults to all supported databases; only
  `"kegg"` involves the network.

- organism:

  Organism shorthand (e.g. `"Hs"`, `"Mm"`).

- max_age_days:

  Skip databases whose cache file is younger than this.

- force:

  If `TRUE`, refresh regardless of age.

- quiet:

  If `TRUE`, suppress progress messages.

## Value

Invisibly, a `data.frame` with one row per database: `database`,
`action` (`"refreshed"`, `"fresh"`, or `"failed"`), `n_sets`, `path`.

## Details

The cache directory is `OMICSCORE_GENESET_CACHE` when set, otherwise a
per-user directory (`tools::R_user_dir("omicsCore", "cache")`). Once a
live KEGG table is in place it is kept current automatically: an
enrichment run that finds it older than `OMICSCORE_GENESET_TTL_DAYS`
(default 30) re-fetches it first, falling back to the stale table if the
network is unavailable. Set `OMICSCORE_GENESET_TTL_DAYS=0` to disable.

KEGG-derived caches are for local use only; KEGG's license does not
permit redistributing them.

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
[`run_enrichment()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_enrichment.md),
[`run_gsva()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_gsva.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  refresh_geneset_cache("kegg")
  geneset_cache_status()
} # }
```
