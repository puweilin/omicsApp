# Run Gene Set Variation Analysis

Computes sample-level pathway activity scores via the `GSVA` package.
The caller may either pass a named list of gene sets directly
(`gene_sets`) or pick one of the supported MSigDB databases
(`database`). RNA-seq raw counts are log2-transformed before scoring.

## Usage

``` r
run_gsva(
  input,
  gene_sets = NULL,
  database = "hallmark",
  organism = "Hs",
  method = c("gsva", "ssgsea"),
  min_size = 10L,
  max_size = 500L,
  kcdf = NULL,
  ...
)
```

## Arguments

- input:

  A validated `omics_input`.

- gene_sets:

  Optional named list of character vectors. When supplied, `database` is
  ignored.

- database:

  Database key when `gene_sets` is `NULL`. See
  [`list_gene_sets()`](https://puweilin.github.io/omicsApp/omicsCore/reference/list_gene_sets.md)
  for supported keys.

- organism:

  Organism shorthand (e.g. `"Hs"`).

- method:

  One of `"gsva"` (default) or `"ssgsea"`.

- min_size, max_size:

  Min / max gene-set size after intersecting with the assay's feature
  set.

- kcdf:

  Kernel CDF used by
  [`GSVA::gsvaParam()`](https://rdrr.io/pkg/GSVA/man/gsvaParam-class.html)
  — typically `"Gaussian"` for log-scale data, `"Poisson"` for raw
  counts.

- ...:

  Reserved for future arguments.

## Value

An
[`analysis_bundle`](https://puweilin.github.io/omicsApp/omicsCore/reference/is_analysis_bundle.md)
with `results$gsva_matrix` (pathways × samples) and
`results$gsva_gene_sets` (the gene-set list actually used).

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
[`run_enrichment()`](https://puweilin.github.io/omicsApp/omicsCore/reference/run_enrichment.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  gsva_bundle <- run_gsva(input, database = "hallmark")
  head(gsva_bundle$results$gsva_matrix[, 1:4])
} # }
```
