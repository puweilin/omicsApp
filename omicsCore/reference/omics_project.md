# Multi-omics project container

An `omics_project` groups multiple
[`omics_input()`](https://puweilin.github.io/omicsApp/omicsCore/reference/omics_input.md)
experiments together with optional cross-omics mapping tables, so that
integration analyses (RNA × Protein, ActivePathways, etc.) can operate
on a single root object.

## Usage

``` r
omics_project(
  name,
  experiments = list(),
  sample_link = NULL,
  feature_link = NULL,
  metadata = list()
)
```

## Arguments

- name:

  Human-readable project label.

- experiments:

  Named list of
  [`omics_input()`](https://puweilin.github.io/omicsApp/omicsCore/reference/omics_input.md)
  objects. Names become experiment tags. May be empty; use
  [`add_experiment()`](https://puweilin.github.io/omicsApp/omicsCore/reference/add_experiment.md)
  to populate.

- sample_link:

  Optional `data.frame` with at least the columns `tag`, `sample_id`,
  `donor_id`. `tag` matches a name in `experiments`. Each row links a
  per-experiment sample to a shared donor.

- feature_link:

  Optional `data.frame` mapping feature IDs across omics (e.g. columns
  `uniprot`, `gene_symbol`).

- metadata:

  Optional named list of arbitrary project-level metadata.

## Value

An object of class `omics_project`.

## Details

Each entry in `experiments` is one omics layer (e.g. `"proteomics"`,
`"rnaseq"`), indexed by an arbitrary tag. `sample_link` maps sample IDs
across layers so paired samples can be aligned; `feature_link` maps
feature IDs (e.g. UniProt accession ↔ gene symbol) so cross-omics
integration can join on shared identifiers.

## See also

Other omics_project:
[`add_experiment()`](https://puweilin.github.io/omicsApp/omicsCore/reference/add_experiment.md),
[`experiment_tags()`](https://puweilin.github.io/omicsApp/omicsCore/reference/experiment_tags.md),
[`is_omics_project()`](https://puweilin.github.io/omicsApp/omicsCore/reference/is_omics_project.md),
[`remove_experiment()`](https://puweilin.github.io/omicsApp/omicsCore/reference/remove_experiment.md)

## Examples

``` r
p <- omics_project(name = "demo")
p
#> <omics_project>
#>   name        : demo 
#>   experiments : 0 
#>   sample_link : <none> 
#>   feature_link: <none> 
```
