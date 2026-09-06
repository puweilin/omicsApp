# Add an experiment to a project

Add an experiment to a project

## Usage

``` r
add_experiment(project, name, input)
```

## Arguments

- project:

  An `omics_project`.

- name:

  Experiment tag (e.g. `"proteomics"`, `"rnaseq"`). Must be unique
  within the project.

- input:

  An `omics_input` to attach.

## Value

The modified `omics_project`.

## See also

Other omics_project:
[`experiment_tags()`](https://puweilin.github.io/omicsApp/omicsCore/reference/experiment_tags.md),
[`is_omics_project()`](https://puweilin.github.io/omicsApp/omicsCore/reference/is_omics_project.md),
[`omics_project()`](https://puweilin.github.io/omicsApp/omicsCore/reference/omics_project.md),
[`remove_experiment()`](https://puweilin.github.io/omicsApp/omicsCore/reference/remove_experiment.md)
