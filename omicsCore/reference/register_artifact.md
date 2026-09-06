# Register an analysis artifact

Register an analysis artifact

## Usage

``` r
register_artifact(registry, artifact_type, label, path = NA_character_)
```

## Arguments

- registry:

  Existing registry or `NULL`.

- artifact_type:

  Artifact category such as `"plot"` or `"table"`.

- label:

  Logical artifact label.

- path:

  Optional file path.

## Value

Updated registry.
