# Create an empty artifact registry

An artifact registry is a `data.frame` tracking files written to disk by
analysis functions. Each row records one artifact's type, label, file
path, and creation timestamp.

## Usage

``` r
new_artifact_registry()
```

## Value

An empty registry `data.frame`.
