# Validate continuous paired design minimum replication

Ensures every pair has at least two observations.

## Usage

``` r
validate_continuous_pairing(meta_df, paired_col, object_name = "meta_df")
```

## Arguments

- meta_df:

  Sample metadata.

- paired_col:

  Pairing column name.

- object_name:

  Object label used in error messages.

## Value

Invisibly returns `TRUE` on success.
