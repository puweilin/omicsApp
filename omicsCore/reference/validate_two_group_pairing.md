# Validate two-group paired design completeness

Ensures every pair has exactly one sample in each of the two target
groups.

## Usage

``` r
validate_two_group_pairing(
  meta_df,
  group_col,
  paired_col,
  control_group,
  case_group,
  object_name = "meta_df"
)
```

## Arguments

- meta_df:

  Sample metadata already filtered to target samples.

- group_col:

  Group column name.

- paired_col:

  Pairing column name.

- control_group:

  Control-group label.

- case_group:

  Case-group label.

- object_name:

  Object label used in error messages.

## Value

Invisibly returns `TRUE` on success.
